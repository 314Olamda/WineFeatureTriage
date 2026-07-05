#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wine_feature_classifier.py — Family triage + 2-10 AA peptide scoring tool
════════════════════════════════════════════════════════════════════════

Simulates the "online tool" step: takes an uploaded, UNANNOTATED UHPLC-RP
feature table (RT, m/z, adduct, CCS, MS fragments — same columns as
generate_template.py's output) and produces a multi-sheet results
workbook:

  Sheet 1 "Summary_Families"   — how many features per chemical family
  Sheet 2 "Oligopeptides_gt50AA" — features consistent with linear
                                    peptides > 50 AA (flagged, not scored
                                    — outside this tool's 2-10 AA focus)
  Sheet 3 "Peptides_2to10AA"   — every feature consistent with a 2-10 AA
                                    linear peptide, PLUS full a/b/y/immonium
                                    de novo scoring against its own MS
                                    fragment list (composite, tier, q_storey)

IMPORTANT — classification method and its limits
─────────────────────────────────────────────────
Family triage (pesticide / lipid / polysaccharide / sugar / protein /
peptide) here uses MASS-WINDOW + CCS-TREND heuristics only. This is a
reasonable FIRST-PASS filter for a template/demo, but it is NOT a
substitute for real identification: pesticides, lipids, and small
peptides can and do overlap in mass, and CCS trend lines vary by
functional class beyond what a single linear model captures. Real
triage should cross-reference formula prediction, isotope pattern,
and spectral library matching (GNPS, MassBank, mzCloud) before
reporting a family assignment. The 2-10 AA peptide layer is more
rigorous — it uses genuine mass-directed composition search plus
a/b/y/immonium fragment scoring — but is still a de novo CALL requiring
experimental confirmation for any candidate reported as HIGH confidence.

Requirements
────────────
pip install pandas numpy scipy pyteomics openpyxl
mass_utils.py, linear_peptide_scoring.py, linear_denovo.py must be
in the same directory.

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from mass_utils import precursor_neutral, H2O, RESIDUE_MASS
from linear_denovo import identify_linear_peptide, compositions_for_linear_mass

# ══════════════════════════════════════════════════════════════════════════
# CONFIGURATION — mass windows for first-pass family triage
# ══════════════════════════════════════════════════════════════════════════

# Rough monoisotopic neutral-mass windows (Da) used ONLY for triage.
# Deliberately permissive at the edges; peptide-range candidates are always
# cross-checked against the real linear-peptide composition search before
# being called a peptide, which is the more reliable test.
FAMILY_MASS_WINDOWS = {
    "Pesticide":       (150, 800),
    "Sugar":           (150, 600),
    "Lipid":           (250, 950),
    "Polysaccharide":  (250, 1500),
    "Protein":         (5000, 200000),
}
OLIGOPEPTIDE_MIN_AA = 51
OLIGOPEPTIDE_MAX_AA = 300
PEPTIDE_MIN_AA = 2
PEPTIDE_MAX_AA = 10          # this tool's scoring focus
PEPTIDE_MAX_AA_EXTENDED = 50  # still "peptide" but outside 2-10 AA scoring scope

MIN_RESIDUE = min(RESIDUE_MASS.values())   # smallest possible residue (Gly)
MAX_RESIDUE = max(RESIDUE_MASS.values())   # largest possible residue (Trp)


def parse_fragments(frag_string: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse 'mz:intensity;mz:intensity;...' into two numpy arrays."""
    if not frag_string or not isinstance(frag_string, str) or frag_string.strip() == "":
        return np.array([]), np.array([])
    mzs, ints = [], []
    for pair in frag_string.split(";"):
        pair = pair.strip()
        if not pair or ":" not in pair:
            continue
        mz_s, int_s = pair.split(":")
        try:
            mzs.append(float(mz_s))
            ints.append(float(int_s))
        except ValueError:
            continue
    return np.array(mzs), np.array(ints)


def neutral_mass_from_row(mz: float, adduct: str, charge: int) -> float:
    """Back-calculate neutral mass from observed m/z, adduct, and charge."""
    adduct_shift = {
        "[M+H]+": 1.007276, "[M+Na]+": 22.989218, "[M+NH4]+": 18.033823,
    }.get(adduct, 1.007276)
    return mz * charge - charge * adduct_shift


def plausible_peptide_aa_range(neutral_mass: float) -> tuple[int, int] | None:
    """
    Estimate the range of possible residue counts (n) consistent with a
    linear peptide of this neutral mass, using min/max single-residue mass
    bounds. Returns None if the mass is below the smallest possible
    dipeptide (2 x Gly + H2O).
    """
    if neutral_mass < 2 * MIN_RESIDUE + H2O:
        return None
    n_min_possible = max(2, int((neutral_mass - H2O) / MAX_RESIDUE))
    n_max_possible = int((neutral_mass - H2O) / MIN_RESIDUE)
    return (max(2, n_min_possible), n_max_possible)


PEPTIDE_11_50_MASS_WINDOW = (11 * MIN_RESIDUE + H2O, 50 * MAX_RESIDUE + H2O)
OLIGOPEPTIDE_MASS_WINDOW = (51 * MIN_RESIDUE + H2O, 300 * MAX_RESIDUE + H2O)

# Minimum composite score (see linear_peptide_scoring.py) required for the
# 2-10 AA tier to actually call something a peptide. Precursor mass alone
# is NOT sufficient: at masses in the 700-1500 Da range there can be
# hundreds of amino-acid compositions matching within 0.02 Da purely by
# coincidence (see compositions_for_linear_mass docstring) -- a sugar,
# lipid, or polysaccharide fragment can accidentally satisfy one of those
# compositions on MASS alone. Requiring the *fragment ions themselves* to
# support b/y/a/immonium coverage above this threshold is what actually
# discriminates a real peptide call from a coincidental mass match.
PEPTIDE_CLASSIFICATION_MIN_SCORE = 0.40


def classify_family(mz: float, adduct: str, charge: int, ccs: float,
                     mz_obs: np.ndarray, int_obs: np.ndarray,
                     max_perms_per_composition: int = 15,
                     max_compositions: int = 15,
                     n_decoys: int = 15,
                     tol: float = 0.02) -> tuple[str, float, tuple | None, list]:
    """
    First-pass family triage for one feature. Tiered strategy:
      1. 2-10 AA: requires BOTH a matching precursor-mass composition AND
         fragment-ion evidence (composite score above threshold) -- mass
         match alone is not proof, since coincidental composition matches
         are common at mid-range masses (see module-level constant above).
      2. Explicit small-molecule family windows (pesticide/sugar/lipid/
         polysaccharide), checked next.
      3. Peptide_11-50AA / Oligopeptide_>50AA / Protein -- mass-window
         heuristic only (no exhaustive search; see module docstring).

    Returns (family_label, neutral_mass, aa_range_or_None, candidates_list).
    candidates_list is the scored candidate list when the 2-10 AA tier ran
    scoring (empty otherwise), so the caller can reuse it for Sheet 3
    instead of re-scoring the same feature twice.
    """
    neutral_mass = neutral_mass_from_row(mz, adduct, charge)
    aa_range = plausible_peptide_aa_range(neutral_mass)

    if (2 * MIN_RESIDUE + H2O <= neutral_mass <= PEPTIDE_MAX_AA * MAX_RESIDUE + H2O
            and len(mz_obs) > 0):
        candidates = identify_linear_peptide(
            mz_obs, int_obs, mz, charge=charge,
            n_min=PEPTIDE_MIN_AA, n_max=PEPTIDE_MAX_AA, tol=tol,
            min_score=0.10, n_decoys=n_decoys, fdr_alpha=0.05,
            max_perms_per_composition=max_perms_per_composition,
            max_compositions=max_compositions,
        )
        if candidates and candidates[0].composite >= PEPTIDE_CLASSIFICATION_MIN_SCORE:
            return "Peptide_2-10AA", neutral_mass, aa_range, candidates

    for family, (lo, hi) in FAMILY_MASS_WINDOWS.items():
        if family == "Protein":
            continue
        if lo <= neutral_mass <= hi:
            return family, neutral_mass, None, []

    lo11, hi11 = PEPTIDE_11_50_MASS_WINDOW
    if lo11 <= neutral_mass <= hi11:
        return "Peptide_11-50AA", neutral_mass, aa_range, []

    lo_oligo, hi_oligo = OLIGOPEPTIDE_MASS_WINDOW
    if lo_oligo <= neutral_mass <= hi_oligo:
        return "Oligopeptide_>50AA", neutral_mass, aa_range, []

    if neutral_mass > FAMILY_MASS_WINDOWS["Protein"][0]:
        return "Protein", neutral_mass, None, []

    return "Unclassified", neutral_mass, None, []


# ══════════════════════════════════════════════════════════════════════════
# MAIN PROCESSING
# ══════════════════════════════════════════════════════════════════════════

def process_feature_table(input_path: str, output_path: str,
                           n_decoys: int = 15, tol: float = 0.02,
                           max_perms_per_composition: int = 15,
                           max_compositions: int = 15) -> None:
    df = pd.read_excel(input_path, sheet_name="Features")

    families, neutral_masses, aa_ranges = [], [], []
    candidates_by_feature: dict[str, list] = {}

    for i, row in df.iterrows():
        mz = float(row["m/z"])
        adduct = row["Adduct"]
        charge = int(row["Charge"])
        ccs = float(row["CCS_A2"])
        mz_obs, int_obs = parse_fragments(row.get("MS_Fragments (mz:intensity)", ""))

        print(f"  classifying {row['Feature_ID']} ({i+1}/{len(df)})...", flush=True)
        fam, nm, aar, candidates = classify_family(
            mz, adduct, charge, ccs, mz_obs, int_obs,
            max_perms_per_composition=max_perms_per_composition,
            max_compositions=max_compositions,
            n_decoys=n_decoys, tol=tol,
        )
        families.append(fam)
        neutral_masses.append(round(nm, 4))
        aa_ranges.append(aar)
        if candidates:
            candidates_by_feature[row["Feature_ID"]] = candidates

    df["Predicted_family"] = families
    df["Neutral_mass_calc"] = neutral_masses
    df["Plausible_AA_range"] = [f"{a[0]}-{a[1]}" if a else "" for a in aa_ranges]

    # ── Sheet 1: family summary ─────────────────────────────────────────
    summary = (
        df["Predicted_family"].value_counts()
        .rename_axis("Family").reset_index(name="N_features")
    )
    summary["Pct_of_total"] = round(100 * summary["N_features"] / len(df), 1)
    total_row = pd.DataFrame([{"Family": "TOTAL", "N_features": len(df), "Pct_of_total": 100.0}])
    summary = pd.concat([summary, total_row], ignore_index=True)

    # ── Sheet 2: oligopeptides > 50 AA (flagged, not scored) ────────────
    oligo_df = df[df["Predicted_family"] == "Oligopeptide_>50AA"][
        ["Feature_ID", "RT_min", "m/z", "Adduct", "Charge", "CCS_A2",
         "Neutral_mass_calc", "Plausible_AA_range"]
    ].copy()
    oligo_df["Note"] = "Consistent with linear peptide >50 AA — outside 2-10 AA scoring scope of this tool"

    # ── Sheet 3: 2-10 AA peptides — reuse scoring already done during
    #    classification (classify_family had to score these to decide the
    #    family, so there's no need to run identify_linear_peptide again) ──
    pep_df = df[df["Predicted_family"] == "Peptide_2-10AA"].copy()
    from linear_peptide_scoring import confidence_tier

    scored_rows = []
    for _, row in pep_df.iterrows():
        fid = row["Feature_ID"]
        candidates = candidates_by_feature.get(fid, [])

        if not candidates:
            scored_rows.append({
                "Feature_ID": fid, "RT_min": row["RT_min"],
                "m/z": row["m/z"], "Adduct": row["Adduct"], "CCS_A2": row["CCS_A2"],
                "Neutral_mass_calc": row["Neutral_mass_calc"],
                "Top_candidate": "", "Composite": np.nan, "Tier": "NO_MATCH",
                "q_storey": np.nan, "rejected_BH": False,
                "N_candidates_scored": 0, "Note": "No scored candidates retained",
            })
            continue

        top = candidates[0]
        ties = [c.sequence for c in candidates
                if abs(c.composite - top.composite) < 1e-6]
        tie_note = f"Isobaric tie: {', '.join(ties)}" if len(ties) > 1 else ""

        scored_rows.append({
            "Feature_ID": fid, "RT_min": row["RT_min"],
            "m/z": row["m/z"], "Adduct": row["Adduct"], "CCS_A2": row["CCS_A2"],
            "Neutral_mass_calc": row["Neutral_mass_calc"],
            "Top_candidate": top.sequence, "Composite": top.composite,
            "Tier": confidence_tier(top.composite),
            "q_storey": round(top.q_storey, 4) if not np.isnan(top.q_storey) else np.nan,
            "rejected_BH": top.rejected_BH,
            "N_candidates_scored": len(candidates),
            "Note": tie_note,
        })

    pep_results_df = pd.DataFrame(scored_rows)

    # ── Sheet 4: ion-level evidence — the audit trail behind each score.
    #    For every feature's top candidate, break down exactly which
    #    a/b/y/immonium ions matched, with theoretical vs. observed m/z
    #    and the mass error (Delta Da / Delta ppm) -- this is what lets
    #    a researcher see WHY a sequence scored the way it did, not just
    #    the collapsed composite number. ─────────────────────────────────
    from linear_peptide_scoring import ion_evidence_table

    evidence_rows = []
    for _, row in pep_df.iterrows():
        fid = row["Feature_ID"]
        candidates = candidates_by_feature.get(fid, [])
        if not candidates:
            continue
        top = candidates[0]
        mz_obs, int_obs = parse_fragments(row.get("MS_Fragments (mz:intensity)", ""))
        if len(mz_obs) == 0:
            continue
        for ion_row in ion_evidence_table(mz_obs, int_obs, top.sequence, tol=tol, charge=1):
            evidence_rows.append({"Feature_ID": fid, "Sequence": top.sequence, **ion_row})

    ion_evidence_df = pd.DataFrame(evidence_rows)

    # ── Write multi-sheet workbook ───────────────────────────────────────
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary_Families", index=False)
        oligo_df.to_excel(writer, sheet_name="Oligopeptides_gt50AA", index=False)
        pep_results_df.to_excel(writer, sheet_name="Peptides_2to10AA", index=False)
        ion_evidence_df.to_excel(writer, sheet_name="Ion_Evidence", index=False)

    # ── Console report ───────────────────────────────────────────────────
    print(f"Processed {len(df)} features from {input_path}\n")
    print("Family breakdown:")
    print(summary.to_string(index=False))
    print(f"\nOligopeptides >50 AA flagged: {len(oligo_df)}")
    print(f"Peptides 2-10 AA scored: {len(pep_results_df)}")
    if len(pep_results_df) > 0:
        tier_counts = pep_results_df["Tier"].value_counts()
        print(f"  Confidence tier breakdown:\n{tier_counts.to_string()}")
    print(f"\nResults written to: {output_path}")


if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "UHPLC_RP_Wine_Template.xlsx"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "WineFeature_Results.xlsx"
    process_feature_table(input_file, output_file)
