#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
linear_peptide_scoring.py — a/b/y + immonium scoring for linear peptides
══════════════════════════════════════════════════════════════════════════
Sibling module to WineCycloPep's score_spectrum_vs_cyclic, adapted for
LINEAR peptides (2-50 AA, optimized for the 2-10 AA range discussed
earlier). Built on mass_utils.py (Pyteomics-backed) rather than a
hand-rolled mass table, so it inherits the same fix already applied to
the cyclic pipeline.

PATCH NOTE (ppm-scaled tolerance)
──────────────────────────────────
The original version matched every ion (b/a/y/immonium, all masses from
~70 Da immonium ions to >1000 Da precursor-range fragments) against a
single flat `tol` in Da (default 0.02). Real instrument mass error scales
with m/z, not as a constant Da offset -- this is exactly what the
Ion_Evidence sheet caught in practice: a MEDIUM-tier call had several ions
matching only at 21-27 ppm error while still passing the fixed 0.02 Da
window, well outside the sub-1 ppm seen on HIGH-tier calls. A flat Da
tolerance is simultaneously too loose at high mass (lets in spurious
matches) and too tight at low mass (can reject real matches at the
low-mass end, e.g. immonium ions).

Fix: `tol_ppm` (default 15 ppm, typical for QTOF/Orbitrap-class
instruments) is now the primary matching window, scaled per-ion as
theo_mz * tol_ppm * 1e-6. An optional `tol_da_floor` sets a minimum
absolute window so very low-mass ions aren't held to an unrealistically
tiny window even in ppm terms (e.g. at 15 ppm, a 70 Da immonium ion gets
only +-0.00105 Da -- probably tighter than your instrument's real
resolution at that mass; a floor of ~0.005-0.01 Da is more realistic).

Both `score_spectrum_vs_linear` and `ion_evidence_table` now expose
tol_ppm/tol_da_floor instead of a single flat `tol`. See the note at the
bottom of this file for what still needs updating in wine_feature_classifier.py
and (likely) linear_denovo.py, which I have not patched here since I
haven't seen that file's source -- happy to do it if you want the fix
propagated end-to-end.

Why the scoring logic differs from the cyclic version
────────────────────────────────────────────────────
Cyclic peptides (WineCycloPep):
  • n possible ring-opening positions -> n overlapping bn-ion ladders
  • NO free C-terminus -> y1 ABSENCE is diagnostic of cyclicity
  • Immonium ions: composition check only (low weight, 0.10)

Linear peptides (this module):
  • ONE fragmentation series -> single b-ion ladder, single y-ion ladder
  • Free C-terminus -> y-ion PRESENCE is diagnostic (opposite logic)
  • a-ions (b - CO) are a genuine third series worth scoring separately --
    for cyclic peptides a-ions are largely redundant with the multi-
    rotation b-ion ladder, but for linear peptides they're an independent
    confirmation, so they get their own scoring dimension here rather
    than being folded into "loss ions" as in the cyclic script.
  • Immonium ions play the same confirmatory (low-weight) role as in
    WineCycloPep -- composition validation, not sequence-order evidence.

This module intentionally covers the SCORING function only -- not a full
Bruker .d reading / de novo composition search pipeline. If you want the
full WineCycloPep-style pipeline (raw .d ingestion, composition search,
FDR/decoy validation, 3D structures) built out for linear peptides, that's
a separate, larger project -- this gives you the core scoring engine to
drop into whatever spectrum source you're already using (e.g. your
MS-DIAL / GNPS exports).

Requirements
────────────
pip install pyteomics numpy
mass_utils.py must be in the same directory (or on PYTHONPATH)

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations
import numpy as np

from mass_utils import (
    RESIDUE_MASS, PROTON, H2O,
    linear_neutral_mass, precursor_neutral,
    immonium_ions,
)

# ══════════════════════════════════════════════════════════════════════════
# ION SERIES — a, b, y for a LINEAR peptide (single series, not rotated)
# ══════════════════════════════════════════════════════════════════════════

def b_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    b-ion series for a linear peptide: N-terminal fragments.
    b_i = sum(residues[0:i]) + proton, for i = 1 .. n-1
    """
    n = len(sequence)
    mz = []
    running = 0.0
    for i in range(n - 1):  # b_n (full length) is not a real fragment
        running += RESIDUE_MASS.get(sequence[i], 0.0)
        mz.append(round((running + PROTON * charge) / charge, 5))
    return mz


def a_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    a-ion series: b_i - CO. Independent confirmatory series for linear
    peptides (unlike the cyclic case, where a-ions are largely redundant
    with the multi-rotation b-ion ladder).
    """
    from pyteomics import mass as _mass
    co = _mass.calculate_mass(formula="CO")
    return [round(b - co / charge, 5) for b in b_ions(sequence, charge)]


def y_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    y-ion series: C-terminal fragments, INCLUDING y1 (unlike the cyclic
    scorer, where y1 presence is treated as evidence AGAINST cyclicity).
    y_i = sum(residues[n-i:n]) + H2O + proton, for i = 1 .. n-1
    """
    n = len(sequence)
    mz = []
    running = 0.0
    for i in range(1, n):
        running += RESIDUE_MASS.get(sequence[n - i], 0.0)
        mz.append(round((running + H2O + PROTON * charge) / charge, 5))
    return mz


# ══════════════════════════════════════════════════════════════════════════
# PPM-SCALED MATCHING TOLERANCE
# ══════════════════════════════════════════════════════════════════════════

DEFAULT_TOL_PPM = 15.0       # typical high-res QTOF/Orbitrap mass accuracy
DEFAULT_TOL_DA_FLOOR = 0.008  # absolute floor for very low-mass ions


def match_window(theo_mz: float, tol_ppm: float = DEFAULT_TOL_PPM,
                  tol_da_floor: float = DEFAULT_TOL_DA_FLOOR) -> float:
    """
    Resolve the matching half-window (Da) for one theoretical ion m/z.

    Primary scaling is ppm (window = theo_mz * tol_ppm * 1e-6), since real
    instrument mass error scales with m/z. tol_da_floor sets a minimum
    absolute window so low-mass ions (immonium ions ~70-140 Da) aren't held
    to an unrealistically tiny window purely because ppm math shrinks fast
    at low mass -- e.g. at 15 ppm, a 70 Da ion gets +-0.00105 Da, likely
    tighter than real achievable resolution at that mass.
    """
    ppm_window = theo_mz * tol_ppm * 1e-6
    return max(ppm_window, tol_da_floor)


# ══════════════════════════════════════════════════════════════════════════
# SCORING
# ══════════════════════════════════════════════════════════════════════════

# Weights sum to 1.0. b/y carry the most sequence-order evidence; a-ions
# and immonium ions are confirmatory. Intensity rewards spectra where the
# matched ions are actually the dominant peaks, not just present above noise.
SCORE_WEIGHTS = {
    "b_coverage": 0.30,
    "y_coverage": 0.30,
    "a_coverage": 0.15,
    "intensity_score": 0.15,
    "immonium_score": 0.10,
}


def score_spectrum_vs_linear(
    mz_obs: np.ndarray,
    int_obs: np.ndarray,
    sequence: str,
    tol_ppm: float = DEFAULT_TOL_PPM,
    tol_da_floor: float = DEFAULT_TOL_DA_FLOOR,
    charge: int = 1,
) -> dict:
    """
    Score an MS2 spectrum against a linear peptide candidate sequence
    using a/b/y ion coverage + immonium confirmation.

    Matching tolerance is now ppm-scaled per ion (see match_window above)
    instead of a single flat Da value -- see module PATCH NOTE.

    Mirrors score_spectrum_vs_cyclic's structure and return-dict shape
    (composite + per-dimension scores + matched/total counts) so it can
    be dropped into the same downstream FDR/decoy pipeline pattern used
    in WineCycloPep, if you extend this to a full de novo search later.
    """
    if len(mz_obs) == 0:
        return {"composite": 0.0}

    total_int = float(np.sum(int_obs))

    b_theo = b_ions(sequence, charge)
    a_theo = a_ions(sequence, charge)
    y_theo = y_ions(sequence, charge)
    imm_theo = list(immonium_ions(sequence).values())

    def _coverage(theo: list[float]) -> tuple[int, int, float]:
        if not theo:
            return 0, 0, 0.0
        matched = sum(
            1 for mz in theo
            if np.any(np.abs(mz_obs - mz) <= match_window(mz, tol_ppm, tol_da_floor))
        )
        return matched, len(theo), matched / len(theo)

    b_matched, b_total, b_coverage = _coverage(b_theo)
    a_matched, a_total, a_coverage = _coverage(a_theo)
    y_matched, y_total, y_coverage = _coverage(y_theo)
    imm_matched, imm_total, immonium_score = _coverage(imm_theo)

    # Intensity score: matched peaks' share of total ion current
    all_theo = b_theo + a_theo + y_theo
    matched_int = 0.0
    for mz in all_theo:
        w = match_window(mz, tol_ppm, tol_da_floor)
        hits = np.where(np.abs(mz_obs - mz) <= w)[0]
        if len(hits) > 0:
            matched_int += float(np.max(int_obs[hits]))
    intensity_score = matched_int / total_int if total_int > 0 else 0.0

    composite = (
        SCORE_WEIGHTS["b_coverage"] * b_coverage +
        SCORE_WEIGHTS["y_coverage"] * y_coverage +
        SCORE_WEIGHTS["a_coverage"] * a_coverage +
        SCORE_WEIGHTS["intensity_score"] * intensity_score +
        SCORE_WEIGHTS["immonium_score"] * immonium_score
    )

    return {
        "composite": round(composite, 4),
        "b_coverage": round(b_coverage, 3),
        "y_coverage": round(y_coverage, 3),
        "a_coverage": round(a_coverage, 3),
        "intensity_score": round(intensity_score, 3),
        "immonium_score": round(immonium_score, 3),
        "b_matched": b_matched, "b_total": b_total,
        "y_matched": y_matched, "y_total": y_total,
        "a_matched": a_matched, "a_total": a_total,
        "immonium_matched": imm_matched, "immonium_total": imm_total,
    }


def confidence_tier(score: float) -> str:
    """Same tier boundaries as WineCycloPep, for consistency across the series."""
    if score >= 0.60:
        return "HIGH"
    if score >= 0.35:
        return "MEDIUM"
    if score >= 0.15:
        return "LOW"
    return "VERY_LOW"


def ion_evidence_table(
    mz_obs: np.ndarray,
    int_obs: np.ndarray,
    sequence: str,
    tol_ppm: float = DEFAULT_TOL_PPM,
    tol_da_floor: float = DEFAULT_TOL_DA_FLOOR,
    charge: int = 1,
) -> list[dict]:
    """
    Build a labeled, per-ion diagnostic table explaining WHY a sequence
    scored the way it did -- the a1/a2/b1/b2/y1/y2/... breakdown with
    theoretical m/z, the closest observed m/z (if matched), the mass
    error (Delta Da and Delta ppm), and whether it counted as a match.

    Now reports the actual per-ion matching window used (`match_window_Da`)
    alongside delta_ppm, so a borderline call is visible at a glance
    instead of requiring a manual check against a hidden flat tolerance.

    This is the audit trail behind score_spectrum_vs_linear's composite
    score -- same ion series (b/a/y/immonium), same ppm-scaled tolerance,
    but reported per-ion instead of collapsed into coverage fractions, so
    a researcher can see exactly which ions did (or didn't) support a call.
    """
    from mass_utils import immonium_ions

    rows: list[dict] = []
    n = len(sequence)

    def _closest_match(theo_mz: float) -> tuple[float | None, float | None, float]:
        """Return (observed_mz, observed_intensity, window_Da) for the
        closest peak within the ppm-scaled window, or (None, None, window)
        if nothing matches."""
        w = match_window(theo_mz, tol_ppm, tol_da_floor)
        if len(mz_obs) == 0:
            return None, None, w
        diffs = np.abs(mz_obs - theo_mz)
        idx = np.argmin(diffs)
        if diffs[idx] <= w:
            return float(mz_obs[idx]), float(int_obs[idx]), w
        return None, None, w

    def _add_ion(label: str, ion_type: str, theo_mz: float) -> None:
        obs_mz, obs_int, window_da = _closest_match(theo_mz)
        matched = obs_mz is not None
        delta_da = round(obs_mz - theo_mz, 5) if matched else None
        delta_ppm = round(1e6 * (obs_mz - theo_mz) / theo_mz, 2) if matched else None
        rows.append({
            "ion_label": label,
            "ion_type": ion_type,
            "theoretical_mz": round(theo_mz, 5),
            "observed_mz": round(obs_mz, 5) if matched else None,
            "delta_Da": delta_da,
            "delta_ppm": delta_ppm,
            "match_window_Da": round(window_da, 5),
            "observed_intensity": round(obs_int, 1) if matched else None,
            "matched": matched,
        })

    # b/a ions: b1, b2, ... b(n-1) and a1, a2, ... a(n-1)
    for i, mz in enumerate(b_ions(sequence, charge), start=1):
        _add_ion(f"b{i}", "b", mz)
    for i, mz in enumerate(a_ions(sequence, charge), start=1):
        _add_ion(f"a{i}", "a", mz)

    # y ions: y1 (C-terminal residue) ... y(n-1), matching the indexing
    # convention used in y_ions() (i-th element = y_i, counted from C-term)
    for i, mz in enumerate(y_ions(sequence, charge), start=1):
        _add_ion(f"y{i}", "y", mz)

    # Immonium ions, one per unique residue, labeled by residue letter
    for aa, mz in immonium_ions(sequence).items():
        _add_ion(f"imm_{aa}", "immonium", mz)

    return rows


# ══════════════════════════════════════════════════════════════════════════
# SELF-TEST
# ══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Synthetic test: perfect spectrum for linear tripeptide "IAA"
    seq = "IAA"
    b_theo = b_ions(seq)
    y_theo = y_ions(seq)
    imm_theo = list(immonium_ions(seq).values())

    print(f"Sequence: {seq} (linear)")
    print(f"Linear neutral mass: {linear_neutral_mass(seq):.5f} Da")
    print(f"b-ions: {b_theo}")
    print(f"a-ions: {a_ions(seq)}")
    print(f"y-ions: {y_theo}")
    print(f"Immonium ions: {imm_theo}")

    # Build a synthetic "perfect" spectrum from theoretical ions
    mz_obs = np.array(b_theo + y_theo + imm_theo)
    int_obs = np.ones_like(mz_obs) * 1000.0

    result = score_spectrum_vs_linear(mz_obs, int_obs, seq)
    print(f"\nScore: {result['composite']} -> {confidence_tier(result['composite'])}")
    print(result)

    # Borderline test: shift one b-ion by 20 ppm to show the ppm scaling
    # catching what a flat 0.02 Da tolerance would have silently passed
    # at high mass and rejected at low mass.
    shifted = mz_obs.copy()
    shifted[0] = shifted[0] * (1 + 20e-6)  # +20 ppm on the first b-ion
    result_shifted = score_spectrum_vs_linear(shifted, int_obs, seq)
    print(f"\n+20ppm shift on first ion -> composite {result_shifted['composite']} "
          f"(window at that mass: {match_window(mz_obs[0]):.5f} Da)")
