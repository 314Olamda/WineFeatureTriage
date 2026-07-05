#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_template.py — Synthetic UHPLC-RP feature table generator
════════════════════════════════════════════════════════════════════

Generates a realistic, UNANNOTATED feature table (RT, m/z, adduct, CCS,
MS fragments) spanning the major chemical families expected in a wine
lees UHPLC-RP untargeted run:

  - Pesticides (grape-relevant, ~15 common actives)
  - Lipids (yeast: ergosterol/PLs/TAGs; grape: free fatty acids)
  - Polysaccharides (yeast mannoprotein fragments; grape pectin fragments)
  - Sugars (grape: glucose/fructose/sucrose/trehalose/raffinose)
  - Proteins (mainly yeast, some grape — deconvolved high-mass species)
  - Oligopeptides > 50 AA (mainly yeast autolysis products)
  - Peptides 2-50 AA (30% of total features), weighted toward 2-10 AA

Peptide-range rows use REAL fragment ion physics — b/y/a/immonium ions
are computed with mass_utils.py / linear_peptide_scoring.py for randomly
generated sequences, so the downstream classifier/scorer tool has
genuine signal to recover, not placeholder numbers.

The exported template has NO compound name, sequence, or family labels
— it simulates a raw feature list as it would arrive before annotation.
A separate "answer key" file is written alongside it for validating the
classifier tool's output, and would NOT normally be shared with whoever
is uploading real unknown data.

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations

import random
import numpy as np
import pandas as pd

from mass_utils import RESIDUE_MASS, PROTON, H2O, linear_neutral_mass
from linear_peptide_scoring import b_ions, a_ions, y_ions
from mass_utils import immonium_ions

random.seed(42)
np.random.seed(42)

# ══════════════════════════════════════════════════════════════════════════
# 1. PESTICIDES — common grape-relevant actives (monoisotopic mass, Da)
# ══════════════════════════════════════════════════════════════════════════
PESTICIDES = {
    "imidacloprid":   255.0524, "azoxystrobin":   403.1305,
    "boscalid":       342.0483, "cyprodinil":     225.1266,
    "fenhexamid":     301.0637, "tebuconazole":   307.1451,
    "chlorpyrifos":   348.9296, "dimethomorph":   387.1450,
    "folpet":         295.9375, "iprodione":      329.0177,
    "metalaxyl":      279.1471, "myclobutanil":   288.1091,
    "penconazole":    283.0824, "pyrimethanil":   199.1109,
    "spinosad_a":     731.4614, "thiamethoxam":   291.0522,
    "trifloxystrobin":408.1347,
}

# ══════════════════════════════════════════════════════════════════════════
# 2. LIPIDS — yeast + grape (monoisotopic mass, Da)
# ══════════════════════════════════════════════════════════════════════════
LIPIDS = {
    "ergosterol_yeast":      396.3392, "oleic_acid":            282.2559,
    "palmitic_acid":         256.2402, "linoleic_acid":         280.2402,
    "linolenic_acid_grape":  278.2246, "PC_34_1_yeast":         759.5778,
    "PE_34_1_yeast":         717.5309, "TAG_52_2_yeast":        858.7838,
    "sitosterol_grape":      414.3862, "campesterol_grape":     400.3705,
}

# ══════════════════════════════════════════════════════════════════════════
# 3. POLYSACCHARIDE FRAGMENTS — mannose/glucose/galacturonic oligomers
# ══════════════════════════════════════════════════════════════════════════
def hexose_oligomer_mass(n: int) -> float:
    """n hexose units linked by glycosidic bonds (n*162.0528 + 18.0106)."""
    return n * 162.0528 + H2O

def galacturonic_oligomer_mass(n: int) -> float:
    """n galacturonic acid units (pectin backbone), n*176.0321 + H2O."""
    return n * 176.0321 + H2O

POLYSACCHARIDES = {
    f"mannooligomer_DP{n}_yeast": hexose_oligomer_mass(n) for n in range(2, 7)
}
POLYSACCHARIDES.update({
    f"galacturonan_DP{n}_grape": galacturonic_oligomer_mass(n) for n in range(2, 6)
})

# ══════════════════════════════════════════════════════════════════════════
# 4. SUGARS — grape
# ══════════════════════════════════════════════════════════════════════════
SUGARS = {
    "glucose": 180.0634, "fructose": 180.0634, "sucrose": 342.1162,
    "trehalose_yeast": 342.1162, "raffinose": 504.1690,
}

# ══════════════════════════════════════════════════════════════════════════
# 5. PROTEINS — deconvolved intact masses (mainly yeast, some grape)
# ══════════════════════════════════════════════════════════════════════════
PROTEINS = {
    "invertase_fragment_yeast": 27400.5, "HSP12_yeast": 11700.2,
    "thaumatin_like_grape": 24100.8, "chitinase_grape": 26300.4,
    "mannoprotein_yeast": 32500.6,
}

# ══════════════════════════════════════════════════════════════════════════
# 6. RANDOM PEPTIDE / OLIGOPEPTIDE GENERATOR
# ══════════════════════════════════════════════════════════════════════════
CANONICAL_AA = list(RESIDUE_MASS.keys())

def random_sequence(n_aa: int) -> str:
    return "".join(random.choice(CANONICAL_AA) for _ in range(n_aa))

def peptide_fragments_string(seq: str, max_ions: int = 12) -> str:
    """Build a realistic 'mz:intensity' fragment string using real b/y/a/immonium ions."""
    ions = {}
    for mz in b_ions(seq):
        ions[round(mz, 4)] = random.uniform(500, 20000)
    for mz in y_ions(seq):
        ions[round(mz, 4)] = random.uniform(500, 25000)
    for mz in a_ions(seq):
        ions[round(mz, 4)] = random.uniform(200, 8000)
    for mz in immonium_ions(seq).values():
        ions[round(mz, 4)] = random.uniform(300, 15000)
    # Cap the ion list length for very long peptides (realistic MS2 depth limit)
    items = sorted(ions.items(), key=lambda x: -x[1])[:max_ions]
    return ";".join(f"{mz:.4f}:{intensity:.0f}" for mz, intensity in items)

def oligopeptide_fragments_string(seq: str, max_ions: int = 8) -> str:
    """Sparse fragment coverage for long oligopeptides — realistic MS2 depth
    drops off steeply above ~15-20 residues, so only terminal ions are given."""
    ions = {}
    for mz in b_ions(seq)[:3]:
        ions[round(mz, 4)] = random.uniform(500, 10000)
    for mz in y_ions(seq)[:3]:
        ions[round(mz, 4)] = random.uniform(500, 12000)
    return ";".join(f"{mz:.4f}:{intensity:.0f}" for mz, intensity in ions.items())


# ══════════════════════════════════════════════════════════════════════════
# 7. ADDUCT / CCS SIMULATION
# ══════════════════════════════════════════════════════════════════════════
ADDUCTS_POS = ["[M+H]+", "[M+Na]+", "[M+NH4]+"]
ADDUCT_MASS_SHIFT = {"[M+H]+": PROTON, "[M+Na]+": 22.9892, "[M+NH4]+": 18.0338}

def simulate_ccs(neutral_mass: float, family: str) -> float:
    """
    Rough empirical CCS simulation (Å²) — NOT a validated predictor, just
    gives template rows plausible, family-differentiated CCS values based
    on typical trend lines reported for each compound class (peptides and
    lipids scale differently with mass due to shape/charge distribution).
    """
    base = {
        "pesticide": 130 + 0.28 * neutral_mass,
        "lipid": 180 + 0.22 * neutral_mass,
        "polysaccharide": 140 + 0.30 * neutral_mass,
        "sugar": 120 + 0.35 * neutral_mass,
        "protein": 300 + 0.09 * neutral_mass,
        "oligopeptide": 150 + 0.16 * neutral_mass,
        "peptide": 110 + 0.24 * neutral_mass,
    }[family]
    return round(base + random.uniform(-3, 3), 2)


# ══════════════════════════════════════════════════════════════════════════
# 8. BUILD THE FEATURE TABLE
# ══════════════════════════════════════════════════════════════════════════
rows = []
feature_id = 1
answer_key = []

def add_row(neutral_mass, family, fragments, rt, charge=1, adduct=None):
    global feature_id
    if adduct is None:
        adduct = random.choice(ADDUCTS_POS)
    mz = round((neutral_mass + ADDUCT_MASS_SHIFT[adduct]) / charge, 4)
    ccs = simulate_ccs(neutral_mass, family)
    fid = f"Feature_{feature_id:04d}"
    rows.append({
        "Feature_ID": fid,
        "RT_min": round(rt, 2),
        "m/z": mz,
        "Adduct": adduct,
        "Charge": charge,
        "CCS_A2": ccs,
        "MS_Fragments (mz:intensity)": fragments,
        "Compound_name": "",       # NO ANNOTATIONS — blind template
        "Sequence_peptide": "",    # NO ANNOTATIONS — blind template
        "Annotation_level": "",
    })
    feature_id += 1
    return fid

# --- Pesticides (~15 rows) ---
for name, mass_val in PESTICIDES.items():
    rt = round(random.uniform(4.0, 12.0), 2)
    frag_mz = round(mass_val * random.uniform(0.35, 0.75), 4)  # generic in-silico loss proxy
    frags = f"{frag_mz:.4f}:{random.uniform(1000,20000):.0f}"
    fid = add_row(mass_val, "pesticide", frags, rt)
    answer_key.append({"Feature_ID": fid, "True_family": "Pesticide", "True_identity": name})

# --- Lipids (~10 rows) ---
for name, mass_val in LIPIDS.items():
    rt = round(random.uniform(10.0, 18.0), 2)
    frag_mz = round(mass_val - 18.0106, 4)  # generic water-loss proxy fragment
    frags = f"{frag_mz:.4f}:{random.uniform(1000,15000):.0f}"
    fid = add_row(mass_val, "lipid", frags, rt)
    answer_key.append({"Feature_ID": fid, "True_family": "Lipid", "True_identity": name})

# --- Polysaccharide fragments (~9 rows) ---
for name, mass_val in POLYSACCHARIDES.items():
    rt = round(random.uniform(0.8, 2.5), 2)  # polar, elute early on RP
    frag_mz = round(mass_val - 162.0528, 4)  # generic hexose-loss proxy
    frags = f"{frag_mz:.4f}:{random.uniform(500,8000):.0f}"
    fid = add_row(mass_val, "polysaccharide", frags, rt)
    answer_key.append({"Feature_ID": fid, "True_family": "Polysaccharide", "True_identity": name})

# --- Sugars (~5 rows) ---
for name, mass_val in SUGARS.items():
    rt = round(random.uniform(0.5, 1.5), 2)
    frag_mz = round(mass_val - 18.0106, 4)
    frags = f"{frag_mz:.4f}:{random.uniform(2000,10000):.0f}"
    fid = add_row(mass_val, "sugar", frags, rt)
    answer_key.append({"Feature_ID": fid, "True_family": "Sugar", "True_identity": name})

# --- Proteins (~5 rows, high charge states typical of intact protein LC-MS) ---
for name, mass_val in PROTEINS.items():
    rt = round(random.uniform(14.0, 20.0), 2)
    charge = random.choice([15, 18, 22, 25])
    frags = ""  # intact-protein MS1 deconvolution typically has no MS2 fragment list here
    fid = add_row(mass_val, "protein", frags, rt, charge=charge, adduct="[M+H]+")
    answer_key.append({"Feature_ID": fid, "True_family": "Protein", "True_identity": name})

# --- Oligopeptides > 50 AA (~8 rows, mainly yeast) ---
for _ in range(8):
    n_aa = random.randint(51, 90)
    seq = random_sequence(n_aa)
    mass_val = linear_neutral_mass(seq)
    rt = round(random.uniform(12.0, 17.0), 2)
    frags = oligopeptide_fragments_string(seq)
    fid = add_row(mass_val, "oligopeptide", frags, rt, charge=random.choice([2, 3]))
    answer_key.append({"Feature_ID": fid, "True_family": "Oligopeptide_>50AA",
                        "True_identity": seq, "n_AA": n_aa})

# --- Peptides 2-50 AA (target: 30% of TOTAL features), weighted to 2-10 AA ---
n_fixed_so_far = feature_id - 1
target_total = int(round(n_fixed_so_far / 0.70))          # so peptides = 30% of total
n_peptide_rows = target_total - n_fixed_so_far
n_2to10 = int(round(n_peptide_rows * 0.75))                # majority of the peptide rows in 2-10 AA
n_11to50 = n_peptide_rows - n_2to10

for _ in range(n_2to10):
    n_aa = random.randint(2, 10)
    seq = random_sequence(n_aa)
    mass_val = linear_neutral_mass(seq)
    rt = round(random.uniform(1.0, 9.0), 2)
    frags = peptide_fragments_string(seq)
    fid = add_row(mass_val, "peptide", frags, rt, charge=1)
    answer_key.append({"Feature_ID": fid, "True_family": "Peptide_2-10AA",
                        "True_identity": seq, "n_AA": n_aa})

for _ in range(n_11to50):
    n_aa = random.randint(11, 50)
    seq = random_sequence(n_aa)
    mass_val = linear_neutral_mass(seq)
    rt = round(random.uniform(6.0, 15.0), 2)
    frags = oligopeptide_fragments_string(seq)
    fid = add_row(mass_val, "peptide", frags, rt, charge=random.choice([1, 2]))
    answer_key.append({"Feature_ID": fid, "True_family": "Peptide_11-50AA",
                        "True_identity": seq, "n_AA": n_aa})

# ══════════════════════════════════════════════════════════════════════════
# 9. SHUFFLE (real feature tables aren't sorted by family) AND EXPORT
# ══════════════════════════════════════════════════════════════════════════
random.shuffle(rows)
df = pd.DataFrame(rows)
df.sort_values("RT_min", inplace=True)  # feature tables ARE typically RT-sorted
df.reset_index(drop=True, inplace=True)

df.to_excel("UHPLC_RP_Wine_Template.xlsx", sheet_name="Features", index=False)

ak_df = pd.DataFrame(answer_key)
ak_df.to_excel("answer_key_INTERNAL_DO_NOT_SHARE.xlsx", index=False)

print(f"Total features generated: {len(df)}")
print(f"  Pesticides       : {len(PESTICIDES)}")
print(f"  Lipids           : {len(LIPIDS)}")
print(f"  Polysaccharides  : {len(POLYSACCHARIDES)}")
print(f"  Sugars           : {len(SUGARS)}")
print(f"  Proteins         : {len(PROTEINS)}")
print(f"  Oligopeptides>50AA: 8")
print(f"  Peptides 2-10 AA : {n_2to10}")
print(f"  Peptides 11-50 AA: {n_11to50}")
print(f"  Peptides 2-50 AA fraction of total (target 30%): {(n_2to10 + n_11to50) / len(df):.1%}")
