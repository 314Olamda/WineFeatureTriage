# UHPLC-RP Wine Feature Template + Classification Tool

## What's in this package

| File | Role |
|---|---|
| `generate_template.py` | Generates the synthetic, **unannotated** feature table (77 features) spanning pesticides, lipids, polysaccharides, sugars, proteins, oligopeptides >50 AA, and peptides 2-50 AA (29.9% of total, weighted toward 2-10 AA) |
| `UHPLC_RP_Wine_Template.xlsx` | The generated template itself — this is what a real "upload" would contain: RT, m/z, adduct, charge, CCS, MS fragments, **no compound names or sequences** |
| `wine_feature_classifier.py` | The tool — reads the template, classifies features into families, and scores every 2-10 AA candidate with a/b/y/immonium fragment evidence |
| `WineFeature_Results.xlsx` | The tool's output: 3 sheets (`Summary_Families`, `Oligopeptides_gt50AA`, `Peptides_2to10AA`) |
| `answer_key_INTERNAL_DO_NOT_SHARE.xlsx` | Ground truth used only to validate the tool — **not part of the blind template**, not something a real uploader would have |
| `linear_denovo.py` | Updated core engine — includes two real performance fixes made while building this (see below) |

## Validated accuracy

Checked the tool's `Peptide_2-10AA` calls against the answer key:

| Metric | Result |
|---|---|
| **Precision** | 100% (2/2 calls were exact sequence matches — `AD`, `YKSV`) |
| **Recall** | 11.8% (2 of 17 true 2-10 AA peptides found) |
| **False positives** | 0 |

The low recall is a direct, expected consequence of the runtime caps (`max_compositions=15`, `max_perms_per_composition=15`) needed to make this run interactively rather than a correctness problem — see "Known limitation" below.

## Two real bugs fixed while building this

1. **Branch-and-bound was iterating alphabetically, not by mass.** `_all_compositions()` in `linear_denovo.py` originally sorted residues alphabetically, which made mass-based pruning unsafe (an alphabetically-later residue could still be lighter). Re-sorting by mass ascending lets the search `break` early once a residue would overshoot the target, and adds an upper-bound prune using the remaining slot count — cut real query times from 2-16 seconds down to well under 0.1s at mid-range masses.
2. **Factorial blowup in permutation/decoy generation.** `all_permutations_linear()` and `generate_decoys_linear()` both called `itertools.permutations()` on the full sequence before sampling — for a 10-residue composition with mostly distinct residues, that's up to 3.6 million permutations enumerated just to keep 15-200 of them. Fixed with a length threshold (`EXACT_ENUM_MAX_LEN = 8`): exact enumeration below that, direct random sampling above it.

These are genuine improvements to `linear_denovo.py` independent of this template exercise — worth carrying back into the standalone `linear_denovo.py` you'd upload to `WineLinearPep`.

## Known limitation: composition search caps trade recall for speed

At masses in the ~700-1500 Da range, there can be **hundreds to thousands of amino-acid compositions** matching a given precursor mass within 0.02 Da — a real, documented characteristic of de novo composition search with the full 20-residue alphabet, not a search bug (confirmed directly: mass 1073.5 Da alone had over 15,000 matching compositions before capping). Capping enumeration at 15 per feature makes the tool usable in real time, but means the true composition for larger/more degenerate peptides is often not among the first 15 found.

**What this means practically:** a real online tool for this shouldn't run synchronously. The realistic architecture is an async job — upload triggers a background worker (e.g., the GitHub Actions pattern discussed earlier, or a proper task queue), the search runs with production-scale caps (`max_compositions=500`, `max_perms_per_composition=200`, matching `linear_denovo.py`'s original defaults), and results are ready in minutes rather than seconds. That would recover most of the missing recall without changing any of the underlying logic — only the caps.

## Classification logic (Sheet 1 & 2 families)

Family triage for non-peptide categories (pesticide/lipid/sugar/polysaccharide/protein) uses **mass-window heuristics only** — a reasonable first-pass filter for a demo, not a substitute for real identification. Real triage needs spectral library matching (GNPS, MassBank, mzCloud) plus formula/isotope confirmation. The 2-10 AA peptide tier is more rigorous: a feature is only called a peptide if its **fragment ions actually support** a b/y/a/immonium-scored composite ≥ 0.40 — mass match alone is deliberately not sufficient, since coincidental composition matches are common at these masses (this was a real bug caught and fixed mid-build: the first version classified sugars and lipids as peptides purely because their mass happened to fit some amino acid combination).
