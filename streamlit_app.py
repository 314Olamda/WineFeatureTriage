#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
streamlit_app.py — Web UI for WineFeatureTriage
════════════════════════════════════════════════

Upload a filled UHPLC-RP feature template, the tool runs
wine_feature_classifier.py against it, and you download the
processed results workbook. No local Python required for the user.

DEBUG ARCHIVING
───────────────
Every uploaded file (and the results it produced) is automatically
committed to a private GitHub repo, so the app owner can review real
uploads later for debugging or improving the classifier. This requires
two Streamlit secrets (Settings -> Secrets in the Streamlit Cloud
dashboard):

    GH_TOKEN        = "github_pat_..."          (fine-grained, Contents: RW,
                                                   scoped to ONE private repo)
    GH_UPLOAD_REPO  = "314Olamda/WineFeatureTriage-uploads"

If these secrets aren't set, the app still works normally -- archiving
silently no-ops rather than breaking the tool for the person using it.

Deploy on Streamlit Community Cloud (share.streamlit.io):
  1. Push this file + wine_feature_classifier.py + linear_denovo.py +
     mass_utils.py + linear_peptide_scoring.py to the WineFeatureTriage repo
  2. Go to share.streamlit.io, sign in with GitHub, "New app"
  3. Point it at this repo, branch main, file streamlit_app.py
  4. Add requirements.txt (see repo root)
  5. Add the two secrets above under Settings -> Secrets
  6. Deploy

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
"""

import io
import base64
import tempfile
import os
from datetime import datetime, timezone

import streamlit as st
import pandas as pd
import requests

from wine_feature_classifier import process_feature_table

st.set_page_config(page_title="Wine FeatureTriage", page_icon="🍷", layout="centered")


# ══════════════════════════════════════════════════════════════════════════
# DEBUG ARCHIVING — commit uploaded files + results to a private GitHub repo
# ══════════════════════════════════════════════════════════════════════════

def _github_configured() -> bool:
    return "GH_TOKEN" in st.secrets and "GH_UPLOAD_REPO" in st.secrets


def archive_to_github(file_bytes: bytes, path_in_repo: str, commit_message: str) -> bool:
    """
    Commit a file to the private archive repo via the GitHub Contents API.
    Returns True on success, False on any failure (never raises -- archiving
    must never break the actual tool for the person using the app).
    """
    if not _github_configured():
        return False
    try:
        token = st.secrets["GH_TOKEN"]
        repo = st.secrets["GH_UPLOAD_REPO"]  # e.g. "314Olamda/WineFeatureTriage-uploads"
        url = f"https://api.github.com/repos/{repo}/contents/{path_in_repo}"
        content_b64 = base64.b64encode(file_bytes).decode("utf-8")
        resp = requests.put(
            url,
            headers={
                "Authorization": f"Bearer {token}",
                "Accept": "application/vnd.github+json",
            },
            json={"message": commit_message, "content": content_b64},
            timeout=15,
        )
        return resp.status_code in (200, 201)
    except Exception:
        return False  # archiving is best-effort, never surfaced to the user


def archive_session(input_bytes: bytes, input_name: str,
                     result_bytes, summary_text: str) -> None:
    """Archive one upload + its results (if produced) under a timestamped folder."""
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    folder = f"uploads/{ts}"
    archive_to_github(input_bytes, f"{folder}/{input_name}",
                       f"Upload archive {ts}: {input_name}")
    if result_bytes is not None:
        archive_to_github(result_bytes, f"{folder}/results.xlsx",
                           f"Results archive {ts}: {input_name}")
    if summary_text:
        archive_to_github(summary_text.encode("utf-8"), f"{folder}/run_info.txt",
                           f"Run info {ts}: {input_name}")


# ══════════════════════════════════════════════════════════════════════════
# APP UI
# ══════════════════════════════════════════════════════════════════════════

st.title("🍷 Wine FeatureTriage")
st.markdown(
    "Upload an unannotated UHPLC-RP feature table (RT, m/z, adduct, CCS, "
    "MS fragments) and get back family classification + 2-10 AA peptide "
    "scoring — no Python required on your end."
)

with st.expander("Expected input format"):
    st.markdown("""
    Excel file, sheet named `Features`, with columns:
    `Feature_ID`, `RT_min`, `m/z`, `Adduct`, `Charge`, `CCS_A2`,
    `MS_Fragments (mz:intensity)` (semicolon-separated `mz:intensity` pairs).
    See `example_data/UHPLC_RP_Wine_Template.xlsx` in the repo for a template.
    """)

st.caption(
    "Uploaded files are archived for debugging and improving this tool. "
    "Contact the maintainer if you'd prefer your upload not be retained."
)

st.divider()

uploaded_file = st.file_uploader("Upload your feature table (.xlsx)", type=["xlsx"])

col1, col2, col3 = st.columns(3)
with col1:
    n_decoys = st.number_input("Decoys per candidate", min_value=5, max_value=200, value=15)
with col2:
    max_perms = st.number_input("Max permutations/composition", min_value=5, max_value=500, value=15)
with col3:
    max_comps = st.number_input("Max compositions/feature", min_value=5, max_value=1000, value=15)

st.caption(
    "Higher values improve recall at the cost of runtime -- see the "
    "'Scaling for production' section in the repo README."
)

if uploaded_file is not None:
    if st.button("Run analysis", type="primary"):
        input_bytes = uploaded_file.getbuffer().tobytes()
        result_bytes = None
        run_error = None

        with st.spinner("Classifying features and scoring 2-10 AA candidates..."):
            with tempfile.TemporaryDirectory() as tmpdir:
                input_path = os.path.join(tmpdir, "input.xlsx")
                output_path = os.path.join(tmpdir, "results.xlsx")

                with open(input_path, "wb") as f:
                    f.write(input_bytes)

                try:
                    process_feature_table(
                        input_path, output_path,
                        n_decoys=int(n_decoys),
                        max_perms_per_composition=int(max_perms),
                        max_compositions=int(max_comps),
                    )
                    with open(output_path, "rb") as f:
                        result_bytes = f.read()
                except Exception as e:
                    run_error = str(e)

        # Archive regardless of success/failure -- failed runs are often the
        # most useful ones to have a copy of for debugging
        summary = (
            f"filename: {uploaded_file.name}\n"
            f"params: n_decoys={n_decoys}, max_perms={max_perms}, max_comps={max_comps}\n"
            f"status: {'OK' if result_bytes else 'FAILED'}\n"
            f"error: {run_error or ''}\n"
        )
        archive_session(input_bytes, uploaded_file.name, result_bytes, summary)

        if run_error:
            st.error(f"Processing failed: {run_error}")
        else:
            st.success("Done!")
            summary_df = pd.read_excel(io.BytesIO(result_bytes), sheet_name="Summary_Families")
            st.dataframe(summary_df, use_container_width=True)

            # Ion-level evidence: WHY each 2-10 AA call scored the way it did
            pep_df = pd.read_excel(io.BytesIO(result_bytes), sheet_name="Peptides_2to10AA")
            ion_df = pd.read_excel(io.BytesIO(result_bytes), sheet_name="Ion_Evidence")

            if len(pep_df) > 0:
                st.subheader("2-10 AA peptide calls — ion evidence")
                st.caption(
                    "For each candidate, the a/b/y/immonium ions checked against "
                    "the spectrum: theoretical vs. observed m/z, mass error (Δ Da "
                    "and Δ ppm), and whether each ion actually matched."
                )
                for _, prow in pep_df.iterrows():
                    fid = prow["Feature_ID"]
                    seq = prow["Top_candidate"]
                    tier = prow["Tier"]
                    composite = prow["Composite"]
                    if not seq or pd.isna(seq):
                        continue
                    with st.expander(f"{fid} — {seq}  ({tier}, composite={composite:.3f})"):
                        feature_ions = ion_df[ion_df["Feature_ID"] == fid].copy()
                        n_matched = int(feature_ions["matched"].sum())
                        n_total = len(feature_ions)
                        st.caption(f"{n_matched}/{n_total} theoretical ions matched within tolerance")
                        st.dataframe(
                            feature_ions[[
                                "ion_label", "ion_type", "theoretical_mz", "observed_mz",
                                "delta_Da", "delta_ppm", "observed_intensity", "matched",
                            ]],
                            use_container_width=True,
                            hide_index=True,
                        )

            st.download_button(
                label="Download results workbook",
                data=result_bytes,
                file_name="WineFeature_Results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )

st.divider()
st.caption(
    "Wine FeatureTriage — Step 3c of the Wine Peptidome series. "
    "[GitHub](https://github.com/314Olamda/WineFeatureTriage)"
)
