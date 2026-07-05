#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
streamlit_app.py — Web UI for WineFeatureTriage
════════════════════════════════════════════════

Upload a filled UHPLC-RP feature template, the tool runs
wine_feature_classifier.py against it, and you download the
processed results workbook. No local Python required for the user.

Deploy on Streamlit Community Cloud (share.streamlit.io):
  1. Push this file + wine_feature_classifier.py + linear_denovo.py +
     mass_utils.py + linear_peptide_scoring.py to the WineFeatureTriage repo
  2. Go to share.streamlit.io, sign in with GitHub, "New app"
  3. Point it at this repo, branch main, file streamlit_app.py
  4. Add a requirements.txt (see bottom of this file) to the repo root
  5. Deploy -- you get a permanent URL like
     https://314olamda-winefeaturetriage.streamlit.app

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
"""

import io
import tempfile
import os

import streamlit as st
import pandas as pd

from wine_feature_classifier import process_feature_table

st.set_page_config(page_title="Wine FeatureTriage", page_icon="🍷", layout="centered")

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
        with st.spinner("Classifying features and scoring 2-10 AA candidates..."):
            with tempfile.TemporaryDirectory() as tmpdir:
                input_path = os.path.join(tmpdir, "input.xlsx")
                output_path = os.path.join(tmpdir, "results.xlsx")

                with open(input_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())

                try:
                    process_feature_table(
                        input_path, output_path,
                        n_decoys=int(n_decoys),
                        max_perms_per_composition=int(max_perms),
                        max_compositions=int(max_comps),
                    )
                    with open(output_path, "rb") as f:
                        result_bytes = f.read()

                    st.success("Done!")

                    # Quick on-page preview of the summary sheet
                    summary_df = pd.read_excel(io.BytesIO(result_bytes), sheet_name="Summary_Families")
                    st.dataframe(summary_df, use_container_width=True)

                    st.download_button(
                        label="Download results workbook",
                        data=result_bytes,
                        file_name="WineFeature_Results.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                except Exception as e:
                    st.error(f"Processing failed: {e}")

st.divider()
st.caption(
    "Wine FeatureTriage — Step 3c of the Wine Peptidome series. "
    "[GitHub](https://github.com/314Olamda/WineFeatureTriage)"
)

# ══════════════════════════════════════════════════════════════════════════
# requirements.txt (create as a separate file in the repo root):
#
# streamlit
# pandas
# numpy
# scipy
# pyteomics
# openpyxl
# ══════════════════════════════════════════════════════════════════════════
