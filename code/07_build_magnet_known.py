#!/usr/bin/env python3
"""Extend MAGNET coverage to KNOWN (FSM/ISM) isoforms via their matched GENCODE transcript.

Novel isoforms (NIC/NNC...) have their own PB.x.y row in the salmon matrix (built by
06_build_magnet_tpm.py). Known isoforms do not — but SQANTI records the reference transcript
they match in `associated_transcript` (an ENST). The salmon index also contains every GENCODE
ENST, so we can pull those ENST rows and use them as the MAGNET expression proxy for FSM/ISM
isoforms.

  FSM = exact full-length match  -> ENST TPM is a faithful proxy.
  ISM = partial (5'/3' truncated) -> ENST TPM is an approximation (full-transcript abundance).

Outputs (data/processed/):
  magnet_enst_tpm.parquet : ref_enst + 366 sample TPM columns (unique ENSTs, float32)
  isoform_ref_map.csv     : isoform, ref_enst, structural_category   (FSM/ISM only)
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
MAG = ROOT.parent / "Magnet_DTE"
TPM_CSV = MAG / "salmon_D_may12_transcript_tpm_matrix.csv"
CLASS_CSV = RAW / "all6_v2_classification_keep_decision.csv"

# ---- 1. known isoform -> reference ENST ------------------------------------
cls = pd.read_csv(CLASS_CSV,
                  usecols=["isoform", "structural_category", "associated_transcript"])
known = cls[cls["structural_category"].isin(
    ["full-splice_match", "incomplete-splice_match"])].copy()
known = known[known["associated_transcript"].fillna("").str.startswith("ENST")]
known = known.rename(columns={"associated_transcript": "ref_enst"})
known[["isoform", "ref_enst", "structural_category"]].to_csv(
    PROC / "isoform_ref_map.csv", index=False)
needed = set(known["ref_enst"])
print(f"known FSM/ISM isoforms: {len(known):,}  unique ref ENSTs: {len(needed):,}")

# ---- 2. pull those ENST rows from the salmon TPM matrix ---------------------
keep = []
reader = pd.read_csv(TPM_CSV, index_col=0, chunksize=20000)
for chunk in reader:
    # salmon ids look like "ENST00000456328.2|ENSG...|...": match on the ENST part
    enst = chunk.index.to_series().astype(str).str.split("|").str[0]
    hit = chunk[enst.isin(needed).values].copy()
    if not hit.empty:
        hit.index = enst[enst.isin(needed).values].values
        keep.append(hit.astype("float32"))

enst_tpm = pd.concat(keep)
enst_tpm = enst_tpm[~enst_tpm.index.duplicated()]   # one row per ENST
enst_tpm.index.name = "ref_enst"
enst_tpm = enst_tpm.reset_index()
out = PROC / "magnet_enst_tpm.parquet"
enst_tpm.to_parquet(out, index=False)
print(f"wrote {out}")
print(f"  ENST rows found:  {len(enst_tpm):,} / {len(needed):,} needed")
print(f"  sample columns:   {enst_tpm.shape[1] - 1:,}")
missing = len(needed) - len(enst_tpm)
if missing:
    print(f"  (note: {missing:,} ref ENSTs not in salmon matrix — dropped)")
