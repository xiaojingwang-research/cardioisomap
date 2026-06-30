#!/usr/bin/env python3
"""Build MAGNET short-read isoform TPM table for the Expression tab box plots.

Source: MAGnet salmon quantification (`Magnet_DTE/salmon_D_may12_transcript_tpm_matrix.csv`,
transcripts x 366 SRA samples) — the salmon index was GENCODE v43 + the 18,637 novel LR
isoforms, so only those novel PacBio isoforms (PB.x.y) carry a row here. We keep those rows
(float32, wide) plus a small Run->etiology map from SraRunTable.xls.

Outputs (data/processed/):
  magnet_isoform_tpm.parquet : transcript_id + 366 sample TPM columns (float32)
  magnet_samples.csv         : Run, etiology  (disease group per sample)
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROC = ROOT / "data" / "processed"
MAG = ROOT.parent / "Magnet_DTE"
TPM_CSV = MAG / "salmon_D_may12_transcript_tpm_matrix.csv"
SRA = MAG / "SraRunTable.xls"

# ---- 1. sample -> etiology -------------------------------------------------
meta = pd.read_excel(SRA)
samp = meta[["Run", "etiology"]].dropna(subset=["Run"]).copy()
# short, plot-friendly group labels
LABEL = {
    "Non-Failing Donor": "Non-Failing",
    "Dilated cardiomyopathy (DCM)": "DCM",
    "Hypertrophic cardiomyopathy (HCM)": "HCM",
    "Peripartum cardiomyopathy (PPCM)": "PPCM",
}
samp["etiology"] = samp["etiology"].map(lambda x: LABEL.get(x, x))
samp.to_csv(PROC / "magnet_samples.csv", index=False)
print(f"samples: {len(samp)}  groups: {samp['etiology'].value_counts().to_dict()}")

# ---- 2. PB isoform TPM rows (stream to keep memory low) --------------------
pb_chunks = []
reader = pd.read_csv(TPM_CSV, index_col=0, chunksize=20000)
for chunk in reader:
    pb = chunk[chunk.index.astype(str).str.startswith("PB.")]
    if not pb.empty:
        pb_chunks.append(pb.astype("float32"))
tpm = pd.concat(pb_chunks)
tpm.index.name = "transcript_id"
tpm = tpm.reset_index()
out = PROC / "magnet_isoform_tpm.parquet"
tpm.to_parquet(out, index=False)
print(f"wrote {out}")
print(f"  isoforms (PB rows): {len(tpm):,}")
print(f"  sample columns:     {tpm.shape[1] - 1:,}")
