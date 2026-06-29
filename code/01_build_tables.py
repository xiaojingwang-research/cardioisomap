#!/usr/bin/env python3
"""
CardioIsoMap — Step 1: build app-ready isoform/gene tables.

Reads the SQANTI classification + novel-verification files and produces compact
parquet tables for the Streamlit app.

IMPORTANT DATA RULE
-------------------
The `FL` column in the SQANTI classification is INCORRECT and must NOT be used
for abundance. Abundance is taken from the `total` column and the per-sample
count columns (sample1311, sample1518, sample1532, sample1535, sample1561,
sample1662).

Outputs (data/processed/):
  - isoforms.parquet   per-isoform metadata + counts + CPM + usage + novel evidence
  - genes.parquet      per-gene summary
  - gene_index.csv     gene symbol + n_isoforms (drives the search box)
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

CLASS_CSV = RAW / "all6_v2_classification_keep_decision.csv"
NOVEL_CSV = RAW / "novel_lr_18637_annotated.csv"

# Per-sample count columns (USE THESE, not FL) and their phenotype.
SAMPLES = ["sample1311", "sample1518", "sample1532", "sample1535", "sample1561", "sample1662"]
SAMPLE_PHENO = {
    "sample1561": "NFH", "sample1662": "NFH",   # non-failing control
    "sample1532": "DCM", "sample1535": "DCM",   # dilated cardiomyopathy
    "sample1311": "IHD", "sample1518": "IHD",   # ischemic heart disease
}
PHENOTYPES = ["NFH", "DCM", "IHD"]

# Columns carried through from the classification table.
META_COLS = [
    "isoform", "chrom", "strand", "length", "exons", "structural_category",
    "subcategory", "associated_gene", "associated_transcript", "ref_length",
    "ref_exons", "coding", "ORF_length", "CDS_length", "predicted_NMD",
    "all_canonical", "RTS_stage", "within_CAGE_peak", "within_polyA_site",
    "polyA_motif_found", "ML_filter", "keep_path", "keep_decision",
]

# Novel-verification columns to merge in (keyed by isoform).
NOVEL_COLS = [
    "novel_isoform", "n_junctions", "n_junct_supported", "n_junct_weak",
    "n_junct_absent", "pct_junct_supported", "all_junct_supported",
    "all_junct_not_absent", "sr_mean_count", "sr_n_detected", "sr_pct_detected",
    "sr_n_count10", "sr_mean_DCM", "sr_mean_NFD", "sr_expressed",
    "DET_log2FC", "DET_padj", "DET_significant", "DET_direction",
    "has_novel_region", "has_novel_junction", "novel_feature",
]


def main() -> None:
    print(f"Reading {CLASS_CSV.name} ...")
    df = pd.read_csv(CLASS_CSV, low_memory=False)
    df = df[df["isoform"].notna()].copy()  # drop any trailing/blank rows
    print(f"  {len(df):,} isoforms x {df['associated_gene'].nunique():,} genes")

    # --- Abundance: raw counts (NA -> 0) ---------------------------------
    for s in SAMPLES:
        df[s] = pd.to_numeric(df[s], errors="coerce").fillna(0.0)
    df["total"] = pd.to_numeric(df["total"], errors="coerce").fillna(0.0)
    # Use the sum of per-sample counts as the authoritative total (consistent CPM/usage).
    df["total_count"] = df[SAMPLES].sum(axis=1)

    # --- CPM per sample (count / library size * 1e6) ---------------------
    lib_sizes = df[SAMPLES].sum(axis=0)
    print("Library sizes (sum of counts per sample):")
    for s in SAMPLES:
        print(f"  {s} ({SAMPLE_PHENO[s]}): {lib_sizes[s]:,.0f}")
    cpm_cols = []
    for s in SAMPLES:
        col = f"cpm_{s}"
        df[col] = df[s] / lib_sizes[s] * 1e6 if lib_sizes[s] else 0.0
        cpm_cols.append(col)
    df["cpm_total"] = df[cpm_cols].sum(axis=1)

    # --- Usage proportion: isoform count / gene-total count --------------
    # Overall usage (from summed counts).
    gene_total = df.groupby("associated_gene")["total_count"].transform("sum")
    df["usage_overall"] = np.where(gene_total > 0, df["total_count"] / gene_total, 0.0)
    # Per-sample usage.
    usage_cols = []
    for s in SAMPLES:
        gtot = df.groupby("associated_gene")[s].transform("sum")
        col = f"usage_{s}"
        df[col] = np.where(gtot > 0, df[s] / gtot, 0.0)
        usage_cols.append(col)

    # --- Phenotype aggregation (mean of the 2 samples per phenotype) -----
    pheno_count_cols, pheno_cpm_cols, pheno_usage_cols = [], [], []
    for ph in PHENOTYPES:
        members = [s for s in SAMPLES if SAMPLE_PHENO[s] == ph]
        df[f"count_{ph}"] = df[members].mean(axis=1)
        df[f"cpm_{ph}"] = df[[f"cpm_{s}" for s in members]].mean(axis=1)
        df[f"usage_{ph}"] = df[[f"usage_{s}" for s in members]].mean(axis=1)
        pheno_count_cols.append(f"count_{ph}")
        pheno_cpm_cols.append(f"cpm_{ph}")
        pheno_usage_cols.append(f"usage_{ph}")

    # --- Merge novel verification ---------------------------------------
    print(f"Reading {NOVEL_CSV.name} ...")
    nv = pd.read_csv(NOVEL_CSV, low_memory=False)
    keep_novel = ["isoform"] + [c for c in NOVEL_COLS if c in nv.columns]
    nv = nv[keep_novel].drop_duplicates("isoform")
    df = df.merge(nv, on="isoform", how="left")
    df["is_novel"] = df["isoform"].isin(set(nv["isoform"]))

    df["keep"] = df["keep_decision"].astype(str).str.lower().eq("keep")

    # --- Assemble isoform output ----------------------------------------
    out_cols = (
        META_COLS + ["total_count", "is_novel", "keep"]
        + SAMPLES + cpm_cols + ["cpm_total"]
        + ["usage_overall"] + usage_cols
        + pheno_count_cols + pheno_cpm_cols + pheno_usage_cols
        + [c for c in NOVEL_COLS if c in df.columns]
    )
    out_cols = [c for c in dict.fromkeys(out_cols) if c in df.columns]
    iso = df[out_cols].copy()
    iso.to_parquet(OUT / "isoforms.parquet", index=False)
    print(f"  wrote isoforms.parquet  ({len(iso):,} rows, {len(out_cols)} cols)")

    # --- Gene-level summary ---------------------------------------------
    g = df.groupby("associated_gene")
    genes = pd.DataFrame({
        "gene": g.size().index,
        "chrom": g["chrom"].agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0]).values,
        "strand": g["strand"].agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0]).values,
        "n_isoforms": g.size().values,
        "n_keep": g["keep"].sum().values,
        "n_novel": g["is_novel"].sum().values,
        "total_count": g["total_count"].sum().values,
    })
    # genomic span from GTF-independent fields isn't available here; leave to structures step.
    # dominant isoform = highest total_count per gene
    dom = (df.sort_values("total_count", ascending=False)
             .drop_duplicates("associated_gene")[["associated_gene", "isoform"]]
             .rename(columns={"associated_gene": "gene", "isoform": "dominant_isoform"}))
    genes = genes.merge(dom, on="gene", how="left")
    # number of significant DET isoforms per gene
    if "DET_significant" in df.columns:
        sig = (df.assign(_sig=df["DET_significant"].astype(str).str.upper().eq("TRUE"))
                 .groupby("associated_gene")["_sig"].sum()
                 .rename("n_DET_sig").reset_index()
                 .rename(columns={"associated_gene": "gene"}))
        genes = genes.merge(sig, on="gene", how="left")
        genes["n_DET_sig"] = genes["n_DET_sig"].fillna(0).astype(int)
    genes = genes.sort_values("gene").reset_index(drop=True)
    genes.to_parquet(OUT / "genes.parquet", index=False)
    print(f"  wrote genes.parquet     ({len(genes):,} genes)")

    # --- Gene index for the search box (keep-only count too) ------------
    idx = genes[["gene", "n_isoforms", "n_keep", "n_novel"]].copy()
    idx.to_csv(OUT / "gene_index.csv", index=False)
    print(f"  wrote gene_index.csv    ({len(idx):,} genes)")

    # --- Sanity checks ---------------------------------------------------
    print("\nSanity checks:")
    cpm_sums = iso[cpm_cols].sum(axis=0)
    print("  CPM column sums (should be ~1e6 each):")
    for c in cpm_cols:
        print(f"    {c}: {cpm_sums[c]:,.0f}")
    print(f"  keep isoforms: {iso['keep'].sum():,} / {len(iso):,}")
    print(f"  novel isoforms: {iso['is_novel'].sum():,}")


if __name__ == "__main__":
    main()
