#!/usr/bin/env python3
"""
CardioIsoMap — Step 2: parse the corrected GTF into a per-exon table used to
draw transcript-structure figures in the app.

Output (data/processed/):
  - exons.parquet   one row per exon: gene, transcript_id, chrom, strand, start, end
                    (start/end are 1-based GTF coordinates)
"""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
GTF = RAW / "all6_v2_corrected.gtf"

TID_RE = re.compile(r'transcript_id "([^"]+)"')


def main() -> None:
    print(f"Reading {GTF.name} ...")
    cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    gtf = pd.read_csv(GTF, sep="\t", header=None, names=cols, comment="#", low_memory=False)
    exons = gtf[gtf["feature"] == "exon"].copy()
    print(f"  {len(exons):,} exon rows")

    exons["transcript_id"] = exons["attr"].str.extract(TID_RE)
    exons = exons.dropna(subset=["transcript_id"])

    # Map transcript_id -> gene symbol via the processed isoform table.
    iso = pd.read_parquet(OUT / "isoforms.parquet", columns=["isoform", "associated_gene"])
    tx2gene = dict(zip(iso["isoform"], iso["associated_gene"]))
    exons["gene"] = exons["transcript_id"].map(tx2gene)

    out = exons[["gene", "transcript_id", "chrom", "strand", "start", "end"]].copy()
    out["start"] = out["start"].astype("int64")
    out["end"] = out["end"].astype("int64")
    out = out.sort_values(["gene", "transcript_id", "start"]).reset_index(drop=True)
    out.to_parquet(OUT / "exons.parquet", index=False)

    n_tx = out["transcript_id"].nunique()
    n_mapped = out["gene"].notna().sum()
    print(f"  wrote exons.parquet  ({len(out):,} exons, {n_tx:,} transcripts)")
    print(f"  exons mapped to a gene symbol: {n_mapped:,} / {len(out):,}")


if __name__ == "__main__":
    main()
