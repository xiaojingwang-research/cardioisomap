#!/usr/bin/env python3
"""Build a GENCODE v43 reference exon table for the Structure tab.

Parses the combined GTF (`gencode_v43_plus_novel_lr_may12.gtf.gz`), keeps only the
GENCODE reference exon rows (transcript_id starts with "ENST"), and restricts to the
gene symbols present in our dataset (gene_index.csv). Output mirrors exons.parquet so
the app can stack reference transcripts under the PacBio isoform models.

Columns: gene, transcript_id, transcript_name, chrom, strand, start, end
"""
import gzip
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROC = ROOT / "data" / "processed"
GTF = ROOT.parent / "gene_example" / "gencode_v43_plus_novel_lr_may12.gtf.gz"

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GN_RE = re.compile(r'gene_name "([^"]+)"')
TN_RE = re.compile(r'transcript_name "([^"]+)"')

# Only keep reference transcripts for genes we actually display.
wanted = set(pd.read_csv(PROC / "gene_index.csv")["gene"].dropna().astype(str))
print(f"target genes: {len(wanted):,}")

rows = []
opener = gzip.open
with opener(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.split("\t", 8)
        if len(parts) < 9 or parts[2] != "exon":
            continue
        attrs = parts[8]
        m_tx = TX_RE.search(attrs)
        if not m_tx or not m_tx.group(1).startswith("ENST"):
            continue  # skip PacBio/novel rows — those are in exons.parquet
        m_gn = GN_RE.search(attrs)
        gene = m_gn.group(1) if m_gn else None
        if gene not in wanted:
            continue
        m_tn = TN_RE.search(attrs)
        rows.append((
            gene,
            m_tx.group(1),
            m_tn.group(1) if m_tn else m_tx.group(1),
            parts[0],          # chrom
            parts[6],          # strand
            int(parts[3]),     # start
            int(parts[4]),     # end
        ))

df = pd.DataFrame(rows, columns=[
    "gene", "transcript_id", "transcript_name", "chrom", "strand", "start", "end"])
out = PROC / "gencode_exons.parquet"
df.to_parquet(out, index=False)
print(f"wrote {out}")
print(f"  exon rows:      {len(df):,}")
print(f"  transcripts:    {df['transcript_id'].nunique():,}")
print(f"  genes covered:  {df['gene'].nunique():,}")
