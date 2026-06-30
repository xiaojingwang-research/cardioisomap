# 🫀 CardioIsoMap

A gene-searchable atlas of cardiac transcript isoforms from PacBio long-read sequencing
(SQANTI3) across human heart phenotypes — **NFH** (non-failing), **DCM** (dilated
cardiomyopathy), and **IHD** (ischemic heart disease).

Search a gene and explore its isoforms, transcript structures, long-read abundance/usage,
independent short-read (MAGNET) expression, and novel-isoform evidence.

## Features
- **Dataset summary** landing page (genes / isoforms before & after filtering, category breakdown).
- Per-gene tabs: **Overview · Structure · Isoforms table · Expression · Novel evidence**.
- Long-read abundance from per-sample read counts and `total` (the SQANTI `FL` column is **not** used).
- **Expression tab**: long-read isoform-usage (left) vs **MAGNET short-read TPM box plots** (right).

### Short-read expression (MAGNET)
The independent **MAGnet** cohort (366 human-heart RNA-seq samples) was **re-quantified with a
Salmon decoy-aware pipeline against a merged reference (GENCODE v43 + the 18,637 novel long-read
isoforms)**. Per-isoform TPM is shown as box plots split by disease group (Non-Failing / DCM /
HCM / PPCM). Novel isoforms use their own salmon row; known FSM/ISM isoforms use the TPM of their
matched GENCODE transcript (`associated_transcript`) — exact for FSM, approximate for ISM.

## Run locally
```bash
pip install -r app/requirements.txt
python3 -m streamlit run app/app.py        # http://localhost:8501
```

### Rebuild the processed data (only if raw inputs change)
Raw inputs live in `data/raw/` (not committed). Then:
```bash
python3 code/01_build_tables.py         # isoforms/genes/gene_index tables
python3 code/02_build_structures.py     # exon table (PacBio) for structure figures
python3 code/05_build_gencode_exons.py  # GENCODE v43 reference exon models (structure overlay)
python3 code/06_build_magnet_tpm.py     # MAGNET TPM for novel isoforms (box plots)
python3 code/07_build_magnet_known.py   # MAGNET TPM for FSM/ISM via matched ENST
```
(`code/03_render_sashimi.R` + `04_rebuild_manifest.py` build the coverage/sashimi PNGs; the
Sashimi tab has since been removed from the app, so those steps are now optional/legacy.)

## Coverage/sashimi images (hosted on a CDN)
The ~12,677 per-gene PNGs (~1.3 GB) are **not** in this repo. They live in a separate public
repo and are served free via **jsDelivr**. The app reads the image URL from a config value:

- Streamlit secret `IMG_BASE`, or
- env var `CARDIOISOMAP_IMG_BASE`

Example:
```
IMG_BASE = "https://cdn.jsdelivr.net/gh/<USER>/cardioisomap-figures@main/sashimi"
```
If unset, the app falls back to local files in `data/processed/sashimi/`.

## Deploy (Streamlit Community Cloud)
1. Push this repo to GitHub (the 1.3 GB images are git-ignored).
2. Push the `data/processed/sashimi/` PNGs to a **separate public repo** `cardioisomap-figures`.
3. On https://share.streamlit.io → New app → this repo, main file `app/app.py`.
4. App settings → Secrets → set `IMG_BASE` to the jsDelivr URL above.

## Layout
```
app/         Streamlit app (app.py) + requirements.txt
code/        01..04 preprocessing / rendering scripts
data/raw/    original inputs (local only, git-ignored)
data/processed/  parquet + csv tables (committed) + sashimi/ PNGs (git-ignored)
docs/        design notes + data dictionary
PROJECT_LOG.md   running project log
```

## Data note
Abundance uses the `total` column + per-sample counts (`sample1311/1518/1532/1535/1561/1662`).
The SQANTI `FL` column is **not** used for abundance.
