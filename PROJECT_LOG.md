# CardioIsoMap — Project Log

> Running memory of decisions, actions, and important notes. Newest entries at top.
> Goal: a gene-searchable web app (Streamlit) presenting a cardiac long-read isoform atlas.

---

## ⚠️ Important data notes (READ BEFORE ANALYSIS)

- **The `FL` column in the SQANTI classification is INCORRECT — do NOT use it for abundance.**
  - For expression/abundance use the **`total`** column (total count across samples) and the
    **per-sample count columns** (`sample1311`, `sample1518`, `sample1532`, `sample1535`,
    `sample1561`, `sample1662`).
  - Any "expression" shown in the app must be derived from `total` + per-sample counts, not `FL`.
- Sample→phenotype mapping (from `supplementary table 1.xlsx`), **3 phenotypes, 2 samples each**:
  | Phenotype | Samples |
  |-----------|---------|
  | NFH (non-failing control) | sample1561, sample1662 |
  | DCM (dilated cardiomyopathy) | sample1532, sample1535 |
  | IHD (ischemic heart disease) | sample1311, sample1518 |
  (The novel file's `DET` differential test = DCM vs NFD/NFH. Short-read split as `sr_mean_DCM`/`sr_mean_NFD`.)

---

## Data inventory (`data/raw/`)

| File | Rows | What it is |
|------|------|-----------|
| `all6_v2_classification_keep_decision.csv` | 62,672 isoforms | ALL isoforms — SQANTI classification + keep/exclude decision |
| `novel_lr_18637_annotated.csv` | 18,637 isoforms | NOVEL isoforms with extra verification (junction/short-read/DET) |
| `all6_v2_corrected.gtf` | 826,816 lines | Transcript structures (exon coordinates), PacBio corrected |

Source platform: **PacBio long-read (SQANTI3-style classification)**. IDs are `PB.x.y`.

### Key summary stats (computed 2026-06-05)
- **All isoforms:** 62,672 across **12,836 genes**.
- **Keep decision:** 48,232 `keep` · 14,440 `exclude`.
- **Structural category (all):**
  full-splice_match 29,595 · novel_not_in_catalog 13,656 · novel_in_catalog 12,045 ·
  incomplete-splice_match 6,076 · intergenic 398 · fusion 329 · antisense 288 · genic 286.
- **Novel file:** all 18,637 are `keep`.
  - `novel_feature`: no_novel_feature 6,703 · novel_region_and_junction 5,036 ·
    novel_junction_only 3,978 · novel_region_only 2,920.
  - `DET_significant` (differential transcript usage DCM vs NFD): TRUE 594 · FALSE 11,211 · NA 6,832.

### How the files relate
- `classification_keep_decision` = master table of every isoform; filter `keep_decision == "keep"`
  to get the curated set.
- `novel_lr_18637` = the novel subset that passed, enriched with verification evidence
  (junction support counts, short-read support, differential testing). Join key = `isoform` / `Name` (`PB.x.y`).
- `gtf` = structures for drawing exon/intron models; key = `transcript_id` (`PB.x.y`), `gene_id` (`PB.x`).

---

## Folder structure
```
CardioIsoMap/
├── PROJECT_LOG.md              # this file (running memory)
├── .gitignore                  # excludes raw/ + sashimi/ + bams/gtf; keeps parquet/csv
├── README.md                   # what it is, run locally, hosting, deploy
├── app/
│   ├── app.py                  # Streamlit app (CDN-aware Coverage tab)
│   └── requirements.txt        # app deps (for Streamlit Cloud deploy)
├── code/
│   ├── 01_build_tables.py      # builds isoforms/genes/gene_index
│   ├── 02_build_structures.py  # parses GTF -> exons table (PacBio isoforms)
│   ├── 03_render_sashimi.R     # batch Gviz coverage/sashimi PNGs
│   ├── 04_rebuild_manifest.py  # rebuild sashimi_manifest.csv from PNGs
│   ├── 05_build_gencode_exons.py # GENCODE v43 reference exon table
│   ├── 06_build_magnet_tpm.py   # MAGNET short-read isoform TPM (novel isoforms)
│   └── 07_build_magnet_known.py # MAGNET TPM for FSM/ISM via matched ENST
├── data/
│   ├── raw/                    # original uploaded files (do not edit)
│   │   ├── all6_v2_classification_keep_decision.csv
│   │   ├── novel_lr_18637_annotated.csv
│   │   ├── all6_v2_corrected.gtf
│   │   └── supplementary table 1.xlsx
│   └── processed/              # generated, app-ready tables
│       ├── isoforms.parquet    # 62,673 isoforms x 77 cols
│       ├── genes.parquet       # 12,836 genes
│       ├── gene_index.csv      # search autocomplete
│       ├── exons.parquet       # 764,144 exons (structure figures)
│       ├── gene_regions.csv    # per-gene plot window (chrom, padded start/end)
│       ├── sashimi/<gene>.png  # pre-rendered coverage/sashimi plots (~1.4 GB all genes)
│       └── sashimi_manifest.csv# gene -> png + render status
└── docs/                       # DESIGN.md, DATA_DICTIONARY.md, DEPLOY.md
```

## Coverage / Sashimi plots (Path B: pre-rendered images)
- Replicates `gene_example/01_plot_myom1.R` (Gviz): 3 phenotype coverage+sashimi tracks
  (NFH blue, DCM red, IHD green) over transcript models (known=red, novel `PB.*`=green).
- **Inputs** live in `../gene_example/`: `mapped_{NFH,DCM,IHD}.bam` (+ .bai) and
  `gencode_v43_plus_novel_lr_may12.gtf.gz`. **Not needed at app runtime** (only the PNGs are).
- `code/03_render_sashimi.R` imports the GTF once and loops over `gene_regions.csv`, writing one
  PNG per gene. Fusion-scale spans (>3 Mb) are skipped. Supports `--genes=`, `--limit=`,
  `--shard=i:n` (parallel), `--max-span=`, `--overwrite`.
- Rate ≈ 1.6 s/gene + ~52 s one-time GTF import; full set rendered with 5 parallel shards (~75 min).
- ⚠ **Storage:** ~1.4 GB for all genes — too big for a plain GitHub repo. For deploy, host images
  in a bucket (S3/GCS/Cloudflare R2) or a separate LFS/releases store and have the app load by path/URL.
- App: **Coverage / Sashimi** tab shows `data/processed/sashimi/<safe_gene>.png`
  (`safe_name = re.sub(r"[^A-Za-z0-9._-]","_", gene)`, matching the R script).

## How to run
```bash
# 1. build processed tables (only when raw data changes)
python3 code/01_build_tables.py
python3 code/02_build_structures.py
# 2. launch the app
python3 -m streamlit run app/app.py     # opens http://localhost:8501
```

## Processed-table schema (key columns)
- `isoforms.parquet`: identity (isoform, chrom, strand, length, exons, structural_category,
  associated_gene), coding info, `keep`/`is_novel` flags, raw per-sample counts (`sample*`),
  `cpm_*` per sample, `usage_*` per sample, per-phenotype means (`count_/cpm_/usage_{NFH,DCM,IHD}`),
  `usage_overall`, and novel-verification columns (junction/short-read/DET).
- `genes.parquet`: gene, chrom, strand, n_isoforms, n_keep, n_novel, dominant_isoform, n_DET_sig.
- `exons.parquet`: gene, transcript_id, chrom, strand, start, end (one row per exon, PacBio).
- `gencode_exons.parquet`: gene, transcript_id, transcript_name, chrom, strand, start, end —
  GENCODE v43 reference exon models (133k transcripts, 11,922 genes) for the Structure tab's
  optional reference overlay. Built by `code/05_build_gencode_exons.py` from
  `../gene_example/gencode_v43_plus_novel_lr_may12.gtf.gz` (ENST rows only, filtered to our genes).

---

## Changelog

### 2026-06-30 — Polish: novel count fix, MAGNET methods, viz tweaks
- **Fixed "Novel isoforms (total)" on the landing page.** It had used the `is_novel` flag (the
  curated `novel_lr_18637` set, all keepers → total == kept == 18,637, which is wrong for a
  "total"). Now defined as **ISM + NIC + NNC** structural categories: **31,777 total / 18,030
  kept**. The verified novel_lr count (18,637) is surfaced in the "kept" tooltip.
- **MAGNET methods note** added to the Expression-tab caption + README + DATA_DICTIONARY:
  "re-quantified with a Salmon decoy-aware pipeline against a merged reference (GENCODE v43 +
  18,637 novel long-read isoforms)."
- MAGNET box plots: show **every sample as a jittered dot** (`points="all"`, size 6) instead of
  outlier markers.
- Unified the structural-category palette: Structure figure now reuses the overview
  `PIE_CAT_COLORS` (single source of truth). Overview phenotype colors → `#ef767a/#456990/#49beaa`.
- Larger fonts (16pt base / 15pt ticks) on landing, overview, and expression charts; landing
  category chart switched to a before→after **overlay** (kept nested in all-isoforms).
- Structure figure made smaller (fits when expanded). README features/tab list refreshed; docs
  note that the Sashimi build steps are now optional/legacy.

### 2026-06-29 (cont.) — Tab reorder, drop Sashimi, side-by-side Expression
- Per-gene tabs reordered to **Overview · Structure · Isoforms table · Expression · Novel
  evidence**; the **Coverage / Sashimi tab was removed** from the app. (The render script
  `03_render_sashimi.R`, the PNGs and manifest still exist on disk, just not surfaced in the UI.)
  → This also removes the need for the jsDelivr/figures-repo CDN at deploy time; the app is now
  a single small repo again. (`IMG_BASE`/sashimi helpers left in code, unused.)
- Expression tab is now **two columns**: long-read isoform-usage (left) vs MAGNET short-read
  TPM box plots (right), instead of stacked top/bottom.

### 2026-06-29 (cont.) — Richer Isoforms-table columns
- Replaced the per-phenotype CPM columns in the Isoforms-table tab with informative SQANTI3
  classification fields: `associated_transcript`, `ref_exons`, `ORF_length`, `predicted_NMD`,
  `all_canonical`, `RTS_stage`, `within_CAGE_peak`, `polyA_motif_found`, `pct_junct_supported`,
  `ML_filter` (kept identity/structure/coding/abundance basics). Added `st.column_config`
  header tooltips + formatting (junc % , usage 3dp). All fields already present in
  isoforms.parquet — no rebuild. Per-sample counts/CPM/DET evidence remain in the CSV download
  and Novel-evidence tab.

### 2026-06-29 (cont.) — Replaced LR heatmap with isoform-usage stacked bars
- User disliked the long-read expression heatmap. Replaced it with **stacked horizontal usage
  bars** (`usage_stack_figure`): one 100%-bar per phenotype (or per sample); segments = isoforms
  sized by share of the gene's reads, so isoform switching across NFH/DCM/IHD is visible at a
  glance. Top-N isoforms colored distinctly, the rest lumped into grey "Other"; Per-phenotype /
  Per-sample toggle + N slider. The original CPM/count/usage heatmap is retained inside a
  collapsible "Show raw expression matrix" expander (info kept, but out of the way). Verified
  on MYH7 (dominant FSM ~90% across all groups).

### 2026-06-29 (cont.) — MAGNET box plots in Expression tab
- Added MAGnet short-read isoform expression to the Expression tab as **box plots** (per
  isoform, split by disease group), above the existing long-read heatmap.
- `code/06_build_magnet_tpm.py` → `data/processed/magnet_isoform_tpm.parquet` (18,637 novel
  PB isoforms × 366 SRA samples, float32, 21 MB) + `magnet_samples.csv` (Run→etiology).
  Source: `Magnet_DTE/salmon_D_may12_transcript_tpm_matrix.csv` (salmon index = GENCODE v43 +
  novel-LR, so only the 18,637 novel isoforms have a PB row — known isoforms are absent).
  Groups from `SraRunTable.xls` `etiology`: DCM 166 · Non-Failing 162 · HCM 27 · PPCM 6.
- App: Expression tab box plot with disease-group multiselect (default Non-Failing + DCM),
  max-isoforms slider, log-y toggle; isoforms ranked by mean TPM. Verified live on MYH7.
- **Extended MAGNET to KNOWN isoforms** (`code/07_build_magnet_known.py`): FSM/ISM isoforms
  have no PB salmon row, but SQANTI's `associated_transcript` gives their matched ENST, and the
  salmon index contains every GENCODE ENST. Pulled those ENST rows →
  `magnet_enst_tpm.parquet` (18,396 unique ENSTs × 366, float32, 31 MB) +
  `isoform_ref_map.csv` (35,671 FSM/ISM isoforms → ref ENST, 1.6 MB). The box plot now unions
  novel (own row) + FSM/ISM (matched-ENST TPM; exact for FSM, approximate for ISM), with the
  source shown on hover. MYH7 coverage 57→71; NPPA now plots its FSM isoform and shows the
  expected ~10× DCM up-regulation.

### 2026-06-29 — Overview chart + GENCODE reference in Structure tab
- **Overview "Isoforms by structural category" chart** redrawn: replaced the stretched
  `st.bar_chart` with a compact **horizontal grouped Plotly bar** (All vs After-filtering,
  sorted, fixed height) — the old vertical chart looked too wide.
- **Structure tab — optional GENCODE v43 reference overlay.** Added
  `code/05_build_gencode_exons.py` → `data/processed/gencode_exons.parquet` (ENST exon models
  from the combined GTF, filtered to our gene symbols; 1.08M exons / 133k transcripts /
  11,922 genes / 8.3 MB). A "Include GENCODE v43 reference transcripts" toggle appends those
  models (grey, labeled by transcript_name e.g. `MYH7-201`) below the PacBio isoforms, with a
  "GENCODE v43" legend entry. Verified MYH7 reference aligns with its FSM isoform.

### 2026-06-08 — Deployment prep (GitHub + Streamlit Cloud + jsDelivr CDN)
- Decision: host the app on **Streamlit Community Cloud** (small app repo) and the **1.3 GB of
  sashimi PNGs in a separate public repo** served free via the **jsDelivr CDN** (no egress fees).
- Rebuilt the sashimi manifest: wrote `code/04_rebuild_manifest.py` (the sharded render's merge had
  left only 79 rows). It scans `data/processed/sashimi/*.png`, maps stems→genes via
  `gene_regions.csv` (same `safe_name` sanitisation), and writes `sashimi_manifest.csv` →
  **12,677 rows, 0 unmatched**; MYOM2/MYH7/TTN/NPPA all present (~159 large fusion-span genes skipped).
- Made the **Coverage tab CDN-aware** in `app/app.py`:
  - `IMG_BASE` from `st.secrets["IMG_BASE"]` or env `CARDIOISOMAP_IMG_BASE`.
  - cached `available_sashimi()` reads `sashimi_manifest.csv`.
  - if `IMG_BASE` set → `st.image(f"{IMG_BASE}/{safe_name(gene)}.png")` + "download full-res" link;
    else falls back to the local `sashimi/<gene>.png` (dev); else "no plot (skipped)" message.
- Added `.gitignore` (excludes `data/raw/`, `data/processed/sashimi/`, `*.bam`, `*.gtf*`,
  `__pycache__/`, `.DS_Store`; **keeps** `data/processed/*.parquet` + `*.csv`, ~11 MB runtime data).
- Added `README.md` and `docs/DEPLOY.md` (exact two-repo `git` push commands + jsDelivr URL +
  Streamlit Cloud `IMG_BASE` secret). `gh` is not installed locally → using plain `git` + GitHub web.
- Verified locally in BOTH modes via headless browser: local-file image render, and CDN mode
  (constructs `https://cdn.jsdelivr.net/gh/<USER>/cardioisomap-figures@main/sashimi/MYOM2.png`).
- **Still user-driven (not done):** create the two GitHub repos, push (~1.3 GB one-time image
  upload), and click Deploy on share.streamlit.io with the `IMG_BASE` secret.

### 2026-06-05
- Organized folder: created `app/`, `code/`, `data/raw/`, `docs/`; moved the 4 uploaded files into `data/raw/`.
- Inspected all files; recorded schema + summary stats above.
- Recorded critical note: **use `total`/per-sample counts, not `FL`, for abundance.**
- Got sample→phenotype mapping from `supplementary table 1.xlsx` (NFH/DCM/IHD, 2 each).
- Web app design finalized (see `docs/DESIGN.md`): Streamlit, gene-searchable.
- Built `code/01_build_tables.py` + `code/02_build_structures.py`; generated 4 processed tables.
  - CPM sums verified = 1e6/sample; 48,232 keep; 18,637 novel; 764,144 exons (all mapped to a gene).
- Built `app/app.py` (5 tabs: Overview, Expression, Isoforms table, Structure, Novel evidence) +
  `requirements.txt`. Decisions applied: keep-only default w/ toggle; CPM + usage; Structure = drawn figure.
- Verified end-to-end in a headless browser: landing page, gene search (MYH7/TTN), all tabs,
  Plotly expression heatmap, and matplotlib structure figure all render.
- Expanded the landing page into a **dataset-summary** view: genes total (12,836) vs after
  filtering (11,760); isoforms total (62,673) vs kept (48,232); excluded (14,441); novel
  (18,637, all kept); % kept (77.0%); plus a before/after structural-category breakdown table+chart.

### 2026-06-05 (cont.) — Coverage/Sashimi feature
- User requested the `gene_example/MYOM2.pdf`-style plot for every gene. Chose **Path B**
  (pre-rendered images), **all genes**, **no BAMs at runtime**.
- Built `data/processed/gene_regions.csv` (per-gene padded window from exons.parquet).
- Wrote `code/03_render_sashimi.R` (Gviz, GTF imported once, shardable). Verified MYOM2 output
  matches the reference PDF; TTN/MYH7/NPPA render correctly.
- Added **Coverage / Sashimi** tab to `app/app.py` (shows `sashimi/<gene>.png` + PNG download).
- Launched full render: 5 parallel shards (~12,834 genes, ~75 min, ~1.4 GB total).

### Next steps / not yet done
- Deploy (user-driven, commands in `docs/DEPLOY.md`): create `cardioisomap` (app) +
  `cardioisomap-figures` (PNGs) repos on GitHub, push both, then New app on share.streamlit.io
  with the `IMG_BASE` jsDelivr secret. `data/processed/*.parquet|*.csv` (~11 MB) are committed.
- Optional polish: gene full names, alias search, protein-domain track, caching tuning.

---

## Open TODOs
- [x] Confirm sample→phenotype mapping (NFH/DCM/IHD — from supplementary table 1).
- [x] Decide main value to display: raw counts + CPM + usage proportion.
- [x] Build preprocessing scripts producing app-ready tables in `data/processed/`.
- [x] Scaffold Streamlit app with gene autocomplete search.
- [x] Structure tab = drawn exon/intron figure (per user).
- [x] Coverage/Sashimi: pre-render all genes + CDN-aware Coverage tab + rebuilt manifest.
- [x] Deploy prep: `.gitignore`, `README.md`, `docs/DEPLOY.md`, two-repo plan.
- [ ] Deploy: create+push GitHub repos, set `IMG_BASE` secret, launch on Streamlit Cloud.
- [ ] Optional: gene full names, alias/Ensembl-ID search, protein-domain track.
