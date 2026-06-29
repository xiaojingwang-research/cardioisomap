# CardioIsoMap — Web App Design

**Stack:** Streamlit (dev locally, deploy free on Streamlit Community Cloud via GitHub).
**Core idea:** user searches a **gene** → sees that gene's isoforms, abundance, and evidence.

## Layout
- Sidebar: gene autocomplete search + filters (structural category, keep-only, novel-only) + download + about.
- Main panel: landing page when empty; per-gene view when a gene is selected.

## Per-gene view (tabs)
1. **Overview** — gene identity (symbol, location, # isoforms, dominant isoform), external links.
2. **Isoform expression** — isoform × sample/condition heatmap or bars.
   Toggle: raw count ↔ usage proportion. (Abundance from `total` + per-sample counts — NOT `FL`.)
3. **Coverage / Sashimi** — pre-rendered Gviz PNG: per-phenotype coverage + splice-junction arcs
   over transcript models (known=red, novel `PB.*`=green). Built in batch by `code/03_render_sashimi.R`;
   only the PNGs are needed at runtime (no BAMs on the server).
4. **Isoforms table** — per-isoform: ID, category, length/exons, coding, keep_decision, counts; CSV download.
4. **Structure** — a **drawn figure** (matplotlib) of exon/intron transcript models per isoform
   (exons = boxes, introns = lines, strand-aware), with PNG download. NOT a table.
5. **Novel evidence** *(for novel isoforms)* — junction support, short-read support, DET stats.

## Per-gene information to show (agreed list)
- Identity: symbol, name, IDs, chrom:location, strand, biotype, external links.
- Summary: # isoforms, dominant isoform, # novel, # significant DET.
- Per-isoform table (the core).
- Expression/usage visualization.
- Structure (if available).
- Provenance/citation footer.

## Resolved decisions
- Phenotype axis = **NFH / DCM / IHD** (2 samples each; see PROJECT_LOG).
- Main value = **raw count + CPM + usage proportion** (toggle in Expression tab).
- Default = **keep-only** isoforms with sidebar toggle to show excluded.
- Structure = drawn figure (see above).

## Still open
- Search currently by gene **symbol** only (alias/Ensembl-ID search = future).
- No gene full names yet (external links cover this).

## Status: IMPLEMENTED — see `app/app.py`, verified locally 2026-06-05.
