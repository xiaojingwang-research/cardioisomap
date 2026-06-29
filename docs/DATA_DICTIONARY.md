# Data Dictionary

Key columns used by the app. (SQANTI3 classification schema + custom keep/verification columns.)

## Identity / structure
| Column | Meaning |
|--------|---------|
| `isoform` / `Name` | Transcript ID, `PB.x.y` (join key across files & GTF `transcript_id`) |
| `chrom`, `strand` | Genomic chromosome and strand |
| `length`, `exons` | Transcript length (bp), exon count |
| `structural_category` | SQANTI category: full-splice_match, novel_in_catalog, novel_not_in_catalog, ISM, fusion, etc. |
| `subcategory` | Finer classification (e.g. at_least_one_novel_splicesite) |
| `associated_gene` | Gene symbol/ID this isoform maps to (**primary search key**) |
| `associated_transcript` | Matched reference transcript (or `novel`) |
| `ref_length`, `ref_exons` | Reference transcript length / exon count |

## Coding / quality
| Column | Meaning |
|--------|---------|
| `coding` | coding / non_coding |
| `ORF_length`, `CDS_length`, `CDS_start/end`, `CDS_genomic_start/end` | ORF/CDS info |
| `predicted_NMD` | Predicted nonsense-mediated decay target |
| `all_canonical` | All splice junctions canonical |
| `RTS_stage` | Possible reverse-transcription template-switching artifact |
| `within_CAGE_peak`, `dist_to_CAGE_peak` | 5' (TSS) support |
| `within_polyA_site`, `polyA_motif_found` | 3' (polyA) support |

## Abundance  ⚠️
| Column | Meaning |
|--------|---------|
| **`total`** | **Total count across samples — USE THIS for abundance** |
| **`sample1311 … sample1662`** | **Per-sample counts — USE THESE** |
| ~~`FL`~~ | ⚠️ **INCORRECT — do not use** |
| `iso_exp`, `gene_exp`, `ratio_exp` | SQANTI-provided (mostly NA here) |
| `n_samples_detected` | # samples isoform detected in |

## Filtering decision
| Column | Meaning |
|--------|---------|
| `ML_filter` | ML classifier call: `Isoform` vs `Artifact` |
| `keep_path` | Reason string (e.g. `ML_Isoform`, `Artifact_rescued`, `EXCLUDE_artifact_no_rescue`) |
| `keep_decision` | Final `keep` / `exclude` |

## Novel-isoform verification (novel file only)
| Column | Meaning |
|--------|---------|
| `novel_isoform` | Flagged novel |
| `n_junctions`, `n_junct_supported`, `n_junct_weak`, `n_junct_absent`, `pct_junct_supported` | Splice-junction support from data |
| `all_junct_supported`, `all_junct_not_absent` | Junction support flags |
| `sr_mean_count`, `sr_n_detected`, `sr_pct_detected`, `sr_n_count10` | Short-read support |
| `sr_mean_DCM`, `sr_mean_NFD`, `sr_expressed` | Short-read mean by condition |
| `DET_log2FC`, `DET_padj`, `DET_significant`, `DET_direction` | Differential transcript usage/expression (DCM vs NFD) |
| `has_novel_region`, `has_novel_junction`, `novel_feature` | Type of novelty |
