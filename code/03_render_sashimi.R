#!/usr/bin/env Rscript
# CardioIsoMap — Step 3: batch-render per-gene coverage / sashimi plots (Gviz).
#
# Produces one PNG per gene: per-phenotype read coverage (NFH/DCM/IHD) with
# splice-junction (sashimi) arcs, over transcript models (known = red, novel PB = green).
# Based on gene_example/01_plot_myom1.R, but imports the GTF once and loops over genes.
#
# Usage:
#   Rscript code/03_render_sashimi.R [--genes=MYOM2,MYH7] [--limit=N] [--max-span=3000000] [--overwrite]
#
# Inputs (defaults point at ../gene_example):
#   BAMs: mapped_{NFH,DCM,IHD}.bam   GTF: gencode_v43_plus_novel_lr_may12.gtf.gz
#   Region table: data/processed/gene_regions.csv  (gene,chrom,plot_start,plot_end,...)
# Output:
#   data/processed/sashimi/<safe_gene>.png  + sashimi_manifest.csv

suppressMessages({ library(Gviz); library(rtracklayer) })
options(ucscChromosomeNames = FALSE)

# ---- paths ------------------------------------------------------------------
ROOT     <- normalizePath(file.path(dirname(sub("--file=", "",
              grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(ROOT) || ROOT == "") ROOT <- normalizePath(".")
GENE_EX  <- file.path(ROOT, "..", "gene_example")
PROC     <- file.path(ROOT, "data", "processed")
OUTDIR   <- file.path(PROC, "sashimi")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GTF_PATH <- file.path(GENE_EX, "gencode_v43_plus_novel_lr_may12.gtf.gz")
BAMS <- list(NFH = file.path(GENE_EX, "mapped_NFH.bam"),
             DCM = file.path(GENE_EX, "mapped_DCM.bam"),
             IHD = file.path(GENE_EX, "mapped_IHD.bam"))
COLOURS <- list(NFH = "#2196F3", DCM = "#F44336", IHD = "#4CAF50")

# ---- args -------------------------------------------------------------------
args <- commandArgs(TRUE)
getarg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
genes_sel <- getarg("--genes", "")
limit     <- as.integer(getarg("--limit", "0"))
max_span  <- as.numeric(getarg("--max-span", "3000000"))
shard     <- getarg("--shard", "")          # "i:n" -> render i-th of n contiguous chunks
overwrite <- any(args == "--overwrite")

safe_name <- function(g) gsub("[^A-Za-z0-9._-]", "_", g)

# ---- region table -----------------------------------------------------------
reg <- read.csv(file.path(PROC, "gene_regions.csv"), stringsAsFactors = FALSE)
if (nzchar(genes_sel)) reg <- reg[reg$gene %in% strsplit(genes_sel, ",")[[1]], ]
reg <- reg[reg$span <= max_span, ]                     # drop fusion-scale spans
reg <- reg[order(-reg$n_tx), ]
if (limit > 0) reg <- head(reg, limit)
if (nzchar(shard)) {                                    # split into n chunks, take i-th
  si <- as.integer(strsplit(shard, ":")[[1]])
  grp <- (seq_len(nrow(reg)) - 1) %% si[2] + 1          # interleave for balanced load
  reg <- reg[grp == si[1], ]
  message(sprintf("Shard %d/%d -> %d genes", si[1], si[2], nrow(reg)))
}
message(sprintf("Genes to render: %d (max_span=%s)", nrow(reg), format(max_span, big.mark=",")))

# ---- GTF once: exon GRanges with known/novel feature ------------------------
message("Importing GTF (once)...")
gtf_raw <- import(GTF_PATH)
exon_all <- gtf_raw[gtf_raw$type == "exon"]
exon_all$transcript <- exon_all$transcript_id
exon_all$gene       <- exon_all$gene_id
exon_all$symbol     <- ifelse(is.na(exon_all$gene_name),
                              exon_all$transcript_id, exon_all$gene_name)
exon_all$feature    <- ifelse(grepl("^PB\\.", exon_all$transcript_id), "lr", "known")
message(sprintf("  %d exon features loaded.", length(exon_all)))

axis_track <- GenomeAxisTrack(labelPos = "below", col = "black", showTitle = FALSE)

render_one <- function(row) {
  chr <- row$chrom; from <- row$plot_start; to <- row$plot_end
  out <- file.path(OUTDIR, paste0(safe_name(row$gene), ".png"))
  if (!overwrite && file.exists(out)) return("skip_exists")

  exon_gr <- subsetByOverlaps(exon_all, GRanges(chr, IRanges(from, to)))
  if (length(exon_gr) == 0) return("no_exons")
  n_tx <- length(unique(exon_gr$transcript_id))

  cov_tracks <- lapply(names(BAMS), function(s)
    AlignmentsTrack(range = BAMS[[s]], genome = "hg38", chromosome = chr,
                    isPaired = FALSE, name = s, type = c("coverage", "sashimi"),
                    fill = COLOURS[[s]], col.coverage = COLOURS[[s]],
                    col.sashimi = COLOURS[[s]], lwd.sashimi = 1.5, sashimiScore = 1,
                    col.axis = "black", ylab = s, showTitle = TRUE,
                    background.title = "white", col.title = COLOURS[[s]]))

  tx_track <- GeneRegionTrack(exon_gr, chromosome = chr, genome = "hg38", name = "",
                              showTitle = TRUE, background.title = "white",
                              col.title = "white", transcriptAnnotation = "transcript",
                              featureAnnotation = "feature", known = "#F44336",
                              lr = "#4CAF50", col = "transparent", cex.group = 0.6,
                              just.group = "left", clip = FALSE)

  all_tracks <- c(list(axis_track), cov_tracks, list(tx_track))
  h <- min(2600, 360 + 16 * n_tx)                       # scale height with #transcripts
  png(out, width = 1500, height = h, res = 130)
  tryCatch(
    plotTracks(all_tracks, chromosome = chr, from = from, to = to,
               sizes = c(0.3, 1, 1, 1, max(2, n_tx * 0.12)),
               main = paste0(row$gene, "  ", chr, ":",
                             format(from, big.mark = ","), "-", format(to, big.mark = ",")),
               cex.main = 1, title.width = 1),
    finally = dev.off())
  "ok"
}

# ---- loop -------------------------------------------------------------------
status <- character(nrow(reg)); t0 <- Sys.time()
for (i in seq_len(nrow(reg))) {
  status[i] <- tryCatch(render_one(reg[i, ]), error = function(e) paste0("error:", conditionMessage(e)))
  if (i %% 25 == 0 || i == nrow(reg))
    message(sprintf("  [%d/%d] %s -> %s (%.1f min)", i, nrow(reg), reg$gene[i],
                    status[i], as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

manifest <- data.frame(gene = reg$gene, png = paste0(safe_name(reg$gene), ".png"),
                       chrom = reg$chrom, plot_start = reg$plot_start,
                       plot_end = reg$plot_end, n_tx = reg$n_tx, status = status)
mf_path <- if (nzchar(shard))
  file.path(PROC, sprintf("sashimi_manifest_shard%s.csv", strsplit(shard, ":")[[1]][1]))
  else file.path(PROC, "sashimi_manifest.csv")
if (nzchar(shard)) {
  write.csv(manifest, mf_path, row.names = FALSE)        # per-shard, merged later
} else if (file.exists(mf_path) && !nzchar(genes_sel) && limit == 0) {
  # full run: overwrite manifest
  write.csv(manifest, mf_path, row.names = FALSE)
} else if (file.exists(mf_path)) {
  old <- read.csv(mf_path, stringsAsFactors = FALSE)
  old <- old[!(old$gene %in% manifest$gene), ]
  write.csv(rbind(old, manifest), mf_path, row.names = FALSE)
} else {
  write.csv(manifest, mf_path, row.names = FALSE)
}
message("Done. Status table:"); print(table(status))
