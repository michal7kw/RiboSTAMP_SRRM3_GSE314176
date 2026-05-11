#####################################################################
# 05h_gaba_depleted_isoform_table.R
#####################################################################
# COMPLEMENTARY DELIVERABLE — isoforms DEPLETED in hippocampal GABAergic
# interneurons (i.e. abundantly expressed in non-GABA cells but rare or
# absent in GABA). The mirror of 05_interneuron_isoform_table.R.
#
# Why this exists: 05_ finds GABA-OVEREXPRESSED isoforms (the 274 list).
# That answers "which isoforms ARE in interneurons?" but does NOT answer
# "which isoforms should be DEPLETED in IN-ribotag samples relative to
# bulk hippocampus?" — i.e. the negative-control / cell-type-contamination
# markers that an IN-ribotag QC needs.
#
# Filter (mirror of 05_, with logFC sign flipped):
#   other_pct >= 10   AND   logFC <= -1
# i.e. detected in >=10% of NON-GABA cells AND >=2x lower in GABA than
# non-GABA. Top hits are textbook glial markers (astrocyte / oligodendrocyte).
#
# Output:
#   results/tables/gaba_depleted_isoforms.tsv
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) NULL
)
if (is.null(SCRIPT_DIR) || SCRIPT_DIR == "" || SCRIPT_DIR == ".") {
  SCRIPT_DIR <- file.path(getwd(), "local")
}
BASE_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
CACHE_DIR <- file.path(BASE_DIR, "cache")
TABLES    <- file.path(BASE_DIR, "results", "tables")
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

# Filter thresholds — symmetric mirror of 05_'s MIN_GABA_DETECTION_PCT / MIN_LOGFC.
# Detection gate is on the OTHER (non-GABA) side because that is where the
# isoform should be abundant; logFC gate is flipped to negative.
MIN_OTHER_DETECTION_PCT <- 10
MIN_NEG_LOGFC           <- 1   # i.e. logFC <= -1

main <- function() {
  meta_path   <- file.path(CACHE_DIR, "cell_meta.rds")
  pooled_path <- file.path(CACHE_DIR, "pooled_normed_counts.rds")
  if (!file.exists(meta_path))   stop("Run local/01_parse_geo_metadata.R first.")
  if (!file.exists(pooled_path)) stop("Run local/05_interneuron_isoform_table.R first ",
                                      "(it builds cache/pooled_normed_counts.rds).")

  meta     <- readRDS(meta_path)
  pooled   <- readRDS(pooled_path)
  counts   <- pooled$counts
  features <- pooled$features

  common_bc <- intersect(colnames(counts), meta$barcode)
  message(sprintf("Cells in matrix and metadata intersection: %d / %d / %d",
                  length(common_bc), ncol(counts), nrow(meta)))
  counts <- counts[, common_bc, drop = FALSE]
  meta   <- meta[match(common_bc, meta$barcode)]

  gaba_mask  <- meta$is_gaba
  other_mask <- !meta$is_gaba
  if (sum(gaba_mask) == 0) {
    stop("No GABA cells found — check is_gaba derivation in 01_parse_geo_metadata.R.")
  }
  message(sprintf("GABA cells: %d  Non-GABA: %d", sum(gaba_mask), sum(other_mask)))

  gaba_counts  <- counts[, gaba_mask,  drop = FALSE]
  other_counts <- counts[, other_mask, drop = FALSE]

  GABA_pct        <- 100 * Matrix::rowSums(gaba_counts  > 0) / ncol(gaba_counts)
  GABA_mean       <- Matrix::rowMeans(gaba_counts)
  other_pct       <- 100 * Matrix::rowSums(other_counts > 0) / ncol(other_counts)
  other_mean      <- Matrix::rowMeans(other_counts)
  logFC           <- log2((GABA_mean + 1e-6) / (other_mean + 1e-6))
  n_GABA_detect   <- Matrix::rowSums(gaba_counts  > 0)
  n_other_detect  <- Matrix::rowSums(other_counts > 0)

  out <- data.table(
    isoform_id       = rownames(counts),
    gene             = features$gene_name[match(rownames(counts), features$feature_id)],
    transcript_id    = features$transcript_id[match(rownames(counts), features$feature_id)],
    biotype          = features$transcript_biotype[match(rownames(counts), features$feature_id)],
    GABA_pct         = GABA_pct,
    GABA_mean        = GABA_mean,
    other_pct        = other_pct,
    other_mean       = other_mean,
    logFC            = logFC,
    n_GABA_detected  = n_GABA_detect,
    n_other_detected = n_other_detect
  )

  keep <- out[other_pct >= MIN_OTHER_DETECTION_PCT & logFC <= -MIN_NEG_LOGFC]
  setorder(keep, logFC, -other_pct)   # most-depleted first

  out_path <- file.path(TABLES, "gaba_depleted_isoforms.tsv")
  fwrite(keep, out_path, sep = "\t")
  message(sprintf("Wrote %d isoforms to %s", nrow(keep), out_path))

  # Biotype breakdown — useful sanity check (glia markers are mostly protein_coding)
  message("\nBiotype breakdown:")
  print(keep[, .N, by = biotype][order(-N)])

  # Positive controls: canonical glial markers should appear here.
  ctrl_genes <- c("Slc1a3", "Sox9", "Hepacam", "Aqp4", "Aldh1l1",
                  "Mbp", "Mog", "Olig2", "Plp1",
                  "Gfap", "Slc1a2",
                  "Cx3cr1", "P2ry12", "Tmem119")
  hits <- ctrl_genes[ctrl_genes %in% unique(keep$gene)]
  miss <- setdiff(ctrl_genes, hits)
  message(sprintf("\nGlial positive-control markers present: %d / %d",
                  length(hits), length(ctrl_genes)))
  message(sprintf("  hits:   %s", paste(hits, collapse = ", ")))
  if (length(miss)) {
    message(sprintf("  missed: %s (probably not in catalog or below detection)",
                    paste(miss, collapse = ", ")))
  }

  message("\nTop 10 most-depleted isoforms (by logFC):")
  print(keep[1:10, .(isoform_id, gene, GABA_pct, other_pct, GABA_mean, other_mean, logFC)])

  invisible(keep)
}

if (sys.nframe() == 0) main()
