#####################################################################
# 05_interneuron_isoform_table.R
#####################################################################
# DELIVERABLE 2 — isoforms expressed in hippocampal GABAergic interneurons,
# with per-subtype breakout (GABA_1 / GABA_2 — the actual subdivisions in
# this dataset; the authors did NOT subtype by Pvalb/Sst/Vip/Lamp5).
#
# Inputs (per mouse, from the 3 long-read normed-counts h5ads):
#   data/metadata/GSM93807{99,800,801}_longread_normed_counts_*__X.mtx
#   data/metadata/...__obs.tsv  (loaded via 01_parse_geo_metadata.R)
#   data/metadata/...__var.tsv
#
# We pool cells across all 3 mice (no batch correction — these are
# biological replicates from the same protocol; pooling improves
# detection statistics for rare subtypes).
#
# Output:
#   results/tables/interneuron_isoforms.tsv
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
META_DIR  <- file.path(BASE_DIR, "data", "metadata")
TABLES    <- file.path(BASE_DIR, "results", "tables")
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)

MICE <- c("GSM9380799", "GSM9380800", "GSM9380801")

# Filter thresholds (documented in LEARNING_04_interneuron_design.md)
MIN_GABA_DETECTION_PCT <- 10
MIN_LOGFC              <- 1

#' Load one mouse's transcript-level matrix (normed counts) + var table.
#' Returns list(counts, features) with counts as features × cells dgCMatrix.
load_mouse_counts <- function(mouse_gsm) {
  prefix <- paste0(mouse_gsm, "_longread_normed_counts_transcript_adata_")
  mtx_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__X\\.mtx$"),
                         full.names = TRUE)[1]
  var_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__var\\.tsv$"),
                         full.names = TRUE)[1]
  obs_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__obs\\.tsv$"),
                         full.names = TRUE)[1]
  if (any(is.na(c(mtx_file, var_file, obs_file)))) {
    stop("Missing exploded normed_counts files for ", mouse_gsm,
         " in ", META_DIR, " — re-run pc_pipeline/01_fetch_processed.sh.")
  }
  m <- Matrix::readMM(mtx_file)
  features <- fread(var_file)
  barcodes <- fread(obs_file)
  rownames(m) <- features$feature_id
  colnames(m) <- barcodes$barcode
  list(counts = as(m, "CsparseMatrix"), features = features)
}

#' Pool counts across mice.
#' - Intersect features (the var table is identical across mice in practice
#'   but we intersect defensively).
#' - cbind cells (barcodes have a "_sampleN" suffix in obs so they're already unique).
load_pooled_counts <- function() {
  cache <- file.path(CACHE_DIR, "pooled_normed_counts.rds")
  if (file.exists(cache)) return(readRDS(cache))

  parts <- lapply(MICE, load_mouse_counts)
  features_per_mouse <- lapply(parts, function(p) p$features$feature_id)
  common_features    <- Reduce(intersect, features_per_mouse)
  message(sprintf("Pooled feature set: %d transcripts (intersection across %d mice)",
                  length(common_features), length(MICE)))

  parts_aligned <- lapply(parts, function(p) p$counts[common_features, , drop = FALSE])
  counts <- do.call(cbind, parts_aligned)

  features <- parts[[1]]$features[parts[[1]]$features$feature_id %in% common_features]
  setkey(features, feature_id)
  features <- features[J(common_features), on = "feature_id"]

  out <- list(counts = counts, features = features)
  saveRDS(out, cache)
  out
}

main <- function() {
  meta_path <- file.path(CACHE_DIR, "cell_meta.rds")
  if (!file.exists(meta_path)) {
    stop("Run local/01_parse_geo_metadata.R first.")
  }
  meta <- readRDS(meta_path)

  pooled <- load_pooled_counts()
  counts <- pooled$counts
  features <- pooled$features

  # Align meta to count matrix columns
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

  GABA_pct       <- 100 * Matrix::rowSums(gaba_counts  > 0) / ncol(gaba_counts)
  GABA_mean      <- Matrix::rowMeans(gaba_counts)
  other_pct      <- 100 * Matrix::rowSums(other_counts > 0) / ncol(other_counts)
  other_mean     <- Matrix::rowMeans(other_counts)
  logFC          <- log2((GABA_mean + 1e-6) / (other_mean + 1e-6))
  n_GABA_detect  <- Matrix::rowSums(gaba_counts > 0)

  subtype_stats <- function(subtype_label) {
    mask <- meta$is_gaba & !is.na(meta$gaba_subtype) & meta$gaba_subtype == subtype_label
    if (sum(mask) == 0) {
      return(list(pct  = rep(NA_real_, nrow(counts)),
                  mean = rep(NA_real_, nrow(counts))))
    }
    sub <- counts[, mask, drop = FALSE]
    list(
      pct  = 100 * Matrix::rowSums(sub > 0) / ncol(sub),
      mean = Matrix::rowMeans(sub)
    )
  }
  g1 <- subtype_stats("GABA_1")
  g2 <- subtype_stats("GABA_2")

  out <- data.table(
    isoform_id      = rownames(counts),
    gene            = features$gene_name[match(rownames(counts), features$feature_id)],
    transcript_id   = features$transcript_id[match(rownames(counts), features$feature_id)],
    GABA_pct        = GABA_pct,
    GABA_mean       = GABA_mean,
    other_pct       = other_pct,
    other_mean      = other_mean,
    logFC           = logFC,
    n_GABA_detected = n_GABA_detect,
    GABA_1_pct = g1$pct, GABA_1_mean = g1$mean,
    GABA_2_pct = g2$pct, GABA_2_mean = g2$mean
  )

  keep <- out[GABA_pct >= MIN_GABA_DETECTION_PCT & logFC >= MIN_LOGFC]
  setorder(keep, -logFC, -GABA_pct)

  fwrite(keep, file.path(TABLES, "interneuron_isoforms.tsv"), sep = "\t")
  message(sprintf("Wrote %d isoforms to %s",
                  nrow(keep), file.path(TABLES, "interneuron_isoforms.tsv")))

  # Positive-control gates
  ctrl_genes <- c("Gad1", "Gad2", "Slc6a1", "Slc32a1")  # canonical GABAergic markers
  hits <- ctrl_genes[ctrl_genes %in% unique(keep$gene)]
  miss <- setdiff(ctrl_genes, hits)
  if (length(miss)) {
    warning(sprintf("Positive-control GABA marker genes MISSING: %s. ",
                    paste(miss, collapse = ", ")),
            "Consider relaxing thresholds or verifying cluster labels.")
  } else {
    message("Positive-control check PASSED — Gad1, Gad2, Slc6a1, Slc32a1 all present.")
  }

  invisible(keep)
}

if (sys.nframe() == 0) main()
