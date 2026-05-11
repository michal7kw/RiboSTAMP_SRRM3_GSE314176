#####################################################################
# 06_editsc_translation_column.R
#####################################################################
# Bonus — add EditsC translation signal to the interneuron table.
#
# Loads the 3 mice's EditsC h5ads (note: each has 12,212 transcripts, FEWER
# than the 24,309 in normed_counts — different feature universe, so we
# left-join by isoform_id and surface NAs where EditsC is undefined).
#
# Adds columns to results/tables/interneuron_isoforms.tsv:
#   editsC_GABA_mean, editsC_GABA_pct, editsC_is_translated (logical)
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

MICE <- c("GSM9380799", "GSM9380800", "GSM9380801")
MIN_EDITSC_DETECTION_PCT <- 5

load_mouse_editsc <- function(mouse_gsm) {
  prefix <- paste0(mouse_gsm, "_longread_EditsC_transcript_adata_")
  mtx_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__X\\.mtx$"),
                         full.names = TRUE)[1]
  var_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__var\\.tsv$"),
                         full.names = TRUE)[1]
  obs_file <- list.files(META_DIR, pattern = paste0("^", prefix, ".*__obs\\.tsv$"),
                         full.names = TRUE)[1]
  if (any(is.na(c(mtx_file, var_file, obs_file)))) {
    stop("Missing exploded EditsC files for ", mouse_gsm,
         " in ", META_DIR, " — re-run pc_pipeline/01_fetch_processed.sh.")
  }
  m <- Matrix::readMM(mtx_file)
  features <- fread(var_file)
  barcodes <- fread(obs_file)
  rownames(m) <- features$feature_id
  colnames(m) <- barcodes$barcode
  list(editsc = as(m, "CsparseMatrix"), features = features)
}

load_pooled_editsc <- function() {
  cache <- file.path(CACHE_DIR, "pooled_editsc.rds")
  if (file.exists(cache)) return(readRDS(cache))

  parts <- lapply(MICE, load_mouse_editsc)
  features_per_mouse <- lapply(parts, function(p) p$features$feature_id)
  common_features    <- Reduce(intersect, features_per_mouse)
  message(sprintf("Pooled EditsC feature set: %d transcripts (intersection across %d mice)",
                  length(common_features), length(MICE)))

  parts_aligned <- lapply(parts, function(p) p$editsc[common_features, , drop = FALSE])
  editsc <- do.call(cbind, parts_aligned)

  saveRDS(editsc, cache)
  editsc
}

main <- function() {
  out_path <- file.path(TABLES, "interneuron_isoforms.tsv")
  if (!file.exists(out_path)) {
    stop("Run local/05_interneuron_isoform_table.R first — missing ", out_path)
  }
  existing <- fread(out_path)

  meta_path <- file.path(CACHE_DIR, "cell_meta.rds")
  if (!file.exists(meta_path)) stop("Run local/01_parse_geo_metadata.R first.")
  meta <- readRDS(meta_path)

  editsc <- load_pooled_editsc()

  common_bc <- intersect(colnames(editsc), meta$barcode)
  editsc <- editsc[, common_bc, drop = FALSE]
  meta   <- meta[match(common_bc, meta$barcode)]
  gaba_mask <- meta$is_gaba

  if (sum(gaba_mask) == 0) {
    stop("No GABA cells found in EditsC matrix — check is_gaba derivation.")
  }

  editsc_gaba <- editsc[, gaba_mask, drop = FALSE]

  editsc_stats <- data.table(
    isoform_id       = rownames(editsc),
    editsC_GABA_mean = Matrix::rowMeans(editsc_gaba),
    editsC_GABA_pct  = 100 * Matrix::rowSums(editsc_gaba > 0) / ncol(editsc_gaba)
  )
  editsc_stats[, editsC_is_translated := editsC_GABA_pct >= MIN_EDITSC_DETECTION_PCT]

  # Left join — many isoforms in the counts matrix have no EditsC entry
  # (different feature set: 24309 vs 12212). NA = EditsC undefined.
  merged <- merge(existing, editsc_stats, by = "isoform_id",
                  all.x = TRUE, sort = FALSE)
  fwrite(merged, out_path, sep = "\t")

  n_with_editsc      <- sum(!is.na(merged$editsC_GABA_mean))
  n_translated       <- sum(merged$editsC_is_translated, na.rm = TRUE)
  n_tx_and_tl        <- sum(merged$GABA_pct >= 10 & merged$editsC_is_translated, na.rm = TRUE)
  message(sprintf("Added EditsC columns. %d/%d isoforms have any EditsC measurement.",
                  n_with_editsc, nrow(merged)))
  message(sprintf("                    %d / %d isoforms show translation evidence (>=%d%% of GABA cells).",
                  n_translated, nrow(merged), MIN_EDITSC_DETECTION_PCT))
  message(sprintf("                    %d isoforms transcribed (>=10%%) AND translated in GABA.",
                  n_tx_and_tl))

  invisible(merged)
}

if (sys.nframe() == 0) main()
