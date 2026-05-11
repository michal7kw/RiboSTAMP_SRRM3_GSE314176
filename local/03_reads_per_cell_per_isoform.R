#####################################################################
# 03_reads_per_cell_per_isoform.R
#####################################################################
# Build a sparse per-cell × per-isoform read count matrix from the
# IsoQuant per-read assignment TSV, keyed by author_barcode (after
# RC correction from 02_barcode_rc_map.R).
#
# Output cache: cache/iso_counts.rds (dgCMatrix).
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
ISOQ_DIR  <- file.path(BASE_DIR, "data", "isoquant_targeted")

build_iso_counts <- function(force = FALSE) {
  cache_file <- file.path(CACHE_DIR, "iso_counts.rds")
  if (!force && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  bc_map_cache <- file.path(CACHE_DIR, "barcode_map.rds")
  if (!file.exists(bc_map_cache)) {
    stop("Run 02_barcode_rc_map.R first — missing ", bc_map_cache)
  }
  bc_map <- readRDS(bc_map_cache)
  setkey(bc_map, long_read_cb)

  assign_files <- list.files(ISOQ_DIR, pattern = "read_assignments\\.tsv(\\.gz)?$",
                             recursive = TRUE, full.names = TRUE)
  all <- rbindlist(lapply(assign_files, function(f) {
    dt <- fread(f)
    cb_col  <- intersect(c("group_id", "read_group", "CB"),         names(dt))[1]
    iso_col <- intersect(c("isoform_id", "transcript_id", "isoform"), names(dt))[1]
    type_col<- intersect(c("assignment_type", "transcript_type"),   names(dt))[1]
    data.table(
      long_read_cb    = as.character(dt[[cb_col]]),
      isoform_id      = as.character(dt[[iso_col]]),
      assignment_type = if (!is.na(type_col)) as.character(dt[[type_col]]) else NA_character_
    )
  }))
  all <- all[long_read_cb != "." & nzchar(long_read_cb)]
  all <- all[!is.na(isoform_id) & nzchar(isoform_id)]

  # Keep only confident (unique) read-to-isoform assignments — IsoQuant's
  # "unique" / "unique_minor_difference" categories. Ambiguous reads
  # cannot be allocated cleanly and would inflate counts.
  keep_types <- c("unique", "unique_minor_difference",
                  "reference_match", "novel_in_catalog")
  if (any(!is.na(all$assignment_type))) {
    all <- all[assignment_type %in% keep_types |
               is.na(assignment_type)]
  }

  # RC-correct long-read CB -> author barcode
  all[, author_barcode := bc_map[all$long_read_cb, on = "long_read_cb"]$author_barcode]
  all <- all[!is.na(author_barcode)]

  # Sparse matrix: rows = isoforms, cols = author_barcodes
  isoforms <- sort(unique(all$isoform_id))
  barcodes <- sort(unique(all$author_barcode))
  counts <- Matrix::sparseMatrix(
    i = match(all$isoform_id,     isoforms),
    j = match(all$author_barcode, barcodes),
    x = 1,
    dims = c(length(isoforms), length(barcodes)),
    dimnames = list(isoforms, barcodes)
  )
  # sparseMatrix collapses duplicate (i,j) by summing, so the above is
  # equivalent to "reads per isoform per cell" already.

  message(sprintf("Built %d isoforms x %d cells matrix.", nrow(counts), ncol(counts)))

  saveRDS(counts, cache_file)
  counts
}

if (sys.nframe() == 0) {
  mat <- build_iso_counts(force = FALSE)
  cat(sprintf("iso_counts.rds: %d x %d, %d non-zero entries.\n",
              nrow(mat), ncol(mat), length(mat@x)))
}
