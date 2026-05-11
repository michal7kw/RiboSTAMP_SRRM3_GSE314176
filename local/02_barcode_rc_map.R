#####################################################################
# 02_barcode_rc_map.R
#####################################################################
# BLOCKING sanity check: confirm the RC orientation between long-read
# BAM CB tags and author obs.csv barcodes, then build a conversion map.
#
# The paper notes that MAS-ISO-seq reads the TSO-adjacent barcode on the
# opposite strand from 10x short-read; in practice, one of the following
# must hold with >80% intersection:
#   (A) author_barcode == long_read_CB
#   (B) author_barcode == reverse_complement(long_read_CB)
# AND the OTHER orientation must have <5% intersection (asymmetry check).
#
# If neither is true, we are mapping the wrong barcode set. EVERYTHING
# downstream (progenitor enrichment, interneuron list) depends on this —
# silent mis-orientation produces random cluster assignments and
# statistically looks like "no signal". Fail fast, fail loud.
#
# See docs/LEARNING_02_barcode_rc_gotcha.md for the full rationale.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
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
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds — documented in LEARNING_02
MIN_MATCH_FRACTION <- 0.80   # correct orientation must exceed this
MAX_WRONG_FRACTION <- 0.05   # wrong orientation must be below this

reverse_complement <- function(x) {
  # Vectorized RC using chartr + rev on each split string.
  complement <- chartr("ACGTNacgtn", "TGCANTGCAN", x)
  vapply(strsplit(complement, "", fixed = TRUE),
         function(cs) paste(rev(cs), collapse = ""),
         character(1))
}

#' Build + assert a barcode mapping between long-read CB tags and
#' author obs.csv barcodes.
#'
#' @return data.table with columns (long_read_cb, author_barcode, rc_applied)
#'   plus attr "orientation" = "identity" | "reverse_complement"
build_barcode_map <- function(force = FALSE) {
  cache_file <- file.path(CACHE_DIR, "barcode_map.rds")
  if (!force && file.exists(cache_file)) {
    message("Loading cached barcode map from ", cache_file)
    return(readRDS(cache_file))
  }

  # 1. Load author barcodes (from obs.csv via 01_parse_geo_metadata.R cache)
  meta_cache <- file.path(CACHE_DIR, "cell_meta.rds")
  if (!file.exists(meta_cache)) {
    stop("Run 01_parse_geo_metadata.R first (missing ", meta_cache, ").")
  }
  author_barcodes <- unique(readRDS(meta_cache)$barcode)

  # 2. Load long-read CB tags from IsoQuant read-assignment TSV
  #    IsoQuant writes one file per BAM; columns include read_id, isoform_id,
  #    transcript_type, and — when --read_group tag:CB — a group column.
  assign_files <- list.files(ISOQ_DIR, pattern = "read_assignments\\.tsv(\\.gz)?$",
                             recursive = TRUE, full.names = TRUE)
  if (length(assign_files) == 0) {
    stop("No IsoQuant read-assignment TSVs under ", ISOQ_DIR,
         " — run hpc/05_isoquant_targeted.sh and sync outputs to data/isoquant_targeted/")
  }

  long_read_cbs <- unique(unlist(lapply(assign_files, function(f) {
    dt <- fread(f)
    # IsoQuant's group column name can vary; pick the first non-standard column
    # matching our --read_group tag:CB output. Defensive search.
    cand <- intersect(c("group_id", "read_group", "CB"), names(dt))
    if (length(cand) == 0) {
      stop("Could not locate CB tag column in ", f,
           " — check IsoQuant --read_group output format.")
    }
    as.character(dt[[cand[1]]])
  })))
  long_read_cbs <- long_read_cbs[long_read_cbs != "." & nzchar(long_read_cbs)]

  message(sprintf("Author barcodes: %d unique. Long-read CB tags: %d unique.",
                  length(author_barcodes), length(long_read_cbs)))

  # 3. Compute intersections in both orientations
  long_read_cbs_rc <- reverse_complement(long_read_cbs)
  inter_identity  <- length(intersect(author_barcodes, long_read_cbs))
  inter_rc        <- length(intersect(author_barcodes, long_read_cbs_rc))

  denom <- min(length(author_barcodes), length(long_read_cbs))
  frac_identity <- inter_identity / denom
  frac_rc       <- inter_rc       / denom

  message(sprintf("Intersection fraction (identity):           %.3f", frac_identity))
  message(sprintf("Intersection fraction (reverse_complement): %.3f", frac_rc))

  # 4. Assert asymmetry
  orientation <- NA_character_
  if (frac_identity >= MIN_MATCH_FRACTION && frac_rc <= MAX_WRONG_FRACTION) {
    orientation <- "identity"
  } else if (frac_rc >= MIN_MATCH_FRACTION && frac_identity <= MAX_WRONG_FRACTION) {
    orientation <- "reverse_complement"
  } else {
    stop(sprintf(paste0(
      "BARCODE ORIENTATION CHECK FAILED.\n",
      "  identity:   %.1f%%  (need >= %.0f%% on the correct orientation)\n",
      "  RC:         %.1f%%  (need <= %.0f%% on the incorrect orientation)\n",
      "Neither orientation shows a clean match. Do not proceed — downstream\n",
      "progenitor/interneuron statistics will be random if orientation is wrong.\n",
      "Re-verify: (a) obs.csv rows correspond to the same 10x chemistry as the\n",
      "long-read library; (b) IsoQuant --read_group tag:CB actually wrote the\n",
      "cell barcode (not the UMI); (c) no accidental trimming of a leading/\n",
      "trailing base during BAM processing."),
      100 * frac_identity, 100 * MIN_MATCH_FRACTION,
      100 * frac_rc,       100 * MAX_WRONG_FRACTION))
  }

  message(sprintf("Orientation accepted: %s", orientation))

  # 5. Build the mapping
  map_dt <- data.table(
    long_read_cb   = long_read_cbs,
    author_barcode = if (orientation == "identity") long_read_cbs else long_read_cbs_rc,
    rc_applied     = (orientation == "reverse_complement")
  )
  attr(map_dt, "orientation")  <- orientation
  attr(map_dt, "frac_correct") <- if (orientation == "identity") frac_identity else frac_rc

  saveRDS(map_dt, cache_file)
  map_dt
}

if (sys.nframe() == 0) {
  bc_map <- build_barcode_map(force = FALSE)
  cat(sprintf("Barcode map cached: %d entries, orientation = %s\n",
              nrow(bc_map), attr(bc_map, "orientation")))
}
