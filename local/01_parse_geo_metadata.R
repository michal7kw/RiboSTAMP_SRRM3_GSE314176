#####################################################################
# 01_parse_geo_metadata.R
#####################################################################
# Pool the cell-type metadata from the three GSE314176 long-read mice
# (GSM9380799 / GSM9380800 / GSM9380801) into a single tibble cached
# at cache/cell_meta.rds. Used by all downstream local R scripts.
#
# The h5ads expose two relevant cluster columns:
#   "Cell Assignment"          fine-grained:  e.g. "Neuron - GABA_1",
#                                                  "Neuron - GABA_2",
#                                                  "Neuron - DG - Immature",
#                                                  "Glia - OPC", "Glia - Astro_4"
#   "Cell Assignments Grouped" coarse:        e.g. "Neuron - GABA",
#                                                  "Neuron - DG - Immature",
#                                                  "Glia - OPC"
# We carry both, and derive the analysis flags below.
#
# Output columns:
#   barcode             cell barcode (from the h5ad index, e.g. "ATCG..._sample3")
#   mouse               GSM accession (GSM9380799 / 800 / 801)
#   sample              author "sample" string from obs (sample1/2/3)
#   cell_assignment     fine cluster (kept verbatim)
#   cell_type           coarse cluster (kept verbatim)
#   leiden              author's leiden cluster integer
#   broad_class         {"progenitor","interneuron","neuron","glia","bbb","immune","other"}
#   is_progenitor       TRUE for OPC + iDG cells (the Track B target population)
#   is_gaba             TRUE for any GABAergic interneuron
#   gaba_subtype        NA / "GABA_1" / "GABA_2" (the actual subdivisions in this dataset)
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# Resolve script dir (rstudioapi / cwd fallback — matches Neuroblastoma scripts)
SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) NULL
)
if (is.null(SCRIPT_DIR) || SCRIPT_DIR == "" || SCRIPT_DIR == ".") {
  SCRIPT_DIR <- file.path(getwd(), "local")
}
BASE_DIR  <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
DATA_META <- file.path(BASE_DIR, "data", "metadata")
CACHE_DIR <- file.path(BASE_DIR, "cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CACHE_FILE <- file.path(CACHE_DIR, "cell_meta.rds")

#' Load the per-mouse obs.tsv (from the normed_counts h5ad, since both
#' counts and EditsC h5ads carry the same author cell metadata) and
#' return a tagged data.table.
load_mouse_obs <- function(mouse_gsm) {
  obs_path <- list.files(
    DATA_META,
    pattern = paste0(mouse_gsm, "_longread_normed_counts_.*__obs\\.tsv$"),
    full.names = TRUE
  )
  if (length(obs_path) == 0) {
    stop("No obs.tsv for ", mouse_gsm, " under ", DATA_META,
         ". Run hpc/01_fetch_processed.sh (or pc_pipeline equivalent) and re-explode.")
  }
  obs <- fread(obs_path[1])
  obs[, mouse := mouse_gsm]
  obs
}

#' Strip the 10x barcode suffix to get the bare 16-mer.
#' Shortread uses "<barcode>-1"; longread uses "<barcode>_sampleN".
strip_barcode_suffix <- function(x) sub("[_\\-].*$", "", x)

#' Load the short-read AVITI obs for the same biological mouse and return
#' a barcode-keyed table with the bare 16-mer + sr cluster annotations.
#'
#' Background: shortread and longread were generated from the same biological
#' mouse but their barcode suffix conventions differ (-1 vs _sampleN), and
#' the libraries captured DIFFERENT subsets of cells (84-94% overlap on the
#' bare 16-mer, but the rare OPC + iDG-Immature populations were NOT
#' captured by the PacBio long-read library at all). So this join is
#' useful for cross-validating cluster labels and for fine-grained
#' "Cell Assignment" subdivisions, but it CANNOT rescue OPC + iDG into
#' the long-read data — those cells are simply absent.
#'
#' Mouse pairing: 796↔799 (mouse1), 797↔800 (mouse2), 798↔801 (mouse3).
load_shortread_obs <- function(mouse_pair) {
  shortread_gsm <- c(GSM9380799 = "GSM9380796",
                     GSM9380800 = "GSM9380797",
                     GSM9380801 = "GSM9380798")[[mouse_pair]]
  obs_path <- list.files(
    DATA_META,
    pattern = paste0(shortread_gsm, "_shortread_.*__obs\\.tsv$"),
    full.names = TRUE
  )
  if (length(obs_path) == 0) return(NULL)
  sr <- fread(obs_path[1])
  sr[, .(bare_barcode      = strip_barcode_suffix(barcode),
         sr_cell_assignment = `Cell Assignment`,
         sr_cell_type       = `Cell Assignments Grouped`)]
}

load_cell_meta <- function(force = FALSE) {
  if (!force && file.exists(CACHE_FILE)) {
    message("Loading cached cell metadata from ", CACHE_FILE)
    return(readRDS(CACHE_FILE))
  }

  parts <- lapply(c("GSM9380799", "GSM9380800", "GSM9380801"), load_mouse_obs)
  all <- rbindlist(parts, use.names = TRUE, fill = TRUE)

  # Normalize the two cluster columns. Keep originals.
  setnames(all, "Cell Assignment",          "cell_assignment", skip_absent = TRUE)
  setnames(all, "Cell Assignments Grouped", "cell_type",       skip_absent = TRUE)

  # Strip whitespace
  all[, cell_assignment := str_squish(cell_assignment)]
  all[, cell_type       := str_squish(cell_type)]

  # ------------------------------------------------------------------
  # Cross-reference with short-read AVITI obs (same biological mice).
  # 84-94% of long-read cells have a matching short-read cell on the
  # bare 16-mer barcode. The cluster labels agree where they overlap.
  # NB: this DOES NOT rescue OPC + iDG cells — those rare populations
  # were captured by the short-read library but NOT by the PacBio
  # long-read library. The join is useful for cross-validation, not
  # cell-type rescue. See docs/LEARNING_06_novelty_layered_defense.md
  # §"Progenitor cells in this dataset" for the implication on Track B.
  # ------------------------------------------------------------------
  all[, bare_barcode := strip_barcode_suffix(barcode)]

  # Build a unified short-read obs table tagged by long-read GSM, then
  # do a single left-merge by (mouse, bare_barcode). Cleaner + avoids the
  # length-mismatch trap of per-mouse data.table updates.
  sr_pairs <- list(GSM9380799 = "GSM9380796",
                   GSM9380800 = "GSM9380797",
                   GSM9380801 = "GSM9380798")
  sr_parts <- list()
  for (lr_gsm in names(sr_pairs)) {
    sr <- load_shortread_obs(lr_gsm)
    if (is.null(sr)) next
    sr[, mouse := lr_gsm]
    sr_parts[[lr_gsm]] <- sr
  }
  if (length(sr_parts)) {
    sr_all <- rbindlist(sr_parts, use.names = TRUE, fill = TRUE)
    all <- merge(all, sr_all,
                 by = c("mouse", "bare_barcode"),
                 all.x = TRUE, sort = FALSE)
  } else {
    all[, sr_cell_assignment := NA_character_]
    all[, sr_cell_type       := NA_character_]
  }

  # Derive analysis flags. is_progenitor is FALSE for everything in this
  # dataset (the long-read library has zero OPC + iDG cells); we keep
  # the column so downstream scripts that test it don't break, but the
  # progenitor enrichment analysis must reframe to "in which clusters
  # is the novel isoform expressed?" rather than OPC+iDG vs other.
  all[, is_progenitor := cell_type %in% c("Glia - OPC", "Neuron - DG - Immature")]

  all[, is_gaba := cell_type == "Neuron - GABA"]

  # Subtype: for GABA we have GABA_1 and GABA_2 in cell_assignment
  all[, gaba_subtype := dplyr::case_when(
    grepl("Neuron - GABA_1", cell_assignment, fixed = TRUE) ~ "GABA_1",
    grepl("Neuron - GABA_2", cell_assignment, fixed = TRUE) ~ "GABA_2",
    is_gaba                                                  ~ "other_GABA",
    TRUE                                                     ~ NA_character_
  )]

  all[, broad_class := dplyr::case_when(
    is_progenitor                                                  ~ "progenitor",
    is_gaba                                                        ~ "interneuron",
    grepl("^Neuron",  cell_type)                                   ~ "neuron",
    grepl("^Glia",    cell_type)                                   ~ "glia",
    grepl("^BBB",     cell_type)                                   ~ "bbb",
    grepl("^Immune",  cell_type)                                   ~ "immune",
    TRUE                                                           ~ "other"
  )]

  message(sprintf("Loaded %d cells across %d unique cluster labels from 3 mice.",
                  nrow(all), length(unique(all$cell_type))))
  message("Broad-class breakdown:")
  print(all[, .N, by = broad_class][order(-N)])
  message("Progenitor cells (OPC + iDG) per mouse:")
  print(all[is_progenitor == TRUE, .N, by = .(mouse, cell_type)])
  message("GABA cells per mouse and subtype:")
  print(all[is_gaba == TRUE, .N, by = .(mouse, gaba_subtype)])

  saveRDS(all, CACHE_FILE)
  all
}

if (sys.nframe() == 0) {
  meta <- load_cell_meta(force = TRUE)
  cat(sprintf("\nCached to %s (%d rows, %d cols).\n",
              CACHE_FILE, nrow(meta), ncol(meta)))
}
