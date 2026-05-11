# =============================================================================
# 05e_interneuron_functional_annotation.R — add gene biotype + GO + SynGO
# annotation columns to the interneuron table (robust version)
#
# Annotation sources, in priority order:
#   1. biotype: from the local var.tsv file (already has transcript_biotype)
#   2. GO BP: try biomaRt; gracefully fall back to org.Mm.eg.db if available;
#      otherwise skip
#   3. SynGO: hardcoded quick-list of well-known synaptic genes
#
# This makes the table filterable by functional context with one click.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

SCRIPT_DIR <- tryCatch(here::here(), error = function(e) getwd())
PROJ_ROOT <- if (basename(SCRIPT_DIR) == "local") dirname(SCRIPT_DIR) else SCRIPT_DIR
META_DIR <- file.path(PROJ_ROOT, "data", "metadata")
RESULTS_DIR <- file.path(PROJ_ROOT, "results", "tables")
CACHE_DIR <- file.path(PROJ_ROOT, "cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load source table
# -----------------------------------------------------------------------------
candidates <- c("interneuron_isoforms_de.tsv",
                "interneuron_isoforms_by_subtype.tsv",
                "interneuron_isoforms.tsv")
INO_PATH <- NULL
for (cand in candidates) {
  p <- file.path(RESULTS_DIR, cand)
  if (file.exists(p)) { INO_PATH <- p; break }
}
if (is.null(INO_PATH)) stop("No interneuron isoform table found")
message(sprintf("Source table: %s", INO_PATH))

ino <- fread(INO_PATH)
genes <- unique(ino$gene)
message(sprintf("Unique genes: %d, isoforms: %d", length(genes), nrow(ino)))

# -----------------------------------------------------------------------------
# Layer 1: biotype from local var.tsv (always works, no network needed)
# -----------------------------------------------------------------------------
message("\nLayer 1: biotype (from local var.tsv)")
var_path <- list.files(META_DIR,
                       pattern = "^GSM9380799_longread_normed_counts.*_var\\.tsv$",
                       full.names = TRUE)[1]
if (file.exists(var_path)) {
  var_dt <- fread(var_path)
  message(sprintf("  loaded var.tsv: %d isoforms, columns: %s",
                  nrow(var_dt), paste(names(var_dt), collapse = ", ")))
  if ("transcript_biotype" %in% names(var_dt) && "feature_id" %in% names(var_dt)) {
    biotype_lookup <- var_dt[, .(isoform_id = feature_id, biotype = transcript_biotype)]
    biotype_lookup <- unique(biotype_lookup, by = "isoform_id")
    ino <- merge(ino, biotype_lookup, by = "isoform_id", all.x = TRUE, sort = FALSE)
    n_biotype <- sum(!is.na(ino$biotype))
    message(sprintf("  annotated %d / %d isoforms with biotype", n_biotype, nrow(ino)))
  } else {
    message("  WARN: var.tsv missing transcript_biotype/feature_id columns")
    ino[, biotype := NA_character_]
  }
} else {
  message("  WARN: var.tsv not found; skipping biotype layer")
  ino[, biotype := NA_character_]
}

# -----------------------------------------------------------------------------
# Layer 2: GO BP — try biomaRt, fall back to org.Mm.eg.db, else skip
# -----------------------------------------------------------------------------
message("\nLayer 2: GO Biological Process")
go_cache <- file.path(CACHE_DIR, "biomart_go_bp.rds")

if (file.exists(go_cache)) {
  message("  using cached GO BP terms")
  go_dt <- readRDS(go_cache)
} else {
  go_dt <- NULL

  # Try biomaRt first
  message("  attempting biomaRt query...")
  ok <- tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE)) stop("biomaRt not installed")
    suppressPackageStartupMessages(library(biomaRt))
    mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
    bm <- getBM(
      attributes = c("external_gene_name", "name_1006", "namespace_1003"),
      filters = "external_gene_name",
      values = genes,
      mart = mart
    )
    bm <- as.data.table(bm)
    bm <- bm[namespace_1003 == "biological_process" & name_1006 != ""]
    go_dt <- bm[, .(go_bp = paste(head(unique(name_1006), 3), collapse = "; ")),
                by = .(gene = external_gene_name)]
    saveRDS(go_dt, go_cache)
    message("  biomaRt succeeded — cached")
    TRUE
  }, error = function(e) {
    message(sprintf("  biomaRt failed: %s", conditionMessage(e)))
    FALSE
  })

  # Fall back to org.Mm.eg.db if biomaRt failed
  if (!ok) {
    message("  trying org.Mm.eg.db fallback...")
    ok2 <- tryCatch({
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
          !requireNamespace("AnnotationDbi", quietly = TRUE) ||
          !requireNamespace("GO.db", quietly = TRUE)) {
        stop("org.Mm.eg.db / AnnotationDbi / GO.db not installed")
      }
      library(org.Mm.eg.db)
      library(AnnotationDbi)
      library(GO.db)
      # Map gene symbols → ENTREZ → GO terms
      sym2entrez <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL",
                           column = "ENTREZID", multiVals = "first")
      sym2entrez <- sym2entrez[!is.na(sym2entrez)]
      if (length(sym2entrez) == 0) stop("no entrez mapping found")
      entrez2go <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = unname(sym2entrez),
                                          keytype = "ENTREZID",
                                          columns = c("GOALL", "ONTOLOGYALL"))
      entrez2go <- as.data.table(entrez2go)
      entrez2go <- entrez2go[ONTOLOGYALL == "BP" & !is.na(GOALL)]
      # Translate GO IDs → term names
      go_terms <- AnnotationDbi::select(GO.db, keys = unique(entrez2go$GOALL),
                                         keytype = "GOID", columns = "TERM")
      entrez2go <- merge(entrez2go, as.data.table(go_terms),
                         by.x = "GOALL", by.y = "GOID", all.x = TRUE)
      # Reverse-map back to symbols
      entrez_to_sym <- data.table(ENTREZID = unname(sym2entrez), gene = names(sym2entrez))
      entrez2go <- merge(entrez2go, entrez_to_sym, by = "ENTREZID", all.x = TRUE)
      go_dt <- entrez2go[, .(go_bp = paste(head(unique(TERM), 3), collapse = "; ")),
                         by = gene]
      saveRDS(go_dt, go_cache)
      message("  org.Mm.eg.db fallback succeeded — cached")
      TRUE
    }, error = function(e) {
      message(sprintf("  org.Mm.eg.db fallback also failed: %s", conditionMessage(e)))
      FALSE
    })

    if (!ok2) {
      message("  GO BP layer skipped (no annotation source available)")
      go_dt <- data.table(gene = character(), go_bp = character())
    }
  }
}

if (!is.null(go_dt) && nrow(go_dt) > 0) {
  ino <- merge(ino, go_dt, by = "gene", all.x = TRUE, sort = FALSE)
  message(sprintf("  annotated %d / %d isoforms with GO BP",
                  sum(!is.na(ino$go_bp)), nrow(ino)))
} else {
  ino[, go_bp := NA_character_]
}

# -----------------------------------------------------------------------------
# Layer 3: SynGO quick-list (synaptic protein flag) — local, always works
# -----------------------------------------------------------------------------
message("\nLayer 3: SynGO synaptic-protein flag (quick-list)")
SYNGO_QUICK_LIST <- c(
  # Pre-synaptic
  "Slc32a1", "Slc17a6", "Slc17a7", "Slc6a1", "Syt1", "Syt2", "Snap25", "Stx1a",
  "Vamp1", "Vamp2", "Syp", "Synj1", "Pclo", "Bsn", "Rims1", "Rims2",
  "Gad1", "Gad2", "Pvalb", "Sst", "Vip", "Lamp5", "Sncg", "Reln", "Cck", "Npy",
  "Calb1", "Calb2",
  # Post-synaptic
  "Dlg1", "Dlg4", "Shank1", "Shank2", "Shank3", "Homer1", "Homer2", "Homer3",
  "Gria1", "Gria2", "Gria3", "Gria4", "Grin1", "Grin2a", "Grin2b", "Grin2c",
  "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabra5", "Gabrb1", "Gabrb2", "Gabrb3",
  "Camk2a", "Camk2b",
  # Synaptic adhesion
  "Nlgn1", "Nlgn2", "Nlgn3", "Nrxn1", "Nrxn2", "Nrxn3",
  "Cdh2", "Cdh4", "Cdh10", "Lrrtm1", "Lrrtm2",
  # Endocytic / vesicle cycling / scaffolding
  "Sgip1", "Stxbp1", "Stxbp5", "Ank3", "Ank2", "Ankrd11",
  # Channels
  "Kcnq2", "Kcnq3", "Kcnq5", "Cacna1a", "Cacna1b", "Cacna1c"
)
ino[, syngo_quick := gene %in% SYNGO_QUICK_LIST]
message(sprintf("  flagged %d isoforms as synaptic", sum(ino$syngo_quick, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# Save final table
# -----------------------------------------------------------------------------
out_path <- file.path(RESULTS_DIR, "interneuron_isoforms_annotated.tsv")
fwrite(ino, out_path, sep = "\t")
message(sprintf("\nWrote final annotated table: %s", out_path))
message(sprintf("Columns (%d):", ncol(ino)))
for (col in names(ino)) message(sprintf("  %s", col))

# Summary
message("\nAnnotation coverage:")
message(sprintf("  with biotype:     %d (%.0f%%)",
                sum(!is.na(ino$biotype)), 100 * sum(!is.na(ino$biotype))/nrow(ino)))
message(sprintf("  with GO BP:       %d (%.0f%%)",
                sum(!is.na(ino$go_bp)), 100 * sum(!is.na(ino$go_bp))/nrow(ino)))
message(sprintf("  SynGO quick-flag: %d (%.0f%%)",
                sum(ino$syngo_quick, na.rm = TRUE),
                100 * sum(ino$syngo_quick, na.rm = TRUE)/nrow(ino)))

if ("biotype" %in% names(ino) && sum(!is.na(ino$biotype)) > 0) {
  message("\nBiotype distribution:")
  print(table(ino$biotype, useNA = "ifany"))
}

message("\n05e_interneuron_functional_annotation.R complete")
