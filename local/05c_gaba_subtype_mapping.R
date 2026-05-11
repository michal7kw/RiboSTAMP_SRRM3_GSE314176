# =============================================================================
# 05c_gaba_subtype_mapping.R — map GABA Leiden sub-clusters to classical
# Pvalb/Sst/Vip/Lamp5 functional subtypes
#
# The published long-read clustering provides only 2 GABA Leiden sub-clusters
# (GABA_1, GABA_2). For Linda's wet-lab plans, classical 4-marker subtypes
# (Pvalb, Sst, Vip, Lamp5) are more useful. This script:
#
# 1. Extracts marker-gene expression for each GABA cell from the long-read
#    counts matrix
# 2. Assigns each cell to its dominant subtype based on marker expression
# 3. Reports the % composition of each Leiden cluster by classical subtype
# 4. Emits per-subtype enrichment columns for the interneuron table
#
# Outputs:
#   results/tables/gaba_subtype_mapping.tsv     — per-cell subtype assignment
#   results/tables/gaba_subtype_composition.tsv — Leiden cluster × subtype matrix
#   results/tables/interneuron_isoforms_by_subtype.tsv — enhanced isoform table
#     with per-subtype Pvalb/Sst/Vip/Lamp5 enrichment columns
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(here)
})

# Resolve project root regardless of how we're called
SCRIPT_DIR <- tryCatch(here::here(), error = function(e) getwd())
PROJ_ROOT <- if (basename(SCRIPT_DIR) == "local") dirname(SCRIPT_DIR) else SCRIPT_DIR
META_DIR <- file.path(PROJ_ROOT, "data", "metadata")
RESULTS_DIR <- file.path(PROJ_ROOT, "results", "tables")
CACHE_DIR <- file.path(PROJ_ROOT, "cache")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
# Subtype markers — extended for hippocampal interneurons. Cortical
# Pvalb/Sst/Vip/Lamp5 typology fits hippocampus poorly because most
# hippocampal GABA cells are calbindin/calretinin/CCK/reelin-positive
# (basket, bistratified, O-LM, axo-axonic, Cajal-Retzius derivatives).
# Including the hippocampus-relevant markers (Calb1, Calb2, Cck, Reln,
# Npy) gives meaningful subtype assignment.
SUBTYPE_MARKERS <- list(
  Pvalb  = c("Pvalb"),                  # fast-spiking basket / axo-axonic
  Sst    = c("Sst"),                    # O-LM, bistratified
  Vip    = c("Vip"),                    # VIP+ interneurons
  Lamp5  = c("Lamp5"),                  # neurogliaform (cortex; rare in hippo)
  Calb1  = c("Calb1"),                  # calbindin — common in hippocampus
  Calb2  = c("Calb2"),                  # calretinin — common in hippocampus
  Cck    = c("Cck"),                    # CCK basket cells (hippocampus)
  Reln   = c("Reln"),                   # Reelin+ (Cajal-Retzius derivatives)
  Npy    = c("Npy")                     # neuropeptide Y (sometimes co-expressed)
)

# To assign a cell to a subtype, the marker count must be ≥ this threshold.
# Long-read sc-RNAseq has very low coverage per cell (~1-5 K reads), so
# requiring >1 transcript count of a marker is too strict for sparse markers.
# Use 1 (any non-zero detection) by default; adjust if too many "Other".
MIN_MARKER_COUNT <- 1

# Three mice (long-read GSMs)
MICE_GSM <- c("GSM9380799", "GSM9380800", "GSM9380801")

# -----------------------------------------------------------------------------
# Helper: load long-read counts matrix + obs for a mouse
# -----------------------------------------------------------------------------
load_lr_data <- function(gsm) {
  message(sprintf("  loading %s ...", gsm))
  obs_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*_obs\\.tsv$", gsm),
                         full.names = TRUE)
  var_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*_var\\.tsv$", gsm),
                         full.names = TRUE)
  mtx_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*\\.mtx$", gsm),
                         full.names = TRUE)
  if (length(obs_path) == 0 || length(var_path) == 0 || length(mtx_path) == 0) {
    stop(sprintf("Missing files for %s in %s", gsm, META_DIR))
  }
  obs <- fread(obs_path[1])
  var <- fread(var_path[1])
  mtx <- Matrix::readMM(mtx_path[1])
  # Conventional: rows = features, cols = cells in our exploded h5ad
  if (nrow(mtx) == nrow(obs) && ncol(mtx) == nrow(var)) {
    mtx <- t(mtx)  # transpose to features × cells
  }
  rownames(mtx) <- var[[1]]
  colnames(mtx) <- obs[[1]]
  list(mtx = mtx, obs = obs, var = var, gsm = gsm)
}

# -----------------------------------------------------------------------------
# Helper: aggregate transcript counts to gene level for marker genes
# -----------------------------------------------------------------------------
sum_isoforms_to_gene <- function(mtx, var, gene_symbols) {
  if (!"gene_name" %in% names(var) && !"gene" %in% names(var)) {
    stop("var.tsv must have a 'gene_name' or 'gene' column")
  }
  gene_col <- if ("gene_name" %in% names(var)) "gene_name" else "gene"
  out <- matrix(0, nrow = length(gene_symbols), ncol = ncol(mtx),
                dimnames = list(gene_symbols, colnames(mtx)))
  for (g in gene_symbols) {
    rows <- which(var[[gene_col]] == g)
    if (length(rows) == 0) next
    out[g, ] <- Matrix::colSums(mtx[rows, , drop = FALSE])
  }
  out
}

# -----------------------------------------------------------------------------
# Step 1: Load all 3 mice and pool GABA cells with marker expression
# -----------------------------------------------------------------------------
message("Step 1: loading and pooling GABA cells across 3 mice")
all_marker_expr <- list()
all_obs <- list()

for (gsm in MICE_GSM) {
  d <- load_lr_data(gsm)
  # Filter to GABA cells
  gaba_mask <- d$obs[["Cell Assignments Grouped"]] == "Neuron - GABA"
  message(sprintf("    GABA cells: %d / %d", sum(gaba_mask, na.rm = TRUE), nrow(d$obs)))
  gaba_obs <- d$obs[gaba_mask, ]
  gaba_mtx <- d$mtx[, gaba_mask, drop = FALSE]

  # Extract marker expression (sum at gene level)
  marker_genes <- unique(unlist(SUBTYPE_MARKERS))
  marker_expr <- sum_isoforms_to_gene(gaba_mtx, d$var, marker_genes)

  all_marker_expr[[gsm]] <- as.data.frame(t(marker_expr))
  all_marker_expr[[gsm]]$mouse <- gsm
  all_marker_expr[[gsm]]$cb <- colnames(marker_expr)
  all_marker_expr[[gsm]]$leiden_subtype <- gaba_obs[["Cell Assignment"]]

  all_obs[[gsm]] <- gaba_obs[, c("barcode", "Cell Assignment", "Cell Assignments Grouped")]
  all_obs[[gsm]][, mouse := gsm]
}

marker_long <- rbindlist(lapply(all_marker_expr, as.data.table), fill = TRUE)
message(sprintf("Total GABA cells pooled: %d", nrow(marker_long)))

# -----------------------------------------------------------------------------
# Step 2: Assign each cell to dominant subtype
# -----------------------------------------------------------------------------
message("\nStep 2: assigning cells to dominant subtype")

# A cell's "subtype" is the marker class with the highest expression, IF that
# expression is ≥ MIN_MARKER_COUNT. Cells with all markers below threshold are
# labeled "Other" (typically rare subtypes like Sncg, Meis2, Igtp+).
assign_subtype <- function(row) {
  marker_vals <- list()
  for (st in names(SUBTYPE_MARKERS)) {
    markers <- SUBTYPE_MARKERS[[st]]
    vals <- as.numeric(row[markers])
    marker_vals[[st]] <- max(vals, na.rm = TRUE)
  }
  marker_vals <- unlist(marker_vals)
  best <- which.max(marker_vals)
  if (length(best) == 0 || marker_vals[best] < MIN_MARKER_COUNT) {
    return("Other")
  }
  names(SUBTYPE_MARKERS)[best]
}

marker_long[, classical_subtype := apply(.SD, 1, assign_subtype),
            .SDcols = unique(unlist(SUBTYPE_MARKERS))]

# Save per-cell mapping
fwrite(marker_long[, .(mouse, cb, leiden_subtype, classical_subtype,
                       Pvalb, Sst, Vip, Lamp5)],
       file.path(RESULTS_DIR, "gaba_subtype_mapping.tsv"), sep = "\t")
message(sprintf("  wrote per-cell mapping: %s", file.path(RESULTS_DIR, "gaba_subtype_mapping.tsv")))

# -----------------------------------------------------------------------------
# Step 3: Composition matrix (Leiden cluster × classical subtype)
# -----------------------------------------------------------------------------
message("\nStep 3: composition (Leiden × classical subtype)")
composition <- marker_long[, .N, by = .(leiden_subtype, classical_subtype)]
composition_wide <- dcast(composition, leiden_subtype ~ classical_subtype, value.var = "N", fill = 0)
# Add row totals + percentages
n_per_leiden <- composition_wide[, total := rowSums(.SD), .SDcols = setdiff(names(composition_wide), "leiden_subtype")]
composition_pct <- copy(composition_wide)
for (col in setdiff(names(composition_pct), c("leiden_subtype", "total"))) {
  composition_pct[[paste0(col, "_pct")]] <- round(100 * composition_wide[[col]] / composition_wide$total, 1)
}

print(composition_wide)
fwrite(composition_pct, file.path(RESULTS_DIR, "gaba_subtype_composition.tsv"), sep = "\t")
message(sprintf("  wrote composition table: %s",
                file.path(RESULTS_DIR, "gaba_subtype_composition.tsv")))

# -----------------------------------------------------------------------------
# Step 4: Augment the existing interneuron_isoforms.tsv with per-subtype enrichment
# -----------------------------------------------------------------------------
message("\nStep 4: per-subtype isoform enrichment")
ino_path <- file.path(RESULTS_DIR, "interneuron_isoforms.tsv")
if (!file.exists(ino_path)) {
  message(sprintf("  WARN: %s not found; skipping per-subtype enrichment join", ino_path))
} else {
  message(sprintf("  loading existing interneuron_isoforms table: %s", ino_path))
  ino <- fread(ino_path)

  # For each isoform, compute mean expression in each classical subtype
  # Need to load the full isoform matrix + the per-cell subtype assignments
  # We already have marker_long with classical_subtype per cell

  # Build a subtype lookup: (mouse, cb) → classical_subtype
  subtype_lookup <- marker_long[, .(mouse, cb, classical_subtype)]
  setkey(subtype_lookup, mouse, cb)

  per_subtype_pct <- list()
  per_subtype_mean <- list()

  for (gsm in MICE_GSM) {
    d <- load_lr_data(gsm)
    obs <- d$obs
    obs[, mouse := gsm]

    gaba_mask <- obs[["Cell Assignments Grouped"]] == "Neuron - GABA"
    if (sum(gaba_mask, na.rm = TRUE) == 0) next

    gaba_cb <- obs$barcode[gaba_mask]
    gaba_mtx <- d$mtx[, gaba_mask, drop = FALSE]

    # Lookup subtype for each GABA cell
    cb_keys <- data.table(mouse = gsm, cb = gaba_cb)
    cb_keys <- merge(cb_keys, subtype_lookup, by = c("mouse", "cb"),
                     all.x = TRUE, sort = FALSE)
    cb_keys[is.na(classical_subtype), classical_subtype := "Unknown"]

    # For each subtype, compute pct + mean per isoform
    for (st in c(names(SUBTYPE_MARKERS), "Other")) {
      st_cells <- cb_keys$classical_subtype == st
      n_st <- sum(st_cells, na.rm = TRUE)
      if (n_st == 0) next
      sub_mtx <- gaba_mtx[, st_cells, drop = FALSE]
      pct <- Matrix::rowSums(sub_mtx > 0) / n_st * 100
      meanval <- Matrix::rowMeans(sub_mtx)
      key <- paste0(gsm, "_", st)
      per_subtype_pct[[key]] <- pct
      per_subtype_mean[[key]] <- meanval
    }
  }

  message("  combining across mice (sum then average)")
  # Combine across mice for each subtype: sum cells, average rates
  combined_pct <- list()
  combined_mean <- list()
  for (st in c(names(SUBTYPE_MARKERS), "Other")) {
    keys_st <- grep(paste0("_", st, "$"), names(per_subtype_pct), value = TRUE)
    if (length(keys_st) == 0) next
    pct_mat <- do.call(cbind, per_subtype_pct[keys_st])
    mean_mat <- do.call(cbind, per_subtype_mean[keys_st])
    # Mean across mice (could weight by cell count, but simpler is fine for ranking)
    combined_pct[[st]] <- rowMeans(pct_mat, na.rm = TRUE)
    combined_mean[[st]] <- rowMeans(mean_mat, na.rm = TRUE)
  }

  # Build augmented table
  message("  augmenting interneuron_isoforms table with per-subtype columns")
  augmented <- copy(ino)
  for (st in names(combined_pct)) {
    pct_col <- paste0(st, "_pct")
    mean_col <- paste0(st, "_mean")
    # Note: rownames of count matrix are feature_id (e.g. "Slc32a1-201"),
    # which corresponds to the isoform_id column of the interneuron table
    augmented[[pct_col]] <- combined_pct[[st]][augmented$isoform_id]
    augmented[[mean_col]] <- combined_mean[[st]][augmented$isoform_id]
  }
  # Replace NAs with 0 (isoform not detected in that subtype)
  for (col in grep("(Pvalb|Sst|Vip|Lamp5|Other)_(pct|mean)$", names(augmented), value = TRUE)) {
    set(augmented, which(is.na(augmented[[col]])), col, 0)
  }

  out_path <- file.path(RESULTS_DIR, "interneuron_isoforms_by_subtype.tsv")
  fwrite(augmented, out_path, sep = "\t")
  message(sprintf("  wrote enhanced table: %s", out_path))

  # Summary stats per subtype: how many isoforms enriched in each?
  message("\nPer-subtype isoform counts (pct ≥ 30% AND mean ≥ 1):")
  for (st in names(SUBTYPE_MARKERS)) {
    pct_col <- paste0(st, "_pct")
    mean_col <- paste0(st, "_mean")
    n_enriched <- sum(augmented[[pct_col]] >= 30 & augmented[[mean_col]] >= 1, na.rm = TRUE)
    message(sprintf("  %-8s: %d isoforms", st, n_enriched))
  }
}

message("\n05c_gaba_subtype_mapping.R complete")
