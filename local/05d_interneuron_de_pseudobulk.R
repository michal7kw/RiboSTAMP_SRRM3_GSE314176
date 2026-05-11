# =============================================================================
# 05d_interneuron_de_pseudobulk.R — formal DE testing for the interneuron table
#
# Replaces the heuristic thresholds (GABA_pct ≥ 10, logFC ≥ 1) with
# DESeq2 pseudobulk DE: aggregate counts per (mouse × cell_type) into
# pseudobulk samples, then run a standard GABA-vs-other DE test with
# Wald tests + BH multiple-testing correction.
#
# This adds FDR-adjusted p-values to the interneuron table, making it
# publication-defensible.
#
# Outputs:
#   results/tables/interneuron_isoforms_de.tsv  — same isoforms + DE columns
#     (log2FC_DESeq2, baseMean, pvalue, padj, etc.)
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(DESeq2)
  library(here)
})

SCRIPT_DIR <- tryCatch(here::here(), error = function(e) getwd())
PROJ_ROOT <- if (basename(SCRIPT_DIR) == "local") dirname(SCRIPT_DIR) else SCRIPT_DIR
META_DIR <- file.path(PROJ_ROOT, "data", "metadata")
RESULTS_DIR <- file.path(PROJ_ROOT, "results", "tables")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

MICE_GSM <- c("GSM9380799", "GSM9380800", "GSM9380801")

# -----------------------------------------------------------------------------
# Load data — pool counts across mice
# -----------------------------------------------------------------------------
message("Loading long-read counts matrices for 3 mice")
all_mtx <- list()
all_obs <- list()
features <- NULL

for (gsm in MICE_GSM) {
  message(sprintf("  %s ...", gsm))
  obs_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*_obs\\.tsv$", gsm),
                         full.names = TRUE)[1]
  var_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*_var\\.tsv$", gsm),
                         full.names = TRUE)[1]
  mtx_path <- list.files(META_DIR,
                         pattern = sprintf("^%s_longread_normed_counts.*\\.mtx$", gsm),
                         full.names = TRUE)[1]
  if (is.na(obs_path) || is.na(mtx_path)) {
    stop(sprintf("Missing files for %s", gsm))
  }
  obs <- fread(obs_path)
  var <- fread(var_path)
  mtx <- Matrix::readMM(mtx_path)
  if (nrow(mtx) == nrow(obs) && ncol(mtx) == nrow(var)) {
    mtx <- t(mtx)
  }
  rownames(mtx) <- var[[1]]
  colnames(mtx) <- obs[[1]]
  obs[, mouse := gsm]
  all_mtx[[gsm]] <- mtx
  all_obs[[gsm]] <- obs
  if (is.null(features)) features <- var
}

# Pool. We assume the feature universe (rownames) is identical across mice;
# verify here.
ref_features <- rownames(all_mtx[[1]])
for (gsm in MICE_GSM) {
  if (!identical(rownames(all_mtx[[gsm]]), ref_features)) {
    stop(sprintf("Feature mismatch for %s", gsm))
  }
}

# -----------------------------------------------------------------------------
# Pseudobulk: aggregate counts per (mouse × cell_type) → 21 pseudobulk samples
#   3 mice × 7 cell types = 21 (or fewer if some cell types are absent in a mouse)
# -----------------------------------------------------------------------------
message("\nBuilding pseudobulk count matrix")
pseudobulk_mat <- list()
pseudobulk_meta <- list()

for (gsm in MICE_GSM) {
  obs <- all_obs[[gsm]]
  mtx <- all_mtx[[gsm]]
  cell_types <- unique(obs[["Cell Assignments Grouped"]])
  for (ct in cell_types) {
    if (is.na(ct)) next
    cells_ct <- which(obs[["Cell Assignments Grouped"]] == ct)
    if (length(cells_ct) < 5) next  # too few cells for a pseudobulk
    sample_id <- sprintf("%s__%s", gsm, gsub("[^A-Za-z0-9]", "_", ct))
    pseudobulk_mat[[sample_id]] <- Matrix::rowSums(mtx[, cells_ct, drop = FALSE])
    pseudobulk_meta[[sample_id]] <- data.frame(
      sample_id = sample_id,
      mouse = gsm,
      cell_type = ct,
      n_cells = length(cells_ct),
      is_gaba = ct == "Neuron - GABA"
    )
  }
}
pb_counts <- do.call(cbind, pseudobulk_mat)
pb_meta <- do.call(rbind, pseudobulk_meta)
rownames(pb_meta) <- pb_meta$sample_id

message(sprintf("  pseudobulk matrix: %d features × %d samples",
                nrow(pb_counts), ncol(pb_counts)))
message(sprintf("  GABA pseudobulk samples: %d", sum(pb_meta$is_gaba)))
message(sprintf("  non-GABA pseudobulk samples: %d", sum(!pb_meta$is_gaba)))

# Round counts to integers (DESeq2 requires integers)
pb_counts_int <- round(as.matrix(pb_counts))
mode(pb_counts_int) <- "integer"

# -----------------------------------------------------------------------------
# DESeq2 — GABA vs other
# -----------------------------------------------------------------------------
message("\nRunning DESeq2 (GABA vs other)")
pb_meta$is_gaba <- factor(pb_meta$is_gaba, levels = c(FALSE, TRUE))
pb_meta$mouse <- factor(pb_meta$mouse)

dds <- DESeqDataSetFromMatrix(
  countData = pb_counts_int,
  colData = pb_meta,
  design = ~ mouse + is_gaba   # control for mouse effects
)
# Filter: features with ≥ 10 counts total (most very-low isoforms aren't testable)
dds <- dds[rowSums(counts(dds)) >= 10, ]
message(sprintf("  features kept after low-count filter: %d", nrow(dds)))

dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, name = "is_gaba_TRUE_vs_FALSE", alpha = 0.05)

message(sprintf("  significant (padj < 0.05): %d", sum(res$padj < 0.05, na.rm = TRUE)))
message(sprintf("  significant + log2FC > 1: %d",
                sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)))

de_dt <- as.data.table(as.data.frame(res), keep.rownames = "isoform_id")
setnames(de_dt,
         c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
         c("DE_baseMean", "DE_log2FC", "DE_lfcSE", "DE_stat", "DE_pvalue", "DE_padj"))

# -----------------------------------------------------------------------------
# Join with existing interneuron_isoforms.tsv
# -----------------------------------------------------------------------------
ino_path <- file.path(RESULTS_DIR, "interneuron_isoforms.tsv")
if (file.exists(ino_path)) {
  message(sprintf("\nJoining with %s", ino_path))
  ino <- fread(ino_path)
  ino <- merge(ino, de_dt, by = "isoform_id", all.x = TRUE, sort = FALSE)
  out_path <- file.path(RESULTS_DIR, "interneuron_isoforms_de.tsv")
  fwrite(ino, out_path, sep = "\t")
  message(sprintf("  wrote: %s", out_path))

  message("\nSummary of DE columns added:")
  message(sprintf("  with DE_padj available: %d / %d isoforms",
                  sum(!is.na(ino$DE_padj)), nrow(ino)))
  message(sprintf("  passing DE_padj < 0.05 AND DE_log2FC > 1: %d",
                  sum(ino$DE_padj < 0.05 & ino$DE_log2FC > 1, na.rm = TRUE)))
  message(sprintf("  passing DE_padj < 0.05 (any direction): %d",
                  sum(ino$DE_padj < 0.05, na.rm = TRUE)))
} else {
  message(sprintf("\nWARN: %s not found — saving DE-only table instead", ino_path))
  out_path <- file.path(RESULTS_DIR, "interneuron_de_results.tsv")
  fwrite(de_dt, out_path, sep = "\t")
  message(sprintf("  wrote: %s", out_path))
}

# Save the dds object for reproducibility / further analysis
saveRDS(dds, file.path(PROJ_ROOT, "cache", "interneuron_dds.rds"))
message(sprintf("  cached DESeq2 object: %s",
                file.path(PROJ_ROOT, "cache", "interneuron_dds.rds")))

message("\n05d_interneuron_de_pseudobulk.R complete")
