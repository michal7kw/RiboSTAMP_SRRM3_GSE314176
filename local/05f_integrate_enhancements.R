# =============================================================================
# 05f_integrate_enhancements.R — merge all Week 2 enhancements into one table
#
# Inputs:
#   results/tables/interneuron_isoforms.tsv               — original (16 cols)
#   results/tables/interneuron_isoforms_by_subtype.tsv    — adds subtype enrichment
#   results/tables/interneuron_isoforms_de.tsv            — adds DESeq2 DE
#   results/tables/interneuron_isoforms_annotated.tsv     — adds biotype/GO/SynGO
#   results/tables/interneuron_isoforms_siblings.tsv      — adds sibling-pair cols
#                                                            (built by 05g_)
#
# Output:
#   results/tables/interneuron_isoforms_FINAL.tsv         — merged, all columns
#
# This is the canonical Linda's deliverable after all Week 2 improvements.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

SCRIPT_DIR <- tryCatch(here::here(), error = function(e) getwd())
PROJ_ROOT <- if (basename(SCRIPT_DIR) == "local") dirname(SCRIPT_DIR) else SCRIPT_DIR
RESULTS_DIR <- file.path(PROJ_ROOT, "results", "tables")

base_path <- file.path(RESULTS_DIR, "interneuron_isoforms.tsv")
subtype_path <- file.path(RESULTS_DIR, "interneuron_isoforms_by_subtype.tsv")
de_path <- file.path(RESULTS_DIR, "interneuron_isoforms_de.tsv")
ann_path <- file.path(RESULTS_DIR, "interneuron_isoforms_annotated.tsv")
sibling_path <- file.path(RESULTS_DIR, "interneuron_isoforms_siblings.tsv")

# Start from the most-enriched layer that includes the original columns
# (the annotated table from 05e includes DE + biotype/GO/SynGO).
if (file.exists(ann_path)) {
  message(sprintf("Loading annotated table: %s", ann_path))
  final_dt <- fread(ann_path)
} else if (file.exists(de_path)) {
  message(sprintf("Loading DE table: %s", de_path))
  final_dt <- fread(de_path)
} else {
  stop("No interneuron table found")
}
message(sprintf("Starting columns: %d", ncol(final_dt)))

# Add per-subtype columns from 05c output (if not already there)
if (file.exists(subtype_path)) {
  message(sprintf("\nMerging subtype enrichment from %s", subtype_path))
  subtype_dt <- fread(subtype_path)
  subtype_cols <- grep("(Pvalb|Sst|Vip|Lamp5|Calb1|Calb2|Cck|Reln|Npy|Other)_(pct|mean)$",
                       names(subtype_dt), value = TRUE)
  if (length(subtype_cols) > 0) {
    keep <- c("isoform_id", subtype_cols)
    subtype_to_merge <- subtype_dt[, ..keep]
    # Drop any of these columns from final_dt if already present, then merge
    overlap <- intersect(names(final_dt), subtype_cols)
    if (length(overlap) > 0) {
      final_dt <- final_dt[, !overlap, with = FALSE]
    }
    final_dt <- merge(final_dt, subtype_to_merge, by = "isoform_id",
                      all.x = TRUE, sort = FALSE)
    message(sprintf("  added %d subtype columns: %s",
                    length(subtype_cols), paste(subtype_cols, collapse = ", ")))
  }
}

# -----------------------------------------------------------------------------
# Add sibling-pair columns from 05g_ output
# -----------------------------------------------------------------------------
if (file.exists(sibling_path)) {
  message(sprintf("\nMerging sibling-pair columns from %s", sibling_path))
  sibling_dt <- fread(sibling_path)
  sibling_cols <- c("sibling_isoform_id", "sibling_GABA_pct", "sibling_GABA_mean",
                    "sibling_other_pct", "sibling_other_mean", "sibling_logFC",
                    "sibling_biotype", "sibling_editsC_GABA_pct",
                    "sibling_editsC_is_translated", "delta_log2_ratio",
                    "pairing_status")
  sibling_cols <- intersect(sibling_cols, names(sibling_dt))
  if (length(sibling_cols) > 0) {
    keep_s <- c("isoform_id", sibling_cols)
    sibling_to_merge <- sibling_dt[, ..keep_s]
    overlap <- intersect(names(final_dt), sibling_cols)
    if (length(overlap) > 0) {
      final_dt <- final_dt[, !overlap, with = FALSE]
    }
    final_dt <- merge(final_dt, sibling_to_merge, by = "isoform_id",
                      all.x = TRUE, sort = FALSE)
    message(sprintf("  added %d sibling columns", length(sibling_cols)))
  }
}

# Reorder columns for readability
preferred_order <- c(
  # Identification
  "isoform_id", "gene", "transcript_id",
  # Primary GABA enrichment (heuristic)
  "GABA_pct", "GABA_mean", "other_pct", "other_mean", "logFC", "n_GABA_detected",
  # Leiden subtype
  "GABA_1_pct", "GABA_1_mean", "GABA_2_pct", "GABA_2_mean",
  # Translation (Ribo-STAMP)
  "editsC_GABA_mean", "editsC_GABA_pct", "editsC_is_translated",
  # Formal DE statistics (DESeq2 pseudobulk)
  "DE_baseMean", "DE_log2FC", "DE_lfcSE", "DE_stat", "DE_pvalue", "DE_padj",
  # Biological subtype enrichment (hippocampus-relevant markers)
  "Pvalb_pct", "Pvalb_mean", "Sst_pct", "Sst_mean",
  "Vip_pct", "Vip_mean", "Lamp5_pct", "Lamp5_mean",
  "Calb1_pct", "Calb1_mean", "Calb2_pct", "Calb2_mean",
  "Cck_pct", "Cck_mean", "Reln_pct", "Reln_mean",
  "Npy_pct", "Npy_mean", "Other_pct", "Other_mean",
  # Functional annotation
  "biotype", "go_bp", "syngo_quick",
  # Sibling-pair (for IN-ribotag-vs-bulk RT-PCR ratio assays — built by 05g_)
  "sibling_isoform_id", "sibling_GABA_pct", "sibling_GABA_mean",
  "sibling_other_pct", "sibling_other_mean", "sibling_logFC",
  "sibling_biotype", "sibling_editsC_GABA_pct",
  "sibling_editsC_is_translated", "delta_log2_ratio", "pairing_status"
)
keep <- intersect(preferred_order, names(final_dt))
extras <- setdiff(names(final_dt), keep)
final_dt <- final_dt[, c(keep, extras), with = FALSE]

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
out_path <- file.path(RESULTS_DIR, "interneuron_isoforms_FINAL.tsv")
fwrite(final_dt, out_path, sep = "\t")
message(sprintf("\nWrote: %s", out_path))
message(sprintf("Final table: %d rows × %d columns", nrow(final_dt), ncol(final_dt)))

# -----------------------------------------------------------------------------
# Print high-confidence subset summaries for the PI
# -----------------------------------------------------------------------------
message("\n=== Summary: high-confidence subsets ===")

# Subset A: heuristic + DE + translation (the strongest)
n_de_sig <- sum(final_dt$DE_padj < 0.05 & final_dt$DE_log2FC > 1, na.rm = TRUE)
message(sprintf("  DESeq2-significant (padj<0.05, log2FC>1):  %d / 274", n_de_sig))

n_de_translated <- sum(final_dt$DE_padj < 0.05 & final_dt$DE_log2FC > 1 &
                       final_dt$editsC_is_translated == TRUE, na.rm = TRUE)
message(sprintf("  ... AND translated (Ribo-STAMP TRUE):       %d", n_de_translated))

# Subset B: subtype-specific (hippocampus markers)
for (st in c("Calb1", "Calb2", "Cck", "Reln", "Vip")) {
  pct_col <- paste0(st, "_pct")
  if (pct_col %in% names(final_dt)) {
    n_st <- sum(final_dt[[pct_col]] >= 30, na.rm = TRUE)
    message(sprintf("  Enriched in %s+ cells (%s ≥ 30%%):         %d",
                    st, pct_col, n_st))
  }
}

# Subset C: protein-coding + translated
n_pc_translated <- sum(final_dt$biotype == "protein_coding" &
                       final_dt$editsC_is_translated == TRUE, na.rm = TRUE)
message(sprintf("  Protein-coding AND translated:              %d", n_pc_translated))

# Subset D: synaptic + translated (highest-confidence wet-lab targets)
n_syn_translated <- sum(final_dt$syngo_quick == TRUE &
                        final_dt$editsC_is_translated == TRUE, na.rm = TRUE)
message(sprintf("  Synaptic flag AND translated:               %d", n_syn_translated))

# Subset E: sibling-pair available (assay-ready for IN-ribotag vs bulk)
if ("pairing_status" %in% names(final_dt)) {
  n_paired <- sum(final_dt$pairing_status == "paired", na.rm = TRUE)
  n_paired_strong <- sum(final_dt$pairing_status == "paired" &
                         final_dt$delta_log2_ratio > 2, na.rm = TRUE)
  message(sprintf("  With sibling pair (assay-ready):            %d", n_paired))
  message(sprintf("  ... AND delta_log2_ratio > 2 (strong swing): %d", n_paired_strong))
}

message("\nFinal table is at: results/tables/interneuron_isoforms_FINAL.tsv")
message("This is the canonical Linda's deliverable with all Week 2 enhancements.")
