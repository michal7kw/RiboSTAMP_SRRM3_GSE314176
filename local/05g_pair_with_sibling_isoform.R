#####################################################################
# 05g_pair_with_sibling_isoform.R
#####################################################################
# Pair each GABA-enriched isoform (from interneuron_isoforms_FINAL.tsv)
# with its "normal" non-GABA sibling — a different isoform of the same
# gene that is NOT GABA-enriched — for use in IN-ribotag-vs-bulk RT-PCR
# / qPCR ratio assays.
#
# Sibling-selection rule: same gene, NOT in the 274 GABA-enriched set,
#   logFC < 1 (i.e. not GABA-enriched), other_mean >= 0.001 (detectable
#   in non-GABA cells), pick the candidate with highest other_mean.
#
# Headline metric:
#   delta_log2_ratio = log2((GABA_iso_GABA_mean + eps)/(sibling_GABA_mean + eps))
#                    − log2((GABA_iso_other_mean + eps)/(sibling_other_mean + eps))
# Higher = bigger predicted RT-PCR swing between IN-ribotag and bulk
# hippocampus, i.e. better assay candidate.
#
# pairing_status taxonomy (mutually exclusive, exhaustive):
#   paired                       sibling found
#   single_isoform_gene          gene has only one isoform in the catalog
#   all_isoforms_gaba_enriched   every isoform of this gene is in the 274,
#                                  or all non-274 isoforms have logFC >= 1
#   no_expressed_sibling         non-GABA candidates exist but all have
#                                  other_mean < 0.001 (undetected in non-GABA)
#
# Inputs:
#   results/tables/interneuron_isoforms_FINAL.tsv   (the 274 deliverable)
#   cache/pooled_normed_counts.rds                   (built by 05_)
#   cache/pooled_editsc.rds                          (built by 06_)
#   cache/cell_meta.rds                              (built by 01_)
#
# Outputs:
#   results/tables/interneuron_isoforms_siblings.tsv             (274 × 15)
#                                                                full table,
#                                                                empty sibling
#                                                                cells where
#                                                                pairing_status
#                                                                != "paired"
#   results/tables/interneuron_isoforms_siblings_paired_only.tsv ( 59 × 15)
#                                                                no-NA wet-lab
#                                                                subset for
#                                                                assay design
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

# Constants — kept consistent with 05_interneuron_isoform_table.R
# and 06_editsc_translation_column.R so siblings carry the same
# definitions of "translated" and "GABA-enriched" as the 274.
PSEUDO                   <- 1e-6   # logFC pseudo-count (matches 05)
SIBLING_LOGFC_CUTOFF     <- 1      # NOT GABA-enriched (matches 05's MIN_LOGFC)
MIN_EXPRESSED_OTHER_MEAN <- 0.001  # sibling must be detectable in non-GABA
MIN_EDITSC_DETECTION_PCT <- 5      # editsC translation flag (matches 06)

main <- function() {
  final_path <- file.path(TABLES, "interneuron_isoforms_FINAL.tsv")
  if (!file.exists(final_path)) {
    stop("Run 05f_integrate_enhancements.R first — missing ", final_path)
  }
  final <- fread(final_path)
  message(sprintf("Loaded %d GABA-enriched isoforms from %s",
                  nrow(final), final_path))

  meta_path   <- file.path(CACHE_DIR, "cell_meta.rds")
  pooled_path <- file.path(CACHE_DIR, "pooled_normed_counts.rds")
  editsc_path <- file.path(CACHE_DIR, "pooled_editsc.rds")
  if (!file.exists(meta_path))   stop("Run 01_parse_geo_metadata.R first.")
  if (!file.exists(pooled_path)) stop("Run 05_interneuron_isoform_table.R first.")

  meta     <- readRDS(meta_path)
  pooled   <- readRDS(pooled_path)
  counts   <- pooled$counts
  features <- pooled$features

  # Align meta to count matrix columns (mirrors 05_, line 105–109)
  common_bc <- intersect(colnames(counts), meta$barcode)
  counts    <- counts[, common_bc, drop = FALSE]
  meta      <- meta[match(common_bc, meta$barcode)]
  gaba_mask  <- meta$is_gaba
  other_mask <- !meta$is_gaba
  if (sum(gaba_mask) == 0) {
    stop("No GABA cells found — check is_gaba in cell_meta.rds.")
  }

  # ---- Subset to candidate pool: all isoforms of genes-of-interest ----
  goi <- unique(final$gene)
  feat_oi <- features[gene_name %in% goi]
  iso_oi  <- feat_oi$feature_id
  message(sprintf("Genes in 274-deliverable: %d. Candidate-pool size: %d isoforms (of %d).",
                  length(goi), length(iso_oi), nrow(features)))
  counts_oi <- counts[iso_oi, , drop = FALSE]

  gaba_counts  <- counts_oi[, gaba_mask,  drop = FALSE]
  other_counts <- counts_oi[, other_mask, drop = FALSE]

  # Per-isoform stats — same code path as 05_ lines 121–125
  GABA_pct   <- 100 * Matrix::rowSums(gaba_counts  > 0) / ncol(gaba_counts)
  GABA_mean  <- Matrix::rowMeans(gaba_counts)
  other_pct  <- 100 * Matrix::rowSums(other_counts > 0) / ncol(other_counts)
  other_mean <- Matrix::rowMeans(other_counts)
  logFC      <- log2((GABA_mean + PSEUDO) / (other_mean + PSEUDO))

  all_stats <- data.table(
    isoform_id = iso_oi,
    gene       = feat_oi$gene_name,
    biotype    = feat_oi$transcript_biotype,
    GABA_pct   = GABA_pct,
    GABA_mean  = GABA_mean,
    other_pct  = other_pct,
    other_mean = other_mean,
    logFC      = logFC
  )

  # ---- editsC stats for the candidate pool (left-join, NA-tolerant) ----
  if (file.exists(editsc_path)) {
    editsc <- readRDS(editsc_path)
    common_bc_e <- intersect(colnames(editsc), meta$barcode)
    editsc      <- editsc[, common_bc_e, drop = FALSE]
    meta_e      <- meta[match(common_bc_e, meta$barcode)]
    gaba_mask_e <- meta_e$is_gaba
    iso_in_editsc <- intersect(iso_oi, rownames(editsc))
    if (length(iso_in_editsc) > 0) {
      editsc_oi <- editsc[iso_in_editsc, gaba_mask_e, drop = FALSE]
      e_pct  <- 100 * Matrix::rowSums(editsc_oi > 0) / ncol(editsc_oi)
      e_mean <- Matrix::rowMeans(editsc_oi)
      editsc_stats <- data.table(
        isoform_id           = iso_in_editsc,
        editsC_GABA_pct      = e_pct,
        editsC_GABA_mean     = e_mean,
        editsC_is_translated = e_pct >= MIN_EDITSC_DETECTION_PCT
      )
      all_stats <- merge(all_stats, editsc_stats, by = "isoform_id",
                         all.x = TRUE, sort = FALSE)
    } else {
      all_stats[, `:=`(editsC_GABA_pct = NA_real_,
                       editsC_GABA_mean = NA_real_,
                       editsC_is_translated = NA)]
    }
  } else {
    message("pooled_editsc.rds not found — sibling editsC columns will be NA.")
    all_stats[, `:=`(editsC_GABA_pct = NA_real_,
                     editsC_GABA_mean = NA_real_,
                     editsC_is_translated = NA)]
  }

  setkey(all_stats, isoform_id)
  by_gene <- split(all_stats$isoform_id, all_stats$gene)
  in_274  <- final$isoform_id

  pick_sibling_for <- function(gaba_iso, gene_sym) {
    iso_in_gene <- by_gene[[gene_sym]]
    if (is.null(iso_in_gene) || length(iso_in_gene) <= 1L) {
      return(list(status = "single_isoform_gene", sib = NULL))
    }
    cands <- setdiff(iso_in_gene, in_274)
    if (length(cands) == 0L) {
      return(list(status = "all_isoforms_gaba_enriched", sib = NULL))
    }
    cs <- all_stats[J(cands)]
    cs <- cs[logFC < SIBLING_LOGFC_CUTOFF]
    if (nrow(cs) == 0L) {
      return(list(status = "all_isoforms_gaba_enriched", sib = NULL))
    }
    cs <- cs[other_mean >= MIN_EXPRESSED_OTHER_MEAN]
    if (nrow(cs) == 0L) {
      return(list(status = "no_expressed_sibling", sib = NULL))
    }
    setorder(cs, -other_mean)
    list(status = "paired", sib = cs[1])
  }

  # Build sibling table row-by-row (274 rows — fine to loop)
  out_rows <- vector("list", nrow(final))
  for (i in seq_len(nrow(final))) {
    iso  <- final$isoform_id[i]
    gene <- final$gene[i]
    res  <- pick_sibling_for(iso, gene)
    if (res$status == "paired") {
      sib <- res$sib
      g_g <- final$GABA_mean[i]
      g_o <- final$other_mean[i]
      s_g <- sib$GABA_mean
      s_o <- sib$other_mean
      delta <- log2((g_g + PSEUDO) / (s_g + PSEUDO)) -
               log2((g_o + PSEUDO) / (s_o + PSEUDO))
      out_rows[[i]] <- data.table(
        isoform_id                   = iso,
        gene                         = gene,
        GABA_iso_GABA_mean           = g_g,
        GABA_iso_other_mean          = g_o,
        sibling_isoform_id           = sib$isoform_id,
        sibling_GABA_pct             = sib$GABA_pct,
        sibling_GABA_mean            = sib$GABA_mean,
        sibling_other_pct            = sib$other_pct,
        sibling_other_mean           = sib$other_mean,
        sibling_logFC                = sib$logFC,
        sibling_biotype              = sib$biotype,
        sibling_editsC_GABA_pct      = sib$editsC_GABA_pct,
        sibling_editsC_is_translated = sib$editsC_is_translated,
        delta_log2_ratio             = delta,
        pairing_status               = "paired"
      )
    } else {
      out_rows[[i]] <- data.table(
        isoform_id                   = iso,
        gene                         = gene,
        GABA_iso_GABA_mean           = final$GABA_mean[i],
        GABA_iso_other_mean          = final$other_mean[i],
        sibling_isoform_id           = NA_character_,
        sibling_GABA_pct             = NA_real_,
        sibling_GABA_mean            = NA_real_,
        sibling_other_pct            = NA_real_,
        sibling_other_mean           = NA_real_,
        sibling_logFC                = NA_real_,
        sibling_biotype              = NA_character_,
        sibling_editsC_GABA_pct      = NA_real_,
        sibling_editsC_is_translated = NA,
        delta_log2_ratio             = NA_real_,
        pairing_status               = res$status
      )
    }
  }
  siblings_dt <- rbindlist(out_rows)

  # Sort: paired rows first (by delta desc), other status categories follow.
  status_order <- c("paired", "no_expressed_sibling",
                    "all_isoforms_gaba_enriched", "single_isoform_gene")
  siblings_dt[, .pri := match(pairing_status, status_order)]
  setorder(siblings_dt, .pri, -delta_log2_ratio, na.last = TRUE)
  siblings_dt[, .pri := NULL]

  out_path <- file.path(TABLES, "interneuron_isoforms_siblings.tsv")
  fwrite(siblings_dt, out_path, sep = "\t")
  message(sprintf("\nWrote %d rows × %d cols → %s",
                  nrow(siblings_dt), ncol(siblings_dt), out_path))

  # Wet-lab-friendly subset: only paired rows, no NA cells.
  paired_only <- siblings_dt[pairing_status == "paired"]
  paired_path <- file.path(TABLES, "interneuron_isoforms_siblings_paired_only.tsv")
  fwrite(paired_only, paired_path, sep = "\t")
  message(sprintf("Wrote %d rows × %d cols → %s  (paired-only, no NAs)",
                  nrow(paired_only), ncol(paired_only), paired_path))

  message("\npairing_status breakdown:")
  print(siblings_dt[, .N, by = pairing_status])

  paired <- siblings_dt[pairing_status == "paired"]
  if (nrow(paired) > 0) {
    message(sprintf("\ndelta_log2_ratio: median=%.2f, top-quartile cutoff=%.2f",
                    median(paired$delta_log2_ratio),
                    quantile(paired$delta_log2_ratio, 0.75)))
    message(sprintf("Paired rows with delta > 0: %d / %d",
                    sum(paired$delta_log2_ratio > 0), nrow(paired)))
    message(sprintf("Paired rows with translated sibling: %d (sibling also actively translated)",
                    sum(paired$sibling_editsC_is_translated %in% TRUE)))
  }

  invisible(siblings_dt)
}

if (sys.nframe() == 0) main()
