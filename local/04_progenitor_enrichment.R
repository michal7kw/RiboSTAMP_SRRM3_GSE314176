#####################################################################
# 04_progenitor_enrichment.R
#####################################################################
# DELIVERABLE 1 — test whether the novel Srrm3 isoform is enriched in
# progenitor-like cells (OPC + iDG at P25).
#
# Inputs (from prior scripts' caches):
#   cache/cell_meta.rds     # barcode, cell_type, is_progenitor, ...
#   cache/iso_counts.rds    # sparse isoforms x cells
#   data/liftover/*         # mm39-lifted novel isoforms + overlap_with_anchor.bed
#
# Outputs:
#   results/tables/novel_srrm3_per_cluster.tsv
#   results/tables/novel_srrm3_fisher_enrichment.tsv
#   results/figures/novel_srrm3_cluster_psi.pdf
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
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
FIGURES   <- file.path(BASE_DIR, "results", "figures")
LIFT_DIR  <- file.path(BASE_DIR, "data", "liftover")
dir.create(TABLES,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

# Short-read anchor coordinate (mm39) — documented in
# ../90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md
SHORTREAD_ANCHOR <- list(
  chrom = "chr5",
  start = 135898574L,
  end   = 135898652L,
  length_bp = 79L,
  strand = "-"
)

identify_novel_srrm3_isoforms <- function(iso_ids) {
  # Pick candidate novel Srrm3 isoforms: IsoQuant assigns them to a gene_id;
  # we filter to (a) SQANTI3 NIC/NNC classes and (b) proximity to the anchor
  # via data/liftover/overlap_with_anchor.bed (populated by hpc/06_liftover).
  overlap_bed <- file.path(LIFT_DIR, "overlap_with_anchor.bed")
  nearest_bed <- file.path(LIFT_DIR, "nearest_novel_to_anchor.bed")

  if (file.exists(overlap_bed) && file.info(overlap_bed)$size > 0) {
    ov <- fread(overlap_bed, header = FALSE)
    candidate_ids <- unique(ov$V4)
    message(sprintf("Found %d novel isoform(s) directly overlapping the short-read anchor.",
                    length(candidate_ids)))
    attr(candidate_ids, "mode") <- "direct_overlap"
    return(candidate_ids)
  }

  if (file.exists(nearest_bed) && file.info(nearest_bed)$size > 0) {
    nr <- fread(nearest_bed, header = FALSE)
    # last column is distance (bp); keep within 50 bp as plausible same-event
    nr[, distance := .SD[[ncol(nr)]]]
    close <- nr[distance >= 0 & distance <= 50]
    candidate_ids <- unique(close[[10]])  # novel feature name column
    message(sprintf("No direct overlap. %d novel isoform(s) within 50 bp of anchor.",
                    length(candidate_ids)))
    attr(candidate_ids, "mode") <- "nearest_within_50bp"
    return(candidate_ids)
  }

  warning("No liftOver outputs found under ", LIFT_DIR,
          " — run hpc/08_liftover_mm10_to_mm39.sh first.")
  character(0)
}

main <- function() {
  meta   <- readRDS(file.path(CACHE_DIR, "cell_meta.rds"))
  counts <- readRDS(file.path(CACHE_DIR, "iso_counts.rds"))

  # Align meta order to counts columns
  common_bc <- intersect(colnames(counts), meta$barcode)
  if (length(common_bc) == 0) {
    stop("No overlap between barcoded count matrix columns and cell metadata. ",
         "Did 02_barcode_rc_map.R succeed?")
  }
  counts <- counts[, common_bc, drop = FALSE]
  meta   <- meta[match(common_bc, meta$barcode)]

  novel_ids <- identify_novel_srrm3_isoforms(rownames(counts))
  if (length(novel_ids) == 0) {
    message("No novel Srrm3 isoform candidate found. Writing empty result and exiting.")
    writeLines("", file.path(TABLES, "novel_srrm3_per_cluster.tsv"))
    return(invisible(NULL))
  }
  novel_ids <- intersect(novel_ids, rownames(counts))

  # TODO: list of canonical annotated Srrm3 isoform IDs — fill in from Gencode vM25.
  canonical_ids <- character(0)
  # Pattern matching as a fallback: IsoQuant carries gene_name in the transcript
  # table; if we have a gene->transcript mapping cache, use it instead.

  # 1) Per-cluster read counts
  per_cluster <- meta[, .(
    n_cells              = .N,
    cells_with_novel     = sum(colSums(counts[novel_ids, , drop = FALSE] > 0) > 0),
    total_novel_reads    = sum(counts[novel_ids, , drop = FALSE]),
    total_canonical_reads = if (length(canonical_ids)) sum(counts[canonical_ids, , drop = FALSE]) else NA_real_
  ), by = cell_type]
  per_cluster[, pct_cells_novel := 100 * cells_with_novel / n_cells]
  per_cluster[, novel_PSI := ifelse(
    is.na(total_canonical_reads) | (total_novel_reads + total_canonical_reads) == 0,
    NA_real_,
    total_novel_reads / (total_novel_reads + total_canonical_reads)
  )]
  fwrite(per_cluster, file.path(TABLES, "novel_srrm3_per_cluster.tsv"), sep = "\t")

  # 2) Cell-type enrichment.
  #
  # The original plan was a 2x2 Fisher test (cells detecting novel x
  # progenitor vs other). But the long-read library has ZERO OPC + iDG
  # cells (the rare progenitor populations exist only in the paired
  # short-read library — see docs/LEARNING_06 §"Progenitor cells in
  # this dataset"). With n_progenitor = 0 the Fisher test is undefined.
  #
  # Reframe: rather than test enrichment in a population that's absent,
  # we test PER-CLUSTER enrichment via Fisher across all clusters present
  # in the long-read data. This identifies which cell types (if any)
  # show preferential expression of the novel isoform.
  detects_novel <- colSums(counts[novel_ids, , drop = FALSE] > 0) > 0
  per_cluster_fisher <- meta[, {
    cluster_mask <- cell_type == .BY$cell_type
    tab <- table(
      cluster = factor(cluster_mask,    levels = c(FALSE, TRUE)),
      novel   = factor(detects_novel,    levels = c(FALSE, TRUE))
    )
    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
      list(odds_ratio = NA_real_, p_value = NA_real_,
           ci_lo = NA_real_, ci_hi = NA_real_,
           n_cluster = sum(cluster_mask), n_novel_in_cluster = tab["TRUE", "TRUE"])
    } else {
      f <- fisher.test(tab)
      list(odds_ratio = as.numeric(f$estimate),
           p_value    = f$p.value,
           ci_lo      = f$conf.int[1], ci_hi = f$conf.int[2],
           n_cluster  = sum(cluster_mask),
           n_novel_in_cluster = tab["TRUE", "TRUE"])
    }
  }, by = cell_type]
  per_cluster_fisher[, p_adj_BH := p.adjust(p_value, method = "BH")]
  setorder(per_cluster_fisher, p_adj_BH, -odds_ratio)
  fwrite(per_cluster_fisher, file.path(TABLES, "novel_srrm3_per_cluster_enrichment.tsv"),
         sep = "\t")
  message("Per-cluster Fisher enrichment written. Top 5 hits by adjusted p:")
  print(head(per_cluster_fisher, 5))

  # 3) Per-cluster PSI plot (if canonical_ids available)
  if (any(!is.na(per_cluster$novel_PSI))) {
    p <- ggplot(per_cluster,
                aes(x = reorder(cell_type, -novel_PSI), y = novel_PSI)) +
      geom_col(fill = "#6A4C93") +
      geom_hline(yintercept = c(0.05, 0.27, 0.57),
                 linetype = "dashed", color = "grey50") +
      annotate("text", x = 1, y = 0.57, label = "Parental (short-read)",
               hjust = 0, vjust = -0.5, size = 3, color = "grey30") +
      annotate("text", x = 1, y = 0.27, label = "Pos (short-read)",
               hjust = 0, vjust = -0.5, size = 3, color = "grey30") +
      annotate("text", x = 1, y = 0.05, label = "KO (short-read)",
               hjust = 0, vjust = -0.5, size = 3, color = "grey30") +
      labs(x = NULL, y = "Novel Srrm3 isoform PSI",
           title = "Novel Srrm3 isoform inclusion across hippocampal clusters (P25)",
           subtitle = "Dashed lines: short-read bulk PSI from 90-1239779069 (WT 57% / Pos 27% / KO 5%)") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(FIGURES, "novel_srrm3_cluster_psi.pdf"),
           p, width = 9, height = 5)
  }

  invisible(list(per_cluster = per_cluster, per_cluster_fisher = per_cluster_fisher))
}

if (sys.nframe() == 0) main()
