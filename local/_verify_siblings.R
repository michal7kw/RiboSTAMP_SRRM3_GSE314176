suppressPackageStartupMessages(library(data.table))

final    <- fread("results/tables/interneuron_isoforms_FINAL.tsv")
siblings <- fread("results/tables/interneuron_isoforms_siblings.tsv")

cat("=== Verification ===\n\n")

cat(sprintf("siblings rows: %d   FINAL rows: %d\n", nrow(siblings), nrow(final)))
stopifnot(nrow(siblings) == nrow(final))

expected <- c("paired", "single_isoform_gene",
              "all_isoforms_gaba_enriched", "no_expressed_sibling")
unexpected <- setdiff(unique(siblings$pairing_status), expected)
cat(sprintf("unexpected pairing_status values: %d\n", length(unexpected)))
stopifnot(length(unexpected) == 0)

paired <- siblings[pairing_status == "paired"]
cat(sprintf("paired rows: %d  with NA sibling: %d  sibling-in-274: %d\n",
            nrow(paired),
            sum(is.na(paired$sibling_isoform_id)),
            sum(paired$sibling_isoform_id %in% final$isoform_id)))
stopifnot(!any(is.na(paired$sibling_isoform_id)))
stopifnot(!any(paired$sibling_isoform_id %in% final$isoform_id))
stopifnot(all(paired$sibling_logFC < 1))

gene_a <- sub("-[0-9]+$", "", paired$isoform_id)
gene_b <- sub("-[0-9]+$", "", paired$sibling_isoform_id)
mismatch <- gene_a != gene_b
cat(sprintf("same-gene mismatches by id-prefix: %d\n", sum(mismatch)))
if (any(mismatch)) {
  print(paired[mismatch, .(isoform_id, gene, sibling_isoform_id)])
}
stopifnot(all(!mismatch))

cat(sprintf("paired rows with delta > 0:  %d / %d\n",
            sum(paired$delta_log2_ratio > 0), nrow(paired)))

cat("\n=== Top 10 paired rows by delta_log2_ratio ===\n")
top <- siblings[pairing_status == "paired"][order(-delta_log2_ratio)][1:10]
print(top[, .(isoform_id, gene, sibling_isoform_id,
              GABA_iso_GABA_mean, GABA_iso_other_mean,
              sibling_GABA_mean,  sibling_other_mean, delta_log2_ratio)])

cat("\n=== Bottom 10 paired rows by delta_log2_ratio ===\n")
bot <- siblings[pairing_status == "paired"][order(delta_log2_ratio)][1:10]
print(bot[, .(isoform_id, gene, sibling_isoform_id, delta_log2_ratio)])

cat("\n=== Canonical GABA-marker pairs ===\n")
ctrl <- siblings[gene %in% c("Slc32a1","Gad2","Slc6a1","Gad1")]
print(ctrl[, .(isoform_id, gene, sibling_isoform_id,
               GABA_iso_GABA_mean, GABA_iso_other_mean,
               sibling_GABA_mean, sibling_other_mean,
               delta_log2_ratio, pairing_status)])

cat("\n=== Meg3 isoforms (multi-isoform-of-same-gene case) ===\n")
print(siblings[gene == "Meg3", .(isoform_id, sibling_isoform_id,
                                  delta_log2_ratio, pairing_status)])

cat("\nFINAL.tsv has sibling_isoform_id: ",
    "sibling_isoform_id" %in% names(final), "\n")
cat("FINAL.tsv has delta_log2_ratio:   ",
    "delta_log2_ratio" %in% names(final), "\n")
cat("FINAL.tsv has pairing_status:     ",
    "pairing_status" %in% names(final), "\n")
cat("FINAL.tsv ncol: ", ncol(final), "\n")

cat("\n=== pairing_status counts (siblings.tsv) ===\n")
print(siblings[, .N, by = pairing_status])
