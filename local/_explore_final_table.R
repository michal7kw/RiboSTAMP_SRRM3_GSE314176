suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

PROJ_ROOT <- here::here()
ino <- fread(file.path(PROJ_ROOT, "results/tables/interneuron_isoforms_FINAL.tsv"))

cat(sprintf("Total isoforms: %d\nColumns: %d\n\n", nrow(ino), ncol(ino)))

# 1. The 9 DESeq2-significant isoforms — what genes are they?
cat("=== Tier-1: DESeq2-significant (padj<0.05 AND log2FC>1) — n=9 ===\n")
tier1 <- ino[DE_padj < 0.05 & DE_log2FC > 1]
print(tier1[, .(isoform_id, gene, GABA_pct, logFC,
                DE_log2FC, DE_padj,
                editsC_is_translated,
                biotype, syngo_quick)])

cat("\n=== Tier-2 (translated + protein-coding) — top 20 by GABA_pct ===\n")
tier2 <- ino[biotype == "protein_coding" & editsC_is_translated == TRUE][order(-GABA_pct)]
print(head(tier2[, .(isoform_id, gene, GABA_pct, logFC,
                     editsC_GABA_pct,
                     DE_padj, biotype, syngo_quick)], 20))

cat(sprintf("\nTier-2 size: %d isoforms (protein-coding + translated)\n", nrow(tier2)))

# 2. Subtype-specific top hits (Cck — most numerous in hippocampus)
cat("\n=== Cck-enriched (Cck_pct >= 50, mean >= 1) — top 10 ===\n")
cck_top <- ino[Cck_pct >= 50 & Cck_mean >= 1][order(-Cck_pct)]
print(head(cck_top[, .(isoform_id, gene, Cck_pct, Cck_mean,
                       editsC_is_translated, biotype, syngo_quick)], 10))

# 3. Synaptic + translated — highest-confidence wet-lab targets
cat("\n=== SynGO + translated — full list (n=5) ===\n")
syn <- ino[syngo_quick == TRUE & editsC_is_translated == TRUE]
print(syn[, .(isoform_id, gene, GABA_pct, logFC,
              DE_padj, editsC_GABA_pct, biotype)])

# 4. Retained-intron isoforms (interesting biotype — n=44)
cat("\n=== Retained-intron isoforms (n=44) — top 10 by GABA_pct ===\n")
ri <- ino[biotype == "retained_intron"][order(-GABA_pct)]
print(head(ri[, .(isoform_id, gene, GABA_pct, logFC,
                  editsC_is_translated, DE_padj)], 10))

# 5. Top hits also in DESeq2-significant set (overlap with formal stats)
cat("\n=== DESeq2-significant overlap with each subset ===\n")
sig_genes <- tier1$gene
cat(sprintf("  Tier-1 gene names: %s\n", paste(sig_genes, collapse = ", ")))
cat(sprintf("  Tier-1 ∩ translated: %d\n",
            sum(tier1$editsC_is_translated == TRUE, na.rm = TRUE)))
cat(sprintf("  Tier-1 ∩ SynGO synaptic: %d\n",
            sum(tier1$syngo_quick == TRUE, na.rm = TRUE)))

# 6. Top genes by combined "tier" score
cat("\n=== Top 15 by combined evidence (logFC × GABA_pct, weighted by translation) ===\n")
ino[, combined_score := logFC * (GABA_pct / 100) * ifelse(editsC_is_translated, 2, 1)]
top_combined <- ino[order(-combined_score)]
print(head(top_combined[, .(isoform_id, gene, GABA_pct, logFC,
                            editsC_is_translated, biotype, combined_score)], 15))

# 7. GO BP overrepresentation (informal)
cat("\n=== Top GO BP terms among the 274 isoforms ===\n")
go_terms <- unlist(strsplit(ino$go_bp[!is.na(ino$go_bp)], "; "))
top_terms <- sort(table(go_terms), decreasing = TRUE)
print(head(top_terms, 20))
