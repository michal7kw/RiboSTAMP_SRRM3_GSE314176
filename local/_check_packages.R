cat("R version:", R.version.string, "\n")
cat("Installed packages check:\n")
for (p in c("data.table", "Matrix", "here", "DESeq2", "biomaRt")) {
  cat(sprintf("  %-12s: %s\n", p, requireNamespace(p, quietly = TRUE)))
}
