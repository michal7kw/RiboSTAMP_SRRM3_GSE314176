#!/bin/bash
# Subset Gencode vM25 to the Srrm3 locus (chr5:135,680,000-135,960,000) so
# IsoQuant's gffutils GTF→db conversion finishes in seconds instead of hours.
set -euo pipefail

GTF=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/mm10/gencode.vM25.annotation.gtf
OUT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/mm10/gencode.vM25.srrm3_locus.gtf

awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {print; next}
  $1=="chr5" && $4 <= 135960000 && $5 >= 135680000 {print}
' "$GTF" > "$OUT"

echo "input lines:  $(wc -l < "$GTF")"
echo "output lines: $(wc -l < "$OUT")"
echo "output size:  $(stat -c%s "$OUT") bytes"
echo "--- genes in subset ---"
awk -F'\t' '$3=="gene"{
  match($9, /gene_name "([^"]+)"/, a); print a[1]
}' "$OUT" | sort -u
