#!/bin/bash
set -euo pipefail
OUT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT
SQ="$OUT/OUT.novel_vs_known.SQANTI-like.tsv"

echo "=== Header columns ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) printf "%d  %s\n", i, $i}' "$SQ"

echo ""
echo "=== Srrm3 novel transcripts in SQANTI-like ==="
awk -F'\t' 'NR==1 || $7 ~ /ENSMUSG00000039860/' "$SQ"

echo ""
echo "=== Top intergenic novel transcripts by length (just to see what they look like) ==="
awk -F'\t' 'NR==1 || $6=="intergenic"' "$SQ" | head -10
