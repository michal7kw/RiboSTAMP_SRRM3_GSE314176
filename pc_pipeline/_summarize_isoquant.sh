#!/bin/bash
set -euo pipefail
OUT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT

echo "=== Per-gene novel-transcript counts (transcript_models.gtf) ==="
awk -F'\t' '$3=="transcript"{
  match($9, /gene_id "([^"]+)"/, g);
  print g[1]
}' "$OUT/OUT.transcript_models.gtf" | sort | uniq -c | sort -rn

echo ""
echo "=== Novel vs known per gene ==="
awk -F'\t' '$3=="transcript"{
  match($9, /gene_id "([^"]+)"/, g);
  match($9, /transcript_id "([^"]+)"/, t);
  is_novel = ($9 ~ /reference_transcript/ ? 0 : 1);
  if (t[1] ~ /^transcript[0-9]+\.chr5/) is_novel = 1;
  if (g[1] != "") print g[1] "\t" t[1] "\t" (is_novel ? "novel" : "known")
}' "$OUT/OUT.transcript_models.gtf" | sort -u | awk -F'\t' '{count[$1"\t"$3]++} END{for (k in count) print k"\t"count[k]}' | sort

echo ""
echo "=== SQANTI-like classification summary ==="
SQ="$OUT/OUT.novel_vs_known.SQANTI-like.tsv"
if [ -f "$SQ" ]; then
  echo "header:"
  head -1 "$SQ"
  echo ""
  echo "structural_category counts:"
  awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="structural_category") col=i; next}
              col {print $col}' "$SQ" | sort | uniq -c | sort -rn
  echo ""
  echo "subcategory counts:"
  awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="subcategory") col=i; next}
              col {print $col}' "$SQ" | sort | uniq -c | sort -rn | head -15
else
  echo "(file not found: $SQ)"
fi

echo ""
echo "=== Novel transcripts in Srrm3 (gene_id starting ENSMUSG00000039860) ==="
awk -F'\t' '$3=="transcript"{
  match($9, /gene_id "([^"]+)"/, g);
  if (g[1] ~ /ENSMUSG00000039860/) print
}' "$OUT/OUT.transcript_models.gtf" | head -60

echo ""
echo "=== Read counts on Srrm3 transcripts ==="
COUNT="$OUT/OUT.transcript_counts.tsv"
if [ -f "$COUNT" ]; then
  echo "header: $(head -1 "$COUNT")"
  awk -F'\t' 'NR>1' "$COUNT" | sort -k2,2 -n -r | head -20
fi
