#!/bin/bash
set -euo pipefail
OUT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT

echo "=== Discovered (novel) transcripts — top 30 by total reads ==="
DC="$OUT/OUT.discovered_transcript_counts.tsv"
echo "header: $(head -1 "$DC")"
awk -F'\t' 'NR>1' "$DC" | sort -k2,2 -n -r | head -30

echo ""
echo "=== Novel Srrm3 transcripts only (full GTF lines + read counts) ==="
SRR3_TIDS=$(awk -F'\t' '$3=="transcript"{
  match($9, /gene_id "([^"]+)"/, g);
  match($9, /transcript_id "([^"]+)"/, t);
  if (g[1] ~ /ENSMUSG00000039860/ && t[1] ~ /transcript[0-9]/) print t[1]
}' "$OUT/OUT.transcript_models.gtf")

echo "Novel Srrm3 transcript IDs: $SRR3_TIDS"
echo ""
for tid in $SRR3_TIDS; do
  count=$(awk -F'\t' -v t="$tid" '$1==t {print $2}' "$DC")
  echo "  $tid : $count reads"
  awk -F'\t' -v t="$tid" '$3=="transcript" && $9 ~ "transcript_id \""t"\""' "$OUT/OUT.transcript_models.gtf" | head -1
  echo ""
done

echo ""
echo "=== Per-sample read counts for novel Srrm3 transcripts ==="
DG="$OUT/OUT.discovered_transcript_grouped_file_name_counts.tsv"
if [ -f "$DG" ]; then
  echo "header: $(head -1 "$DG")"
  for tid in $SRR3_TIDS; do
    grep -P "^${tid}\t" "$DG" || true
  done
fi
