#!/bin/bash
set -euo pipefail
F=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT/OUT.read_assignments.tsv.gz

echo "=== Header columns ==="
zcat "$F" | grep -E "^#|^read_id" | head -3 | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) printf "  %d: %s\n", i, $i}'

echo ""
echo "=== Sample rows (first 3 non-header) ==="
zcat "$F" | grep -v "^#" | head -4 | tail -3 | awk -F'\t' '{
  for(i=1;i<=NF;i++) printf "  %d: %s\n", i, $i;
  print "---"
}'

echo ""
echo "=== assignment_type frequency ==="
zcat "$F" | grep -v "^#" | awk -F'\t' 'NR>1 {print $6}' | sort | uniq -c | sort -rn

echo ""
echo "=== Are there inconsistent reads with novel-model hints in assignment_events or additional_info? ==="
zcat "$F" | grep -v "^#" | awk -F'\t' 'NR>1 && $6 ~ /inconsistent/' | head -5 | awk -F'\t' '{print "events:", $7; print "addinfo:", $9; print "---"}'
