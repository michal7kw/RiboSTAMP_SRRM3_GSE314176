#!/bin/bash
# Quick sanity check on per_read.tsv outputs
PERREAD_DIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read

for f in "$PERREAD_DIR"/SRR3648045*.per_read.tsv; do
  [ -f "$f" ] || continue
  echo "--- $(basename "$f") ---"
  awk -F'\t' '
    NR == 1 { next }
    {
      n[$3]++
      if ($5 != "") cb++
      tot++
    }
    END {
      printf "  total reads:    %d\n", tot
      for (k in n) printf "  %-16s %d (%.2f%%)\n", k, n[k], 100*n[k]/tot
      printf "  cb_matched:     %d (%.2f%%)\n", cb, 100*cb/tot
    }
  ' "$f"
  echo
done
