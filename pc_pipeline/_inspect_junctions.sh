#!/bin/bash
set -euo pipefail
J=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sqanti3/ribostamp_srrm3_junctions.txt_tmp

echo "=== junctions header ==="
head -1 "$J" | tr '\t' '\n' | nl

echo ""
echo "=== junctions for our 3 novel transcripts ==="
head -1 "$J"
for tid in transcript11.chr5.nnic transcript157.chr5.nnic transcript185.chr5.nnic; do
  grep -P "^${tid}\t" "$J" || true
done

echo ""
echo "=== RTS_stage breakdown across all 58 transcripts ==="
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++) if($i=="RTS_stage") c=i; print "RTS_stage column:", c; next }
  c { print $c }
' "$J" | sort | uniq -c | sort -rn
