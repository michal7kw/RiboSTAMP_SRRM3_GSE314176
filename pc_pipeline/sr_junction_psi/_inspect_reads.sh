#!/bin/bash
set -eo pipefail
THIS_DIR="$(dirname "$(readlink -f "$0")")"
source "${THIS_DIR}/00_config.sh"
activate_env

BAM=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/GSM8465528.bam

echo "=== Read length distribution (first 1000 reads) ==="
samtools view "${BAM}" | head -1000 | awk '{print length($10)}' | sort -n | uniq -c | head -10

echo ""
echo "=== Sample 5 reads — show CIGAR ==="
samtools view "${BAM}" chr5:135860000-135910000 | head -5 | awk '{print $6, length($10), $1}'

echo ""
echo "=== STAR Log.final.out (alignment summary) ==="
LOG_FINAL=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_star_logs/GSM8465528_Log.final.out
if [[ -s "$LOG_FINAL" ]]; then
    cat "$LOG_FINAL"
fi

echo ""
echo "=== Splice junctions detected by STAR (top 10 by chr5 read count) ==="
SJ=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_star_logs/GSM8465528_SJ.out.tab
if [[ -s "$SJ" ]]; then
    awk '$1=="chr5"' "$SJ" | sort -k7 -nr | head -10
fi

echo ""
echo "=== Junctions in Srrm3 region detected by STAR ==="
if [[ -s "$SJ" ]]; then
    awk '$1=="chr5" && $2 >= 135860000 && $3 <= 135910000' "$SJ"
fi
