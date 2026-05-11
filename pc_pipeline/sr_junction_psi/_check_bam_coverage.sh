#!/bin/bash
set -eo pipefail
THIS_DIR="$(dirname "$(readlink -f "$0")")"
source "${THIS_DIR}/00_config.sh"
activate_env

BAM=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/GSM8465528.bam

echo "=== BAM stats: ${BAM} ==="
samtools flagstat "${BAM}" | head -5

echo ""
echo "=== Reads at Srrm3 gene region (chr5:135,860,000-135,910,000 = ~50 kb) ==="
samtools view -c "${BAM}" chr5:135860000-135910000

echo ""
echo "=== Reads at cassette ± 200 bp (chr5:135,869,094-135,873,280 = ~4 kb) ==="
samtools view -c "${BAM}" chr5:135869094-135873280

echo ""
echo "=== Reads spanning N op anywhere in Srrm3 region (= splice junction reads) ==="
samtools view "${BAM}" chr5:135860000-135910000 | awk '$6 ~ /N/' | wc -l

echo ""
echo "=== Total mapped reads in BAM ==="
samtools view -c -F 4 "${BAM}"

echo ""
echo "=== Top 10 most expressed genes in chr5 (any region with high coverage) ==="
samtools view -F 4 "${BAM}" | awk '{print $3}' | sort | uniq -c | sort -nr | head -5
