#!/bin/bash
ls -lt /mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/*.bam 2>/dev/null | head -10
echo ""
echo "---active processes---"
ps -ef | grep -E "STAR|prefetch|fasterq|02_align" | grep -v grep | head -10
