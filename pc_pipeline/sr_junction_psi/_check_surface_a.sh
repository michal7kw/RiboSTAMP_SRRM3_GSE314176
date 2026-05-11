#!/bin/bash
# Quick status check for the Surface A background pipeline.
echo "=== Active processes ==="
ps -ef | grep -E 'STAR|prefetch|fasterq-dump|surface_a|02_align_bulk|04_validate' | grep -v grep || echo "  (none)"

echo ""
echo "=== Latest log lines ==="
tail -30 /tmp/surface_a.log 2>/dev/null || echo "  /tmp/surface_a.log not found"

echo ""
echo "=== Disk usage ==="
df -BG /mnt/e 2>/dev/null | tail -1

echo ""
echo "=== STAR index status ==="
STAR_DIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/STAR_mm10_index
if [[ -f "${STAR_DIR}/Genome" ]]; then
    du -sh "${STAR_DIR}" 2>/dev/null
    echo "  Genome file: $(stat -c '%y' "${STAR_DIR}/Genome" 2>/dev/null)"
else
    echo "  Index not yet built (check if STAR is currently running)"
    if [[ -d "${STAR_DIR}" ]]; then
        echo "  Partial: $(ls "${STAR_DIR}" 2>/dev/null | head -5)"
    fi
fi

echo ""
echo "=== Cell-line BAM status ==="
BAMDIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk
if [[ -d "${BAMDIR}" ]]; then
    ls -la "${BAMDIR}"/*.bam 2>/dev/null | head -10
    echo "  Total BAMs: $(ls "${BAMDIR}"/*.bam 2>/dev/null | wc -l) / 9"
else
    echo "  No BAMs yet"
fi
