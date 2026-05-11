#!/bin/bash
# Deep status check.
echo "=== All processes spawned by Surface A tree ==="
ps -ef --forest | grep -E '(_run_surface_a|02_align|prefetch|fasterq|STAR|samtools|gzip)' | grep -v grep || echo "  (none)"

echo ""
echo "=== FASTQ_temp directory contents ==="
FQDIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_fastq_temp
if [[ -d "$FQDIR" ]]; then
    ls -la "$FQDIR" 2>&1
    for d in "$FQDIR"/*/; do
        [[ -d "$d" ]] || continue
        echo "  -- ${d}:"
        ls -la "$d" 2>&1
    done
fi

echo ""
echo "=== STAR _star_logs ==="
SLDIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_star_logs
if [[ -d "$SLDIR" ]]; then
    ls -la "$SLDIR" 2>&1
fi

echo ""
echo "=== Latest log (last 50 lines) ==="
LOG=$(ls -t /mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/logs/02_align_bulk_*.log 2>/dev/null | head -1)
if [[ -n "$LOG" ]]; then
    echo "Log: $LOG"
    tail -50 "$LOG"
fi
