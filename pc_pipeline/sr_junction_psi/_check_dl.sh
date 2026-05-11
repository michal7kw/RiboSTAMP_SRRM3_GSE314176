#!/bin/bash
# Check download progress for the currently-running prefetch
echo "=== Active prefetch process ==="
ps -ef | grep prefetch | grep -v grep | head -3

echo ""
echo "=== Current time ==="
date '+%Y-%m-%d %H:%M:%S'

echo ""
echo "=== FASTQ_temp dir for currently-processing GSM ==="
for d in /mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_fastq_temp/*/; do
    [[ -d "$d" ]] || continue
    echo "  $d"
    du -sh "${d}" 2>/dev/null
    ls -la "${d}" 2>/dev/null | head -10
done

echo ""
echo "=== Most recent log entries ==="
LOG=$(ls -t /mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/logs/02_align_bulk_*.log 2>/dev/null | head -1)
if [[ -n "$LOG" ]]; then
    echo "Log: $LOG"
    tail -15 "$LOG"
fi

echo ""
echo "=== /tmp/STAR_* tmpdirs (orphans?) ==="
ls -la /tmp/STAR_* 2>/dev/null | head -10 || echo "(none)"
