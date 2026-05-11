#!/bin/bash
# Kill any running pipeline, clean stale prefetch state, relaunch.
# Note: pkill returns non-zero when nothing matches, so we avoid set -e during kills.
THIS_DIR="$(dirname "$(readlink -f "$0")")"

echo "=== Step 1: kill anything running ==="
pkill -9 -f _run_surface_a 2>/dev/null || true
pkill -9 -f 02_align_bulk 2>/dev/null || true
pkill -9 -f STAR 2>/dev/null || true
pkill -9 -f prefetch 2>/dev/null || true
pkill -9 -f fasterq-dump 2>/dev/null || true
sleep 2
remaining=$(pgrep -af 'run_surface_a|02_align_bulk|STAR|prefetch|fasterq-dump' || true)
if [[ -n "$remaining" ]]; then
    echo "Still running:"
    echo "$remaining"
fi

echo ""
echo "=== Step 2: clean stale prefetch state ==="
FQ_BASE=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_fastq_temp
for GSMDIR in "$FQ_BASE"/*/; do
    [[ -d "$GSMDIR" ]] || continue
    # Remove any incomplete SRR subdir (lock files etc.)
    for SRRDIR in "$GSMDIR"/SRR*; do
        [[ -d "$SRRDIR" ]] || continue
        echo "  cleaning $SRRDIR"
        rm -rf "$SRRDIR"
    done
done

# Also clean orphan STAR tmpdirs
for d in /tmp/STAR_*; do
    [[ -d "$d" ]] && rm -rf "$d"
done

echo ""
echo "=== Step 3: show remaining FASTQ state ==="
for GSMDIR in "$FQ_BASE"/*/; do
    [[ -d "$GSMDIR" ]] || continue
    echo "  $GSMDIR:"
    ls -la "$GSMDIR" 2>/dev/null | head -10
done

echo ""
echo "=== Step 4: relaunch ==="
bash "$THIS_DIR/_launch_surface_a.sh"
