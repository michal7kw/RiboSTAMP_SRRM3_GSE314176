#!/bin/bash
# Quick diagnostic on pbmm2 alignment progress.
echo "=== I/O bytes read per pbmm2 ==="
for pid in 16020 16021 16024; do
  if [ -d /proc/$pid ]; then
    rbytes=$(awk '/^read_bytes:/ {print $2}' /proc/$pid/io 2>/dev/null)
    rchar=$(awk '/^rchar:/ {print $2}' /proc/$pid/io 2>/dev/null)
    state=$(awk '/^State:/ {print $2}' /proc/$pid/status 2>/dev/null)
    rss=$(awk '/^VmRSS:/ {print $2, $3}' /proc/$pid/status 2>/dev/null)
    echo "  PID $pid: state=$state, RSS=$rss"
    echo "    rchar (logical bytes read):  $(numfmt --to=iec ${rchar:-0})"
    echo "    read_bytes (disk bytes read): $(numfmt --to=iec ${rbytes:-0})"
    # Check what file the input fd is reading
    in_fd=$(ls -la /proc/$pid/fd/ 2>/dev/null | grep "segmented.bam" | head -1)
    echo "    open input: $in_fd"
  fi
done
