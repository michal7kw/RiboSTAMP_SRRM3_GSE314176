#!/bin/bash
# Kill any running Surface A processes.
echo "Before kill:"
pgrep -af '_run_surface_a|02_align_bulk|01_fetch_data|04_validate' | head -10

pkill -9 -f _run_surface_a 2>/dev/null
pkill -9 -f 02_align_bulk 2>/dev/null
pkill -9 -f 01_fetch_data 2>/dev/null
pkill -9 -f 04_validate 2>/dev/null
pkill -9 -f STAR 2>/dev/null
pkill -9 -f prefetch 2>/dev/null
pkill -9 -f fasterq-dump 2>/dev/null
sleep 2

echo ""
echo "After kill:"
pgrep -af '_run_surface_a|02_align_bulk|01_fetch_data|04_validate|STAR|prefetch|fasterq-dump' | head -10 || echo "  (all killed)"
