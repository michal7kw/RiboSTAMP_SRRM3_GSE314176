#!/usr/bin/env python3
"""Show STAR splice junctions in Srrm3 region — read SJ.out.tab cleanly."""
import os
import sys

SJ_PATH = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_star_logs/GSM8465528_SJ.out.tab"

# STAR SJ.out.tab uses 1-based both-inclusive intron coords; our junction
# coords were 0-based (Inc1: 135869294→135869720, Inc2: 135869798→135873080).
# Convert STAR's coords by subtracting 1 from start.
INC1_START = 135869294
INC1_END = 135869720
INC2_START = 135869798
INC2_END = 135873080
SKIP_START = INC1_START
SKIP_END = INC2_END
TOL = 5

print("=" * 80)
print("All splice junctions in Srrm3 region (chr5:135,860,000-135,910,000)")
print("=" * 80)
print(f"{'start':>12} {'end':>12} {'strand':>6} {'motif':>5} {'anno':>4} {'uniq':>5} {'multi':>5} {'overhang':>8}")
print("-" * 80)

junctions = []
with open(SJ_PATH) as f:
    for line in f:
        f_chrom, f_start, f_end, f_strand, f_motif, f_anno, f_uniq, f_multi, f_overhang = line.strip().split("\t")
        if f_chrom != "chr5":
            continue
        start = int(f_start)
        end = int(f_end)
        if start < 135860000 or end > 135910000:
            continue
        # STAR uses 1-based intron coords (start = first intronic base, end = last intronic base)
        # Our junction coords are 0-based: Inc1 5' end is the LAST exonic base on left, Inc1 3' end is the FIRST exonic base on right.
        # So STAR's start ≈ our_left + 1, STAR's end ≈ our_right - 1
        our_left = start - 1
        our_right = end + 1
        match = []
        if abs(our_left - INC1_START) <= TOL and abs(our_right - INC1_END) <= TOL:
            match.append("Inc1")
        if abs(our_left - INC2_START) <= TOL and abs(our_right - INC2_END) <= TOL:
            match.append("Inc2")
        if abs(our_left - SKIP_START) <= TOL and abs(our_right - SKIP_END) <= TOL:
            match.append("Skip")
        match_str = ",".join(match) if match else "(other Srrm3 junction)"
        junctions.append((start, end, f_strand, f_motif, f_anno, int(f_uniq), int(f_multi), f_overhang, match_str))

if not junctions:
    print("(no junctions in Srrm3 region)")
else:
    for j in junctions:
        print(f"{j[0]:>12} {j[1]:>12} {j[2]:>6} {j[3]:>5} {j[4]:>4} {j[5]:>5} {j[6]:>5} {j[7]:>8}  {j[8]}")

print()
print(f"Cassette anchors (mm10):")
print(f"  Inc1: {INC1_START} → {INC1_END}")
print(f"  Inc2: {INC2_START} → {INC2_END}")
print(f"  Skip: {SKIP_START} → {SKIP_END}")
