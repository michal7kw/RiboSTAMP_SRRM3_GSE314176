#!/bin/bash
SJ=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sr_junction_psi/bams/cell_line_bulk/_star_logs/GSM8465528_SJ.out.tab

echo "=== STAR splice junctions in Srrm3 region (chr5:135,860,000-135,910,000) ==="
echo "Columns: chrom start end strand intron-motif annotated unique-reads multi-reads max-overhang"
awk '$1=="chr5" && $2 >= 135860000 && $3 <= 135910000' "$SJ" | sort -k7 -nr

echo ""
echo "=== Specifically at our cassette junctions ==="
echo "Inc1 (135869294 → 135869720), Inc2 (135869798 → 135873080), Skip (135869294 → 135873080)"
echo ""
echo "Looking for matches within ±5 bp:"
awk '
    $1=="chr5" && (
        ($2 >= 135869289 && $2 <= 135869299 && $3 >= 135869715 && $3 <= 135869725) ||
        ($2 >= 135869793 && $2 <= 135869803 && $3 >= 135873075 && $3 <= 135873085) ||
        ($2 >= 135869289 && $2 <= 135869299 && $3 >= 135873075 && $3 <= 135873085)
    )
' "$SJ"

echo ""
echo "(Note: STAR's SJ.out.tab uses 1-based, both endpoints inclusive intron coords)"
