#!/bin/bash
# Quick verification of step 02 output BAMs.
set +e
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

DIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/raw_pacbio

echo "=== samtools quickcheck ==="
for BAM in "$DIR"/SRR*.unaligned.bam; do
  printf "%s: " "$(basename "$BAM")"
  if samtools quickcheck "$BAM" 2>/dev/null; then
    echo "valid BAM"
  else
    echo "CORRUPTED"
  fi
done

echo ""
echo "=== first-record check (QNAME + aux tags) ==="
for BAM in "$DIR"/SRR*.unaligned.bam; do
  echo "--- $(basename "$BAM") ---"
  samtools view "$BAM" 2>/dev/null | head -1 2>/dev/null | awk -F'\t' '{
    print "  QNAME:    " $1
    print "  SEQ_len:  " length($10)
    print "  QUAL_len: " length($11)
    for (i = 12; i <= NF; i++) print "  AUX[" i "]: " $i
  }'
done
