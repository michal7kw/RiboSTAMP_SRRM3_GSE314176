#!/bin/bash
# Look at a sample read's full sequence to understand polyA + barcode structure
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

BAM=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/aligned_locus/SRR36480452.locus.bam

echo "=== first read - full info ==="
samtools view "$BAM" 2>/dev/null | head -1 | awk -F'\t' '
{
  print "QNAME: " $1
  print "FLAG:  " $2
  print "POS:   " $4
  print "MAPQ:  " $5
  print "CIGAR_len: " length($6)
  print "SEQ_len:   " length($10)
  print "QUAL_len:  " length($11)
  print
  print "SEQ first 100bp:  " substr($10, 1, 100)
  print "SEQ last 100bp:   " substr($10, length($10) - 99)
  print
  # polyA stretches anywhere in SEQ
  seq = $10
  while (match(seq, /A{8,}/)) {
    print "polyA at pos " RSTART " (length " RLENGTH "), context: " substr(seq, RSTART-25 < 1 ? 1 : RSTART-25, 25 + RLENGTH + 25)
    seq = substr(seq, RSTART + RLENGTH)
  }
  # polyT stretches
  seq = $10
  while (match(seq, /T{8,}/)) {
    print "polyT at pos " RSTART " (length " RLENGTH "), context: " substr(seq, RSTART-25 < 1 ? 1 : RSTART-25, 25 + RLENGTH + 25)
    seq = substr(seq, RSTART + RLENGTH)
  }
}'
