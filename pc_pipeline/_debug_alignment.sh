#!/bin/bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

BAM=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/aligned_locus/SRR36480452.locus.bam

echo "=== first 5 records: POS + CIGAR (truncated) ==="
samtools view "$BAM" 2>/dev/null | head -5 | awk -F'\t' '{
  print "POS="$4, "MAPQ="$5, "CIGAR_len="length($6), "CIGAR_first80="substr($6, 1, 80)
}'

echo
echo "=== POS histogram (10 kb bins, first 100K records) ==="
samtools view "$BAM" 2>/dev/null | head -100000 | awk -F'\t' '{
  bin = int($4 / 10000) * 10000
  n[bin]++
}
END {
  for (k in n) printf "  pos %d - %d : %d reads\n", k, k+9999, n[k]
}' | sort -n -k2

echo
echo "=== reads spanning the cassette region (29950-30130) ==="
# Reads whose POS is < 29950 (should start before cassette+50bp upstream)
# AND whose POS+SEQ_LEN extends past 30128 (cassette+50bp downstream)
samtools view "$BAM" 2>/dev/null | head -200000 | awk -F'\t' '
{
  pos = $4
  # crude span estimate from CIGAR: count consumes-ref ops
  cigar = $6
  span = 0
  num = ""
  for (i = 1; i <= length(cigar); i++) {
    c = substr(cigar, i, 1)
    if (c >= "0" && c <= "9") { num = num c }
    else {
      if (c == "M" || c == "D" || c == "N" || c == "=" || c == "X") span += num
      num = ""
    }
  }
  end = pos + span - 1
  if (pos < 29950 && end > 30128) spanning++
  total++
}
END {
  printf "  total examined : %d\n", total
  printf "  spanning cassette: %d (%.2f%%)\n", spanning, 100*spanning/total
}'

echo
echo "=== samples reads' span around the cassette (first 5 spanning the locus midpoint) ==="
samtools view "$BAM" chr5:135839721-135899798 2>/dev/null | head -1000 | awk -F'\t' '
$4 < 28000 && $4 > 0 {
  pos = $4
  print "POS="pos, "CIGAR_len="length($6), "first40="substr($6,1,40)"...", "last20="substr($6,length($6)-20)
  shown++
  if (shown >= 5) exit
}'
