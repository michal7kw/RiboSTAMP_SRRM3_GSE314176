#!/bin/bash
# One-off — verify what sam-dump actually outputs from a PacBio SRA.
set -e
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

SRA=/tmp/SRR36480452/SRR36480452.sra
[ -f "$SRA" ] || { echo "missing $SRA"; exit 1; }

echo "=== first record raw (text) ==="
sam-dump --unaligned "$SRA" 2>/dev/null | head -1 > /tmp/first_record.sam
echo "Number of fields:"
awk -F'\t' '{print NF}' /tmp/first_record.sam
echo
echo "Each field with name:"
awk -F'\t' '{
    print "QNAME:  " $1
    print "FLAG:   " $2
    print "RNAME:  " $3
    print "POS:    " $4
    print "MAPQ:   " $5
    print "CIGAR:  " $6
    print "RNEXT:  " $7
    print "PNEXT:  " $8
    print "TLEN:   " $9
    print "SEQ_len:  " length($10)
    print "QUAL_len: " length($11)
    print "QUAL_first40: " substr($11, 1, 40)
    for (i = 12; i <= NF; i++) {
        print "AUX[" i "]: " substr($i, 1, 80)
    }
}' /tmp/first_record.sam
