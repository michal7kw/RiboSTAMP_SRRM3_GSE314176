#!/bin/bash
# Probe how many PacBio-format constraints skera enforces when given an
# SRA-derived BAM + synthesized aux fields. If this works at small scale,
# it's a free path; if it keeps requiring more synth, the AWS pay path is
# the right answer.
set -e
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate longread_iso

SRA=/tmp/SRR36480452/SRR36480452.sra
TEST_BAM=/tmp/test_synthtag.bam
TEST_OUT=/tmp/test_synthtag.split.bam
ADAPTER=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/adapters/MAS-Seq_Adapter_v1/mas16_primers.fasta

RG_ID="m_synthetic_pacbio_RG"

echo "=== build BAM with PacBio-style QNAME + @RG + synth aux tags ==="
{
  echo -e "@HD\tVN:1.6\tSO:unknown\tpb:5.0.0"
  printf "@RG\tID:%s\tPL:PACBIO\tPM:REVIO\tPU:%s\tSM:1\tDS:READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100\n" "$RG_ID" "$RG_ID"
  sam-dump --unaligned "$SRA" 2>/dev/null | head -10000 \
    | awk -v rg="$RG_ID" 'BEGIN{OFS="\t"} {
        if ($0 ~ /^@/) next
        seq_len = length($10)
        # Replace QNAME with PacBio-style m<movie>/<zmw>/ccs
        $1 = rg "/" NR "/ccs"
        # Strip the original RG:Z: tag if present, build clean record
        out = ""
        for (i = 1; i <= 11; i++) out = out (i == 1 ? "" : "\t") $i
        # Append our tags + reset RG to match @RG header
        printf "%s\tRG:Z:%s\trq:f:0.99\tnp:i:15\tqs:i:0\tqe:i:%d\n", out, rg, seq_len
      }'
} | samtools view -b -h - > "$TEST_BAM"

echo "BAM size: $(stat -c%s $TEST_BAM)"
echo ""
echo "=== verify header + first record format ==="
samtools view -H "$TEST_BAM" | head -3
echo "--- first record (truncated) ---"
samtools view "$TEST_BAM" | head -1 | cut -c-200
echo "--- aux tags of first record ---"
samtools view "$TEST_BAM" | head -1 | tr '\t' '\n' | tail -7

echo ""
echo "=== run skera split ==="
rm -f "$TEST_OUT" "$TEST_OUT".*
skera split --num-threads 4 --log-level INFO "$TEST_BAM" "$ADAPTER" "$TEST_OUT" 2>&1 | head -30

echo ""
echo "=== output ==="
ls -lh "$TEST_OUT"* 2>&1 | head
SUMMARY="${TEST_OUT%.bam}.summary.json"
[ -f "$SUMMARY" ] && cat "$SUMMARY" | python -m json.tool | head -25
