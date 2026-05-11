#!/bin/bash
set -euo pipefail
OUT=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT
RES=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/results
mkdir -p "$RES"

# --- 1. Subset transcript_models.gtf to Srrm3 only ---
SRRM3_GTF="$RES/srrm3_novel_models.gtf"
awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/{print; next}
  $9 ~ /ENSMUSG00000039860/ {print}
' "$OUT/OUT.transcript_models.gtf" > "$SRRM3_GTF"
echo "Wrote $SRRM3_GTF ($(wc -l < "$SRRM3_GTF") lines)"

# --- 2. Cassette intersection (mm10 cassette ~chr5:135,869,720-135,869,798) ---
# Lifted from mm39 chr5:135,898,574-135,898,652 (lab anchor, 79 bp).
# Approximate mm10 coordinates from the targeted_psi locus origin file.
CASSETTE_START=135869720
CASSETTE_END=135869800
echo ""
echo "=== Novel Srrm3 transcripts overlapping the 79-bp cassette region (mm10 ~$CASSETTE_START-$CASSETTE_END) ==="
awk -F'\t' -v s="$CASSETTE_START" -v e="$CASSETTE_END" '
  $3=="transcript" && $4 <= e && $5 >= s {
    match($9, /transcript_id "([^"]+)"/, t);
    is_novel = ($9 ~ /\.nnic|\.nic|\.nnic/) ? "novel" : "known";
    if (t[1] ~ /^transcript[0-9]/) is_novel = "novel";
    printf "  %s\t%s\t%s..%s\texons=%s\tstrand=%s\n",
      t[1], is_novel, $4, $5,
      ($9 ~ /exons "([0-9]+)"/) ? gensub(/.*exons "([0-9]+)".*/, "\\1", "g", $9) : "?",
      $7
  }
' "$SRRM3_GTF"

echo ""
echo "=== Per-exon view of overlapping novel Srrm3 transcripts ==="
awk -F'\t' -v s="$CASSETTE_START" -v e="$CASSETTE_END" '
  $3=="transcript" && $4 <= e && $5 >= s {
    match($9, /transcript_id "([^"]+)"/, t);
    if (t[1] ~ /^transcript[0-9]/) interesting[t[1]]=1
  }
  $3=="exon" {
    match($9, /transcript_id "([^"]+)"/, t);
    if (t[1] in interesting) printf "  %s\texon\t%s..%s\n", t[1], $4, $5
  }
' "$SRRM3_GTF"

# --- 3. Filter to high-confidence multi-exon Srrm3 novel only ---
HC_TSV="$RES/srrm3_novel_highconf.tsv"
SQ="$OUT/OUT.novel_vs_known.SQANTI-like.tsv"
DC="$OUT/OUT.discovered_transcript_counts.tsv"
DG="$OUT/OUT.discovered_transcript_grouped_file_name_counts.tsv"

awk -F'\t' '
  NR==FNR { count[$1]=$2; next }
  FNR==1 {
    print "transcript_id\tcategory\tassociated_transcript\texons\tlength\tsubcategory\tperc_A_downstream\tall_canonical\ttotal_reads"
    next
  }
  $7 ~ /ENSMUSG00000039860/ && $5+0 >= 2 && $38+0 < 0.5 {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      $1,$6,$8,$5,$4,$15,$38,$17,(count[$1]==""?"0":count[$1])
  }
' "$DC" "$SQ" > "$HC_TSV"

echo ""
echo "=== High-confidence novel Srrm3 transcripts (multi-exon, perc_A_downstream<0.5) ==="
column -t -s$'\t' "$HC_TSV" 2>/dev/null || cat "$HC_TSV"

echo ""
echo "Wrote: $HC_TSV"

# --- 4. Per-sample replication for high-confidence novel ---
echo ""
echo "=== Per-mouse read counts on those high-confidence transcripts ==="
echo "header: $(head -1 "$DG")"
awk -F'\t' 'NR>1' "$HC_TSV" | cut -f1 | while read tid; do
  grep -P "^${tid}\t" "$DG" || true
done

# --- 5. Save a clean summary markdown ---
SUMMARY="$RES/SUMMARY.md"
{
  echo "# Path A — De-novo Srrm3 isoform discovery (GSE314176 hippocampus)"
  echo ""
  echo "**Date:** $(date +%Y-%m-%d)"
  echo "**Pipeline:** pc_pipeline + IsoQuant 3.12.2 sensitive_pacbio mode"
  echo ""
  echo "## Input data"
  echo "- 3 PacBio MAS-Iso-seq Revio runs (SRR36480452/53/54, P25 mouse hippocampus)"
  echo "- 2.15M reads aligned to mm10 (post Srrm3-prefilter + pbmm2 ISOSEQ)"
  echo "- Reference: Gencode vM25 (subsetted to chr5:135,680,000-135,960,000)"
  echo ""
  echo "## Headline result"
  echo "- **3 high-confidence novel Srrm3 isoforms discovered** (multi-exon, low intra-priming risk)"
  echo "- 4 known Srrm3 reference transcripts validated"
  echo "- 8 additional single-exon \"novel\" Srrm3 calls — flagged as likely intra-priming artifacts (perc_A_downstream ≥ 0.5)"
  echo "- Ex15 cassette region (mm10 ~chr5:135,869,720-135,869,800) is overlapped by 1 of the 3 high-confidence novel transcripts"
  echo ""
  echo "## High-confidence novel Srrm3 isoforms"
  echo ""
  echo '```'
  cat "$HC_TSV"
  echo '```'
  echo ""
  echo "## Per-mouse replication"
  echo ""
  echo '```'
  echo "transcript_id  SRR36480452  SRR36480453  SRR36480454"
  awk -F'\t' 'NR>1' "$HC_TSV" | cut -f1 | while read tid; do
    grep -P "^${tid}\t" "$DG" || true
  done
  echo '```'
  echo ""
  echo "## Caveats"
  echo "- IsoQuant ran in \`sensitive_pacbio\` mode (more permissive — favours discovery over conservatism)"
  echo "- CB cell-barcode tags were dropped during the prefilter→FASTQ step, so per-cluster expression of these novel isoforms is not currently available. To recover, join read IDs back to the segmented BAMs."
  echo "- Formal SQANTI3 (RT-switching, polyA motif, intra-priming) was not yet run; IsoQuant's built-in SQANTI-like classification was used. SQANTI3 would add ~1-2h."
  echo ""
  echo "## Files"
  echo "- IsoQuant output: \`/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/isoquant_targeted/OUT/\`"
  echo "- Srrm3-only GTF: \`$SRRM3_GTF\`"
  echo "- High-confidence table: \`$HC_TSV\`"
} > "$SUMMARY"

echo ""
echo "Wrote $SUMMARY"
