#!/bin/bash
set -euo pipefail
CL=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/sqanti3/ribostamp_srrm3_classification.txt_tmp

echo "=== Header columns ==="
head -1 "$CL" | tr '\t' '\n' | nl

echo ""
echo "=== Our 3 high-confidence novel transcripts ==="
head -1 "$CL"
for tid in transcript11.chr5.nnic transcript157.chr5.nnic transcript185.chr5.nnic; do
  grep -P "^${tid}\t" "$CL" || echo "  (not found: $tid)"
done

echo ""
echo "=== Key columns extracted (associated_gene, structural_category, RTS_stage, perc_A_downstream_TTS) ==="
awk -F'\t' '
  NR==1 {
    for(i=1;i<=NF;i++) {
      if($i=="isoform") c_iso=i;
      if($i=="associated_gene") c_gene=i;
      if($i=="structural_category") c_sc=i;
      if($i=="subcategory") c_sub=i;
      if($i=="RTS_stage") c_rts=i;
      if($i=="all_canonical") c_can=i;
      if($i=="perc_A_downstream_TTS") c_aprime=i;
      if($i=="polyA_motif") c_pam=i;
      if($i=="polyA_motif_found") c_pamf=i;
      if($i=="coding") c_cod=i;
      if($i=="ORF_length") c_orf=i;
      if($i=="predicted_NMD") c_nmd=i;
    }
    printf "%s\n", "isoform\tgene\tcategory\tsubcategory\tRTS_stage\tall_canonical\tperc_A_downstream\tpolyA_motif\tpolyA_found\tcoding\tORF_length\tNMD"
    next
  }
  $c_gene ~ /Srrm3|ENSMUSG00000039860/ || $c_iso ~ /transcript(11|157|185)\.chr5\.nnic/ {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      $c_iso,$c_gene,$c_sc,$c_sub,$c_rts,$c_can,$c_aprime,$c_pam,$c_pamf,$c_cod,$c_orf,$c_nmd
  }
' "$CL" | column -t -s $'\t'
