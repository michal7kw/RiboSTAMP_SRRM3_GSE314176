#!/bin/bash
TSV=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read/all_samples.per_read.tsv

echo "=== column headers ==="
head -1 "$TSV"
echo
echo "=== classification breakdown (n_total = $(($(wc -l < "$TSV") - 1))) ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="class"||$i=="classification") cli=i; next}{n[$cli]++} END{for(k in n) printf "%-15s %10d\n", k, n[k]}' "$TSV"

echo
echo "=== CB match rate ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="cb_matched") cbi=i; next}{
    n_total++;
    if($cbi != "" && $cbi != "NA") n_match++;
}
END{
    printf "matched: %d / %d = %.2f%%\n", n_match, n_total, 100*n_match/n_total
}' "$TSV"

echo
echo "=== CB match rate among informative reads only ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="cb_matched") cbi=i; if($i=="class"||$i=="classification") cli=i}; next}
$cli != "UNINFORMATIVE" {
    n_total++;
    if($cbi != "" && $cbi != "NA") n_match++;
}
END{
    printf "matched: %d / %d = %.2f%%\n", n_match, n_total, 100*n_match/n_total
}' "$TSV"

echo
echo "=== orientation breakdown ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="cb_orientation") oi=i; next}{n[$oi]++} END{for(k in n) printf "%-15s %10d\n", (k==""?"(empty)":k), n[k]}' "$TSV"
