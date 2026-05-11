#!/bin/bash
# Count cells per cluster in the long-read dataset:
#   1. Cells in the published GSE314176 long-read whitelist (per mouse + pooled)
#   2. Cells we actually got reads from at Srrm3 (after CB extraction + matching)

OBS_DIR=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata
PER_READ=/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/per_read/all_samples.per_read.tsv

echo "============================================================"
echo "CELLS IN THE PUBLISHED LONG-READ WHITELIST (per mouse, by cluster)"
echo "============================================================"
for f in "${OBS_DIR}"/GSM*_obs.tsv; do
    mouse=$(basename "$f" | sed -E 's/_obs.tsv$//')
    n=$(awk -F'\t' 'NR>1' "$f" | wc -l)
    echo
    echo ">>> ${mouse}  (n=${n} cells total)"
    awk -F'\t' '
        NR==1 { for(i=1;i<=NF;i++) if($i=="cell_type"||$i=="celltype"||$i=="cluster") c=i; next }
        { n[$c]++ }
        END { for(k in n) printf "%-25s %6d\n", k, n[k] | "sort -k2 -nr" }
    ' "$f"
done

echo
echo "============================================================"
echo "POOLED across 3 mice — published whitelist"
echo "============================================================"
awk -F'\t' '
    FNR==1 { for(i=1;i<=NF;i++) if($i=="cell_type"||$i=="celltype"||$i=="cluster") c=i; next }
    { n[$c]++; total++ }
    END {
        for(k in n) printf "%-25s %6d\n", k, n[k] | "sort -k2 -nr"
        close("sort -k2 -nr")
        printf "\nTotal cells: %d\n", total
    }
' "${OBS_DIR}"/GSM*_obs.tsv

echo
echo "============================================================"
echo "CELLS DETECTED at Srrm3 (≥1 read with matched CB)"
echo "============================================================"
awk -F'\t' '
    NR==1 { next }
    $5 != "" && $5 != "NA" { key = $1 ":" $5; seen[key]++ }
    END { for(k in seen) print k }
' "${PER_READ}" \
| awk -F':' '{print $1, $2}' OFS='\t' \
| sort -u > /tmp/detected_cells.tsv
n_detected=$(wc -l < /tmp/detected_cells.tsv)
echo "Total cells detected with ≥1 Srrm3 read: ${n_detected}"

# Now intersect with each mouse's obs.tsv to get cluster assignment
SRR_TO_GSM_452="GSM9380801"   # mouse3
SRR_TO_GSM_453="GSM9380800"   # mouse2
SRR_TO_GSM_454="GSM9380799"   # mouse1

awk -F'\t' '
    BEGIN {
        srr["SRR36480452"]="GSM9380801"
        srr["SRR36480453"]="GSM9380800"
        srr["SRR36480454"]="GSM9380799"
    }
    FNR==NR {
        bc = $1; sub(/.*\//, "", FILENAME); fname = FILENAME
        next
    }
' /dev/null

# A simpler join: for each mouse load (barcode -> cluster) map, then look up each detected cell
for f in "${OBS_DIR}"/GSM*_obs.tsv; do
    mouse=$(basename "$f" | sed -E 's/_obs.tsv$//')
    # which SRR maps to this mouse?
    case "$mouse" in
      GSM9380801) srr=SRR36480452 ;;
      GSM9380800) srr=SRR36480453 ;;
      GSM9380799) srr=SRR36480454 ;;
      *) continue ;;
    esac

    awk -F'\t' -v srr="$srr" -v mouse="$mouse" -v obs="$f" '
        BEGIN {
            # Read the obs file: bare_barcode -> cell_type
            cmd = "awk -F\"\t\" \"NR==1 {for(i=1;i<=NF;i++) if(\\$i==\\\"barcode\\\"||\\$i==\\\"obs_name\\\"||\\$i==\\\"cell_id\\\"||\\$i==\\\"_index\\\") b=i; for(i=1;i<=NF;i++) if(\\$i==\\\"cell_type\\\"||\\$i==\\\"celltype\\\"||\\$i==\\\"cluster\\\") c=i; next} {bc=\\$b; sub(/_.*/, \\\"\\\", bc); print bc \\\"\t\\\" \\$c}\" " obs
            while ((cmd | getline line) > 0) {
                split(line, a, "\t")
                cluster[a[1]] = a[2]
            }
            close(cmd)
        }
        # /tmp/detected_cells.tsv: SRR<TAB>barcode
        $1 == srr {
            ct = (cluster[$2] == "" ? "(unmapped)" : cluster[$2])
            n[ct]++
        }
        END {
            print ""
            print ">>> " mouse " (" srr ") — cells detected at Srrm3"
            tot=0
            for(k in n) {
                printf "%-25s %6d\n", k, n[k] | "sort -k2 -nr"
                tot += n[k]
            }
            close("sort -k2 -nr")
            printf "Total: %d\n", tot
        }
    ' /tmp/detected_cells.tsv
done

echo
echo "============================================================"
echo "POOLED across 3 mice — cells detected at Srrm3"
echo "============================================================"
# Build the pooled cluster map: prefix barcode with mouse to disambiguate
awk -F'\t' '
    BEGIN {
        # SRR -> mouse map
        srr["SRR36480452"]="GSM9380801"
        srr["SRR36480453"]="GSM9380800"
        srr["SRR36480454"]="GSM9380799"
    }
    # First, read each obs.tsv into per-mouse maps
' /dev/null

# Simpler: combine all 3 obs.tsvs into one (mouse, bare_barcode) -> cluster map, then look up each detected cell
{
    for f in "${OBS_DIR}"/GSM*_obs.tsv; do
        mouse=$(basename "$f" | sed -E 's/_obs.tsv$//')
        case "$mouse" in
          GSM9380801) srr=SRR36480452 ;;
          GSM9380800) srr=SRR36480453 ;;
          GSM9380799) srr=SRR36480454 ;;
          *) continue ;;
        esac
        awk -F'\t' -v srr="$srr" '
            NR==1 { for(i=1;i<=NF;i++) {
                       if($i=="barcode"||$i=="obs_name"||$i=="cell_id"||$i=="_index") b=i
                       if($i=="cell_type"||$i=="celltype"||$i=="cluster") c=i
                   }
                   next
            }
            { bc=$b; sub(/_.*/,"",bc); print srr "\t" bc "\t" $c }
        ' "$f"
    done
} > /tmp/cluster_map.tsv

awk -F'\t' '
    NR==FNR { ct[$1 "\t" $2] = $3; next }
    { key = $1 "\t" $2; ctype = (ct[key] == "" ? "(unmapped)" : ct[key]); n[ctype]++; total++ }
    END {
        for(k in n) printf "%-25s %6d\n", k, n[k] | "sort -k2 -nr"
        close("sort -k2 -nr")
        printf "\nTotal cells detected at Srrm3: %d\n", total
    }
' /tmp/cluster_map.tsv /tmp/detected_cells.tsv

echo
echo "(See data/targeted_psi/results/per_cluster_psi.tsv for the same numbers as input to PSI aggregation.)"
