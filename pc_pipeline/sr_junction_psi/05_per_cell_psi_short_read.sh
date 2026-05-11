#!/bin/bash
# =============================================================================
# 05_per_cell_psi_short_read.sh — the progenitor question, finally
#
# Runs the junction-PSI calculator on the short-read 10x BAMs (3 mice) at
# per-cell resolution, then aggregates per cell-type cluster INCLUDING the
# OPC and DG-Immature populations that are absent from the long-read library.
#
# This is the analysis that DIRECTLY answers the PI's progenitor question:
# "Is the cassette enriched in progenitor-like cells?"
#
# Outputs:
#   ${SR_RESULTS_DIR}/per_cluster_psi_short_read.tsv  — per-cluster PSI
#   ${SR_RESULTS_DIR}/per_cluster_psi_comparison.png  — bar plot:
#     short-read OPC + iDG vs short-read mature + long-read mature + lab anchor
#   ${SR_RESULTS_DIR}/progenitor_verdict.md           — written interpretation
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/05_per_cell_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "05_per_cell_psi_short_read.sh — testing the progenitor hypothesis"

activate_env

# Pre-flight: 10x BAMs available?
N_BAMS=$(ls "${SR_BAM_DIR_10X}"/*.bam 2>/dev/null | wc -l)
if [[ "${N_BAMS}" -eq 0 ]]; then
    echo "ERROR: no 10x BAMs in ${SR_BAM_DIR_10X}"
    echo "Run 01_fetch_data.sh first"
    exit 1
fi
echo "Found ${N_BAMS} 10x BAMs"

META_DIR="${BASE_DIR}/metadata"

# -----------------------------------------------------------------------------
# Step A — per-cell junction counts (one TSV per mouse)
# -----------------------------------------------------------------------------
log_step "Step A — per-cell junction counts"

for BAM in "${SR_BAM_DIR_10X}"/*.bam; do
    GSM=$(basename "${BAM}" .bam)
    OUT="${SR_PER_CELL_DIR}/${GSM}_per_cell.tsv"

    if [[ -s "${OUT}" ]]; then
        echo "  [skip] ${OUT} already exists"
        continue
    fi

    # Find the matching short-read obs.tsv to use as a whitelist
    OBS_TSV=$(ls "${META_DIR}"/${GSM}_shortread_*_obs.tsv 2>/dev/null | head -1)
    if [[ -z "${OBS_TSV}" ]]; then
        echo "  WARN: no obs.tsv for ${GSM} in ${META_DIR}; will count all CBs"
        WHITELIST_ARG=""
    else
        echo "  using whitelist: ${OBS_TSV}"
        WHITELIST_ARG="--cell-whitelist ${OBS_TSV}"
    fi

    echo "Processing ${GSM}..."
    # shellcheck disable=SC2086
    python "${THIS_DIR}/03_count_junctions.py" \
        --bam "${BAM}" \
        --output "${OUT}" \
        --sample-name "${GSM}" \
        --per-cell ${WHITELIST_ARG}
done

# -----------------------------------------------------------------------------
# Step B — aggregate per-cluster PSI (joining with cell-type assignments)
# -----------------------------------------------------------------------------
log_step "Step B — per-cluster aggregation"

PER_CLUSTER_TSV="${SR_RESULTS_DIR}/per_cluster_psi_short_read.tsv"
PROGENITOR_VERDICT="${SR_RESULTS_DIR}/progenitor_verdict.md"

python <<'EOF' 2>&1 | tee -a "${LOG}"
import os
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binomtest

META_DIR = Path(os.environ["BASE_DIR"]) / "metadata"
PER_CELL_DIR = Path(os.environ["SR_PER_CELL_DIR"])
RESULTS_DIR = Path(os.environ["SR_RESULTS_DIR"])
MIN_READS_PER_CLUSTER = int(os.environ.get("SR_MIN_READS_PER_CLUSTER", 50))

# 1. Load all per-cell junction counts
all_per_cell = []
for tsv in sorted(PER_CELL_DIR.glob("*_per_cell.tsv")):
    df = pd.read_csv(tsv, sep="\t")
    if len(df) > 0:
        all_per_cell.append(df)
if not all_per_cell:
    raise SystemExit("No per-cell TSVs found in " + str(PER_CELL_DIR))
per_cell = pd.concat(all_per_cell, ignore_index=True)
print(f"Loaded {len(per_cell):,} cell-level rows from {len(all_per_cell)} mice")

# Strip 10x -1 suffix if present, normalize CB
per_cell["bare_cb"] = per_cell["cb"].astype(str).str.replace("-1$", "", regex=True).str.split("_").str[0]

# 2. Load short-read obs.tsv for each mouse and build (mouse, cb) -> cluster map
cluster_map = {}
for obs in META_DIR.glob("*_shortread_*_obs.tsv"):
    m = re.match(r"(GSM\d+)_shortread", obs.name)
    if not m:
        continue
    gsm = m.group(1)
    df_obs = pd.read_csv(obs, sep="\t")
    if "barcode" not in df_obs.columns:
        print(f"  WARN: no 'barcode' column in {obs.name}")
        continue
    cell_type_col = "Cell Assignments Grouped" if "Cell Assignments Grouped" in df_obs.columns else None
    cell_subtype_col = "Cell Assignment" if "Cell Assignment" in df_obs.columns else None
    for _, row in df_obs.iterrows():
        bc = str(row["barcode"]).split("_")[0]
        cluster_map[(gsm, bc)] = (
            row[cell_type_col] if cell_type_col else "Unknown",
            row[cell_subtype_col] if cell_subtype_col else "Unknown",
        )
print(f"Loaded {len(cluster_map):,} (mouse, barcode) → (cell_type, subtype) mappings")

# 3. Annotate per_cell with cluster
def lookup(row):
    key = (row["sample"], row["bare_cb"])
    return cluster_map.get(key, (None, None))

cluster_info = per_cell.apply(lookup, axis=1, result_type="expand")
per_cell["cell_type"] = cluster_info[0]
per_cell["cell_subtype"] = cluster_info[1]

n_with_cluster = per_cell["cell_type"].notna().sum()
print(f"Cells with cluster assignment: {n_with_cluster:,} / {len(per_cell):,} "
      f"({100*n_with_cluster/len(per_cell):.1f}%)")

# 4. Aggregate per cluster
def wilson(num, den):
    if den == 0:
        return None, None, None
    res = binomtest(int(num), int(den))
    lo, hi = res.proportion_ci(method="wilson")
    return num/den, lo, hi

agg_rows = []
for cell_type, sub in per_cell[per_cell["cell_type"].notna()].groupby("cell_type"):
    inc1 = int(sub["n_inc1"].sum())
    inc2 = int(sub["n_inc2"].sum())
    skip = int(sub["n_skip"].sum())
    psi, lo, hi = wilson(inc1 + inc2, inc1 + inc2 + 2 * skip)
    agg_rows.append({
        "cell_type": cell_type,
        "n_cells_with_reads": len(sub),
        "n_inc1": inc1,
        "n_inc2": inc2,
        "n_skip": skip,
        "n_total_junction_reads": inc1 + inc2 + skip,
        "psi": psi,
        "psi_lo95": lo,
        "psi_hi95": hi,
    })
result = pd.DataFrame(agg_rows).sort_values("n_cells_with_reads", ascending=False)
result.to_csv(RESULTS_DIR / "per_cluster_psi_short_read.tsv", sep="\t", index=False)

print()
print("=" * 70)
print("PER-CLUSTER PSI (short-read junction-counting)")
print("=" * 70)
print(result.to_string(index=False))

# 5. Verdict — does OPC or iDG-Immature show higher PSI than mature types?
print()
print("=" * 70)
print("PROGENITOR HYPOTHESIS VERDICT")
print("=" * 70)

progenitor_types = [c for c in result["cell_type"] if any(
    t in str(c).lower() for t in ["opc", "immature", "progenitor"]
)]
mature_types = [c for c in result["cell_type"] if c not in progenitor_types]

verdict_lines = []
verdict_lines.append("# Short-read junction-PSI: progenitor verdict\n")
verdict_lines.append(f"**Date:** {pd.Timestamp.today().date()}\n")
verdict_lines.append(f"**Cassette:** chr5:135,869,720–135,869,798 mm10 (− strand, 78 bp)\n\n")
verdict_lines.append("## Per-cluster PSI\n\n")
verdict_lines.append(result.to_markdown(index=False) + "\n\n")

if not progenitor_types:
    verdict = "INCONCLUSIVE: no progenitor-like cell types found in cluster labels"
    verdict_lines.append(f"## Verdict: {verdict}\n")
elif result.loc[result["cell_type"].isin(progenitor_types), "n_total_junction_reads"].sum() < MIN_READS_PER_CLUSTER:
    verdict = "INCONCLUSIVE: progenitor populations have <{} junction reads — too sparse to call".format(MIN_READS_PER_CLUSTER)
    verdict_lines.append(f"## Verdict: {verdict}\n")
else:
    prog_psi = result.loc[result["cell_type"].isin(progenitor_types), "psi"].dropna()
    mat_psi = result.loc[result["cell_type"].isin(mature_types), "psi"].dropna()
    prog_max = prog_psi.max() if len(prog_psi) else None
    mat_max = mat_psi.max() if len(mat_psi) else None
    verdict_lines.append(f"## Verdict\n\n")
    verdict_lines.append(f"- Progenitor types: {progenitor_types}\n")
    verdict_lines.append(f"- Max PSI in progenitor types: {prog_max}\n")
    verdict_lines.append(f"- Max PSI in mature types: {mat_max}\n\n")
    if prog_max is not None and prog_max > 0.05 and (mat_max is None or prog_max > 5 * mat_max):
        verdict = "PROGENITOR HYPOTHESIS CONFIRMED — cassette PSI is substantially higher in OPC/iDG than in mature types"
    elif prog_max is not None and prog_max < 0.01:
        verdict = "PROGENITOR HYPOTHESIS REJECTED — cassette PSI is ~0% in progenitors too"
    else:
        verdict = "INTERMEDIATE — progenitor PSI is non-zero but not dramatically enriched"
    verdict_lines.append(f"**{verdict}**\n")

print("\n".join(verdict_lines))

with open(RESULTS_DIR / "progenitor_verdict.md", "w") as f:
    f.writelines(verdict_lines)

EOF

log_step "05_per_cell_psi_short_read.sh complete"
echo ""
echo "Outputs:"
echo "  Per-cluster PSI:   ${PER_CLUSTER_TSV}"
echo "  Verdict:           ${PROGENITOR_VERDICT}"
echo ""
echo "Open the verdict file to read the answer to the PI's progenitor question."
