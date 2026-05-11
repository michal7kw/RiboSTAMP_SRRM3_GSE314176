#!/bin/bash
# =============================================================================
# 09_per_cell_progenitor_psi.sh — the answer to the PI's progenitor question
#
# Inputs: per-mouse Srrm3-locus BAMs (output of 08_subset_locus.sh) with CB
# tags from STARsolo, plus the authors' obs.tsv with cell-type assignments.
#
# Steps:
#   A — for each mouse, count Inc1/Inc2/Skip junction reads per CB at the
#       cassette locus (uses 03_count_junctions.py --per-cell)
#   B — join per-CB counts with cluster labels from obs.tsv
#   C — aggregate per cluster (focus on OPC + DG-Immature) with Wilson 95% CI
#   D — write per-cluster TSV + verdict markdown
#
# Outputs:
#   ${SR_PER_CELL_DIR}/{GSM}_per_cell.tsv             (raw per-CB counts)
#   ${SR_RESULTS_DIR}/per_cluster_psi_short_read.tsv  (per-cluster PSI)
#   ${SR_RESULTS_DIR}/progenitor_verdict.md           (final verdict)
# =============================================================================

set -eo pipefail

THIS_DIR="$(dirname "$(readlink -f "$0")")"
# shellcheck disable=SC1091
source "${THIS_DIR}/00_config.sh"

LOG="${SR_LOG_DIR}/09_per_cell_progenitor_psi_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

log_step "09_per_cell_progenitor_psi.sh — testing the progenitor hypothesis"

activate_env

# -----------------------------------------------------------------------------
# Step A — per-cell junction counts
# -----------------------------------------------------------------------------
log_step "Step A — per-cell junction counts (per mouse)"

for GSM in "${SR_10X_GSMS[@]}"; do
    BAM="${SR_BAM_DIR_10X}/${GSM}_srrm3.bam"
    OUT="${SR_PER_CELL_DIR}/${GSM}_per_cell.tsv"

    if [[ ! -s "${BAM}" ]]; then
        echo "WARN: subset BAM missing for ${GSM} (${BAM}); skip"
        continue
    fi
    if [[ -s "${OUT}" ]]; then
        echo "  [skip] ${OUT} exists"
        continue
    fi

    OBS_TSV=$(ls "${BASE_DIR}/metadata/${GSM}_shortread_"*"_obs.tsv" 2>/dev/null | head -1)
    if [[ -z "${OBS_TSV}" ]]; then
        echo "ERROR: no obs.tsv for ${GSM}"
        exit 1
    fi
    echo "Processing ${GSM} (whitelist: ${OBS_TSV})"
    python "${THIS_DIR}/03_count_junctions.py" \
        --bam "${BAM}" \
        --output "${OUT}" \
        --sample-name "${GSM}" \
        --per-cell \
        --cell-whitelist "${OBS_TSV}"
done

# -----------------------------------------------------------------------------
# Step B+C+D — per-cluster aggregation + verdict
# -----------------------------------------------------------------------------
log_step "Step B/C/D — per-cluster aggregation + verdict"

PER_CLUSTER_TSV="${SR_RESULTS_DIR}/per_cluster_psi_short_read.tsv"
PROGENITOR_VERDICT="${SR_RESULTS_DIR}/progenitor_verdict.md"

python <<EOF 2>&1 | tee -a "${LOG}"
import os
import re
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest

META_DIR = Path("${BASE_DIR}") / "metadata"
PER_CELL_DIR = Path("${SR_PER_CELL_DIR}")
RESULTS_DIR = Path("${SR_RESULTS_DIR}")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PROGENITOR_TYPES = {"Glia - OPC", "Neuron - DG - Immature"}

BC_RE = re.compile(r"^([ACGT]+)")
def strip_suffix(s: str) -> str:
    m = BC_RE.match(str(s))
    return m.group(1) if m else str(s)

# 1. Load per-cell junction counts (one TSV per GSM)
per_cell_dfs = []
for tsv in sorted(PER_CELL_DIR.glob("GSM*_per_cell.tsv")):
    df = pd.read_csv(tsv, sep="\t")
    if len(df) > 0:
        per_cell_dfs.append(df)
if not per_cell_dfs:
    raise SystemExit("No per-cell TSVs found in " + str(PER_CELL_DIR))
per_cell = pd.concat(per_cell_dfs, ignore_index=True)
per_cell["bare_cb"] = per_cell["cb"].astype(str).map(strip_suffix)
print(f"Loaded {len(per_cell):,} cell-level rows from {len(per_cell_dfs)} mice")

# 2. Build (mouse-GSM, bare_barcode) -> cluster assignment map
cluster_map = {}
for obs in META_DIR.glob("GSM*_shortread_*_obs.tsv"):
    m = re.match(r"(GSM\d+)_shortread", obs.name)
    if not m:
        continue
    gsm = m.group(1)
    df_obs = pd.read_csv(obs, sep="\t")
    if "barcode" not in df_obs.columns:
        continue
    cell_type_col = "Cell Assignments Grouped" if "Cell Assignments Grouped" in df_obs.columns else None
    cell_subtype_col = "Cell Assignment" if "Cell Assignment" in df_obs.columns else None
    for _, row in df_obs.iterrows():
        bc = strip_suffix(row["barcode"])
        cluster_map[(gsm, bc)] = (
            row[cell_type_col] if cell_type_col else "Unknown",
            row[cell_subtype_col] if cell_subtype_col else "Unknown",
        )
print(f"Loaded {len(cluster_map):,} (mouse, barcode) -> (cell_type, subtype) mappings")

# 3. Annotate per_cell with cluster
def lookup(row):
    return cluster_map.get((row["sample"], row["bare_cb"]), (None, None))
ci = per_cell.apply(lookup, axis=1, result_type="expand")
per_cell["cell_type"] = ci[0]
per_cell["cell_subtype"] = ci[1]

n_with_cluster = per_cell["cell_type"].notna().sum()
print(f"Cells with cluster assignment: {n_with_cluster:,} / {len(per_cell):,} "
      f"({100*n_with_cluster/max(1,len(per_cell)):.1f}%)")

# 4. Aggregate per cluster (rMATS-style PSI)
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
    den = inc1 + inc2 + 2 * skip
    psi, lo, hi = wilson(inc1 + inc2, den) if den else (None, None, None)
    agg_rows.append({
        "cell_type": cell_type,
        "n_cells_with_reads": int(len(sub)),
        "n_inc1": inc1,
        "n_inc2": inc2,
        "n_skip": skip,
        "n_junction_reads": inc1 + inc2 + skip,
        "psi": psi,
        "psi_lo95": lo,
        "psi_hi95": hi,
    })
result = pd.DataFrame(agg_rows).sort_values("n_cells_with_reads", ascending=False)
result.to_csv(RESULTS_DIR / "per_cluster_psi_short_read.tsv", sep="\t", index=False)

print()
print("=" * 70)
print("PER-CLUSTER PSI (short-read junction-counting, rMATS-style)")
print("=" * 70)
print(result.to_string(index=False))

# 5. Verdict
print()
print("=" * 70)
print("PROGENITOR HYPOTHESIS VERDICT")
print("=" * 70)

prog_mask = result["cell_type"].isin(PROGENITOR_TYPES)
mat_mask = ~prog_mask
prog_df = result[prog_mask]
mat_df = result[mat_mask]

prog_inc = int(prog_df["n_inc1"].sum() + prog_df["n_inc2"].sum())
prog_skip = int(prog_df["n_skip"].sum())
prog_den = prog_inc + 2 * prog_skip
prog_psi, prog_lo, prog_hi = wilson(prog_inc, prog_den) if prog_den else (None, None, None)

mat_inc = int(mat_df["n_inc1"].sum() + mat_df["n_inc2"].sum())
mat_skip = int(mat_df["n_skip"].sum())
mat_den = mat_inc + 2 * mat_skip
mat_psi, mat_lo, mat_hi = wilson(mat_inc, mat_den) if mat_den else (None, None, None)

print(f"Progenitors pooled (OPC + DG-Immature):  inc={prog_inc}  skip={prog_skip}  PSI={prog_psi}")
print(f"Mature pooled (everything else):         inc={mat_inc}  skip={mat_skip}  PSI={mat_psi}")

# Verdict logic
if prog_den == 0:
    verdict = "INCONCLUSIVE — no informative junction reads from progenitor cells"
elif prog_psi is not None and mat_psi is not None:
    if prog_psi > 0.05 and prog_psi > 5 * (mat_psi if mat_psi else 0.001):
        verdict = "PROGENITOR HYPOTHESIS CONFIRMED"
    elif prog_psi < 0.01:
        verdict = "PROGENITOR HYPOTHESIS REJECTED"
    else:
        verdict = "INTERMEDIATE — progenitor PSI elevated but not dramatically"
else:
    verdict = "INCONCLUSIVE — insufficient data"

verdict_lines = []
verdict_lines.append("# Surface B (full) — progenitor PSI verdict\n\n")
verdict_lines.append(f"**Date:** {pd.Timestamp.today().date()}\n")
verdict_lines.append(f"**Method:** STARsolo on raw 10x v3 short-read FASTQ (3 mice) → "
                     "per-CB junction counting at the Srrm3 cassette → cluster aggregation\n\n")
verdict_lines.append(f"**Cassette (mm10):** chr5:135,869,720-135,869,798 (- strand, 78 bp)\n\n")

verdict_lines.append("## Per-cluster PSI\n\n")
verdict_lines.append(result.to_markdown(index=False) + "\n\n")

verdict_lines.append("## Pooled progenitor vs mature\n\n")
verdict_lines.append("| Pool | n_inc | n_skip | PSI | 95% CI |\n|---|---:|---:|---:|---|\n")
def fmt_pct(p, lo, hi):
    if p is None: return "—"
    return f"{p*100:.3f}% ({(lo or 0)*100:.3f}%, {(hi or 0)*100:.3f}%)"
verdict_lines.append(f"| Progenitors (OPC + DG-Immature) | {prog_inc} | {prog_skip} | "
                     f"{fmt_pct(prog_psi, prog_lo, prog_hi)} |  |\n")
verdict_lines.append(f"| Mature (all other types) | {mat_inc} | {mat_skip} | "
                     f"{fmt_pct(mat_psi, mat_lo, mat_hi)} |  |\n\n")

verdict_lines.append(f"## Verdict\n\n**{verdict}**\n\n")

verdict_lines.append("## Comparison to other contexts\n\n")
verdict_lines.append("| Context | PSI | Source |\n|---|---:|---|\n")
verdict_lines.append("| Proliferating N2A — Parental (lab anchor, rMATS) | 57% | "
                     "\`90-1239779069/SRRM3_novel_exon/docs/01_ANALYSIS_REPORT.md\` |\n")
verdict_lines.append("| Proliferating N2A — Parental (our calc) | 50.6% | "
                     "\`data/lab_bam_validation/lab_n2a_psi_per_group.tsv\` |\n")
verdict_lines.append("| Proliferating N2A — KO | 5% / 2.7% | same |\n")
verdict_lines.append("| Differentiated cell-line neurons (NT/B15/B60) | 0% | "
                     "\`data/sr_junction_psi/results/cell_line_validation_verdict_final.md\` |\n")
verdict_lines.append("| Adult hippocampus mature types (long-read, 7) | <0.05% | "
                     "\`data/targeted_psi/results/per_cluster_psi.tsv\` |\n")
verdict_lines.append(f"| **OPC + DG-Immature progenitors (THIS analysis)** | "
                     f"**{fmt_pct(prog_psi, prog_lo, prog_hi)}** | this verdict |\n")

with open(RESULTS_DIR / "progenitor_verdict.md", "w") as f:
    f.writelines(verdict_lines)

print()
print(f"Wrote {RESULTS_DIR / 'per_cluster_psi_short_read.tsv'}")
print(f"Wrote {RESULTS_DIR / 'progenitor_verdict.md'}")
print()
print(f"FINAL VERDICT: {verdict}")
EOF

log_step "09_per_cell_progenitor_psi.sh complete"
echo ""
echo "Open ${PROGENITOR_VERDICT} for the final answer."
