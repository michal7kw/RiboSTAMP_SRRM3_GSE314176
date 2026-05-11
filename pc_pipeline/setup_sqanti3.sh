#!/bin/bash
# =============================================================================
# pc_pipeline/setup_sqanti3.sh — install SQANTI3 in its own conda env.
#
# Why separate from longread_iso: SQANTI3 is not on bioconda. Upstream ships
# its own conda environment YAML and Python module layout that conflicts with
# clean conda dependency resolution. Easiest path is to follow upstream's
# recommended install: clone, create the SQANTI3 env from their YAML, and
# call sqanti3_qc.py from that env.
#
# Run once. After this, hpc/07_sqanti3_qc.sh and pc_pipeline/07_sqanti3_qc.sh
# should activate the SQANTI3.env (not longread_iso) before invoking the tool.
# =============================================================================

set -euo pipefail

# Source conda init explicitly — this script may be invoked from a non-
# interactive shell (e.g. nohup'd background, or `bash setup_sqanti3.sh`
# from a `bash -lc` parent) where ~/.bashrc isn't sourced and conda/mamba
# aren't on PATH otherwise.
CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found at $CONDA_BASE — set CONDA_BASE env var." >&2
  exit 1
fi

SQANTI3_HOME="${SQANTI3_HOME:-$HOME/SQANTI3}"
SQANTI3_VERSION="${SQANTI3_VERSION:-v5.4}"   # adjust if upstream releases newer

echo "=== Step 1: clone SQANTI3 ==="
if [ ! -d "${SQANTI3_HOME}" ]; then
  git clone --branch "${SQANTI3_VERSION}" --depth 1 \
      https://github.com/ConesaLab/SQANTI3.git "${SQANTI3_HOME}"
else
  echo "Already cloned at ${SQANTI3_HOME}"
fi

cd "${SQANTI3_HOME}"

echo ""
echo "=== Step 2: create sqanti3 env from upstream YAML ==="
# Upstream's SQANTI3.conda_env.yml declares 'name: sqanti3' (lowercase).
if conda env list | awk '{print $1}' | grep -qx 'sqanti3'; then
  echo "Env sqanti3 already exists; skipping."
else
  mamba env create -f SQANTI3.conda_env.yml
fi

echo ""
echo "=== Step 3: smoke test ==="
conda activate sqanti3
python "${SQANTI3_HOME}/sqanti3_qc.py" --version || true

echo ""
echo "Done. To use SQANTI3:"
echo "  conda activate sqanti3"
echo "  python ${SQANTI3_HOME}/sqanti3_qc.py ..."
echo ""
echo "The pc_pipeline/07_sqanti3_qc.sh script picks this up via SQANTI3_ENV + SQANTI3_HOME in 00_config.sh."
