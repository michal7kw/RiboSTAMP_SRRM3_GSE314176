#!/bin/bash
# =============================================================================
# pc_pipeline/03 — skera split MAS-Iso-seq concatemers into S-reads.
# Each task uses NUM_THREADS; loop is sequential to avoid disk contention.
# Critical: see docs/LEARNING_05_masseq_skera.md for adapter-kit selection.
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

# Adapter FASTA — pbskera does NOT bundle adapter FASTAs in the conda package.
# They are downloaded by setup_references.sh from PacBio's public dataset
# bucket: https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-MAS_adapters/
# PacBio's directory naming is misleading: v1=mas16, v2=mas12, v3=mas8.
# We use mas16 because the Revio reads are ~13.5 kb (= 16 segments × ~845 bp).
ADAPTER_FASTA="${ADAPTER_FASTA:-${BASE_DIR}/reference/adapters/MAS-Seq_Adapter_v1/mas16_primers.fasta}"
if [ ! -f "${ADAPTER_FASTA}" ]; then
  echo "ERROR: MAS-16 adapter FASTA not found at ${ADAPTER_FASTA}." >&2
  echo "Run pc_pipeline/setup_references.sh (which now mirrors all three" >&2
  echo "PacBio MAS adapter sets) — or set ADAPTER_FASTA in 00_config.sh." >&2
  exit 1
fi
echo "Adapter FASTA: ${ADAPTER_FASTA}"

for SRR in ${SRR_LIST}; do
  IN_BAM="${RAW_DIR}/${SRR}.unaligned.bam"
  OUT_BAM="${SPLIT_DIR}/${SRR}.segmented.bam"

  if [ ! -f "${IN_BAM}" ]; then
    echo "ERROR: missing ${IN_BAM}; run 02_fetch_pacbio.sh first." >&2
    exit 1
  fi
  if [ -s "${OUT_BAM}" ]; then
    log_step "${SRR}: already split; skipping."
    continue
  fi

  log_step "${SRR}: skera split (using ${NUM_THREADS} threads)"
  skera split \
      --num-threads "${NUM_THREADS}" \
      "${IN_BAM}" \
      "${ADAPTER_FASTA}" \
      "${OUT_BAM}"

  SUMMARY="${OUT_BAM%.bam}.summary.json"
  if [ -f "${SUMMARY}" ]; then
    echo "--- ${SRR} skera summary ---"
    python -m json.tool "${SUMMARY}" | head -30
  fi
done

log_step "All splits done"
ls -lh "${SPLIT_DIR}"/*.segmented.bam
