#!/bin/bash
# =============================================================================
# pc_pipeline/setup_references.sh — one-time download of mm10 + Gencode vM25
# + the two liftOver chain files. Re-run is safe (skips files that exist).
# =============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
activate_env

mkdir -p "${REF_DIR}" "${CHAIN_DIR}"

log_step "mm10 reference FASTA (~3 GB compressed)"
if [ ! -f "${MM10_REF}" ]; then
  cd "${REF_DIR}"
  wget -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
  gunzip mm10.fa.gz
  samtools faidx mm10.fa
  cd - >/dev/null
else
  echo "Already present: ${MM10_REF}"
fi

log_step "Gencode vM25 GTF (~40 MB compressed)"
if [ ! -f "${GENCODE_GTF}" ]; then
  cd "${REF_DIR}"
  wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
  gunzip gencode.vM25.annotation.gtf.gz
  cd - >/dev/null
else
  echo "Already present: ${GENCODE_GTF}"
fi

log_step "MAS adapter FASTAs (PacBio MAS-Iso-seq) — required by skera split"
ADAPTER_DIR="${BASE_DIR}/reference/adapters"
mkdir -p "${ADAPTER_DIR}"
if [ ! -f "${ADAPTER_DIR}/MAS-Seq_Adapter_v1/mas16_primers.fasta" ]; then
  cd "${ADAPTER_DIR}"
  wget -q -r -np -nH --cut-dirs=4 -A "*.fasta" \
    "https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-MAS_adapters/"
  cd - >/dev/null
  echo "Mirrored adapter FASTAs:"
  find "${ADAPTER_DIR}" -name "*.fasta" -exec ls -lh {} \;
else
  echo "Already present: MAS-Seq adapters under ${ADAPTER_DIR}"
fi

log_step "liftOver chains"
for url in \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz|${CHAIN_MM10_TO_MM39}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToMm10.over.chain.gz|${CHAIN_MM39_TO_MM10}" ; do
  src="${url%|*}"
  dst="${url#*|}"
  if [ ! -f "${dst}" ]; then
    wget -O "${dst}" "${src}"
  else
    echo "Already present: ${dst}"
  fi
done

log_step "Done"
ls -lh "${REF_DIR}" "${CHAIN_DIR}"
