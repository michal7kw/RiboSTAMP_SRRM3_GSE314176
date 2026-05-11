#!/bin/bash
# Internal helper — fetch the MAS-16 adapter FASTA from PacBio's published
# repos. The pbskera conda package does NOT ship the adapter FASTA; users
# must provide it as a positional arg to `skera split`. PacBio publishes
# the canonical sequences in the Kinnex-single-cell-RNA repo.

set -euo pipefail

OUT_DIR="/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/reference/adapters"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

URLS=(
  "https://raw.githubusercontent.com/PacificBiosciences/skera/master/adapters/mas16_primers.fasta"
  "https://raw.githubusercontent.com/PacificBiosciences/skera/master/data/mas16_primers.fasta"
  "https://raw.githubusercontent.com/PacificBiosciences/Kinnex-single-cell-RNA/master/scripts/mas16_primers.fasta"
  "https://raw.githubusercontent.com/PacificBiosciences/Kinnex-single-cell-RNA/master/mas16_primers.fasta"
  "https://raw.githubusercontent.com/PacificBiosciences/skera/main/adapters/mas16_primers.fasta"
  "https://raw.githubusercontent.com/PacificBiosciences/MAS-Seq-Adapter/master/mas16_primers.fasta"
)

for url in "${URLS[@]}"; do
  echo "trying: $url"
  if wget -q "$url" -O mas16_primers.fasta.tmp && [ -s mas16_primers.fasta.tmp ]; then
    # Verify it looks like a FASTA (starts with > and has nucleotide lines)
    if head -1 mas16_primers.fasta.tmp | grep -q '^>'; then
      mv mas16_primers.fasta.tmp mas16_primers.fasta
      echo "  -> success: $(wc -l < mas16_primers.fasta) lines, $(wc -c < mas16_primers.fasta) bytes"
      head -6 mas16_primers.fasta
      exit 0
    fi
  fi
  rm -f mas16_primers.fasta.tmp
done

echo "ERROR: could not fetch mas16_primers.fasta from any known URL." >&2
exit 1
