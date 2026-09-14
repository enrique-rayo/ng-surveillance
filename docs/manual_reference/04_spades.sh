#!/usr/bin/env bash
set -euo pipefail

# Step 4 - Assembly (SPAdes)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FWD="/path/to/trimmomatic/${SAMPLE}_forward_paired.fq.gz"
REV="/path/to/trimmomatic/${SAMPLE}_reverse_paired.fq.gz"
OUTDIR="/path/to/output/spades/${SAMPLE}"
THREADS=16
MEM_GB=64
# =====================

spades.py -1 "${FWD}" -2 "${REV}" --careful -o "${OUTDIR}" -t "${THREADS}" -m "${MEM_GB}"

echo "Done. Assembly: ${OUTDIR}/contigs.fasta"
