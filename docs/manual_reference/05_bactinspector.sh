#!/usr/bin/env bash
set -euo pipefail

# Step 5 - Species confirmation (BactInspector)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FASTA="/path/to/spades/${SAMPLE}/contigs.fasta"
OUTDIR="/path/to/output/bactinspector/${SAMPLE}"
# =====================

mkdir -p "${OUTDIR}"
bactinspector check_species -f "${FASTA}" -o "${OUTDIR}"

echo "Done. Result in ${OUTDIR}"
