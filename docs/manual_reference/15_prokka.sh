#!/usr/bin/env bash
set -euo pipefail

# Step 15 - Annotation (Prokka)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FASTA="/path/to/spades/${SAMPLE}/contigs.fasta"
OUTDIR="/path/to/output/prokka"
GENUS_FLAG=""   # set to "--genus Neisseria" to enable a genus hint
# =====================

prokka --outdir "${OUTDIR}/${SAMPLE}_prokka" --prefix "${SAMPLE}" --force ${GENUS_FLAG} "${FASTA}"

echo "Done. Output: ${OUTDIR}/${SAMPLE}_prokka"
