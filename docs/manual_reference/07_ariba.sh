#!/usr/bin/env bash
set -euo pipefail

# Step 7 - AMR gene calling (ARIBA)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
R1="/path/to/trimmomatic/${SAMPLE}_forward_paired.fq.gz"
R2="/path/to/trimmomatic/${SAMPLE}_reverse_paired.fq.gz"
ARIBA_DB="/path/to/ariba_db"
OUTDIR="/path/to/output/ariba/${SAMPLE}"
# =====================

ariba run --verbose "${ARIBA_DB}" "${R1}" "${R2}" "${OUTDIR}"

echo "Done. Output: ${OUTDIR}"
