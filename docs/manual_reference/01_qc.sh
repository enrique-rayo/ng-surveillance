#!/usr/bin/env bash
set -euo pipefail

# Step 1 - Read QC (FastQC + MultiQC)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
R1="/path/to/${SAMPLE}_R1_001.fastq.gz"
R2="/path/to/${SAMPLE}_R2_001.fastq.gz"
OUTDIR="/path/to/output/fastqc"
# =====================

mkdir -p "${OUTDIR}"
fastqc "${R1}" -o "${OUTDIR}"
fastqc "${R2}" -o "${OUTDIR}"

echo "Done. Once every sample is processed, run 'multiqc .' inside ${OUTDIR}."
