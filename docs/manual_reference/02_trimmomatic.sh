#!/usr/bin/env bash
set -euo pipefail

# Step 2 - Trimming (Trimmomatic)

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
R1="/path/to/${SAMPLE}_R1_001.fastq.gz"
R2="/path/to/${SAMPLE}_R2_001.fastq.gz"
OUTDIR="/path/to/output/trimmomatic"
THREADS=16
# =====================

mkdir -p "${OUTDIR}"
trimmomatic PE -threads "${THREADS}" \
  "${R1}" "${R2}" \
  "${OUTDIR}/${SAMPLE}_forward_paired.fq.gz" "${OUTDIR}/${SAMPLE}_forward_unpaired.fq.gz" \
  "${OUTDIR}/${SAMPLE}_reverse_paired.fq.gz" "${OUTDIR}/${SAMPLE}_reverse_unpaired.fq.gz" \
  LEADING:25 TRAILING:25 MINLEN:36 SLIDINGWINDOW:4:20

echo "Done. Output: ${OUTDIR}/${SAMPLE}_forward_paired.fq.gz / _reverse_paired.fq.gz"
