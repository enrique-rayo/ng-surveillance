#!/usr/bin/env bash
set -euo pipefail

# Step 9 - Reference-mapped variant calling (Snippy, per sample)
# Run this for every sample before moving to step 10.

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FWD="/path/to/trimmomatic/${SAMPLE}_forward_paired.fq.gz"
REV="/path/to/trimmomatic/${SAMPLE}_reverse_paired.fq.gz"
REF_FASTA="/path/to/snippy_reference/FA1090.fasta"
OUTDIR="/path/to/output/snippy"
THREADS=16
# =====================

snippy --outdir "${OUTDIR}/${SAMPLE}_snippy" --ref "${REF_FASTA}" \
  --R1 "${FWD}" --R2 "${REV}" \
  --prefix "${SAMPLE}" --mincov 5 --minfrac 0.75 --unmapped --cpus "${THREADS}"

VCF="${OUTDIR}/${SAMPLE}_snippy/${SAMPLE}.vcf"
bgzip -f "${VCF}"
bcftools index -f "${VCF}.gz"

echo "Done. Output: ${OUTDIR}/${SAMPLE}_snippy"
