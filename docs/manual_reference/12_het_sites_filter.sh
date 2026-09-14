#!/usr/bin/env bash
set -euo pipefail

# Step 12 - Heterozygous-site QC: per-sample filtering (run for every sample)
# Point this at either the reference-mapped (step 9) or assembly-mapped
# (step 11) VCFs - run the whole step 12-14 sequence once for each branch.

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
RAW_VCF="/path/to/snippy/${SAMPLE}_snippy/${SAMPLE}.vcf.gz"
OUTDIR="/path/to/output/het_sites"
# =====================

mkdir -p "${OUTDIR}"
bcftools annotate -x ^INFO/DP -i 'INFO/AF>0.1 && INFO/AF<0.9 && INFO/DP>10' \
  "${RAW_VCF}" -o "${OUTDIR}/${SAMPLE}_hets.bcf" -O b

echo "Done. Output: ${OUTDIR}/${SAMPLE}_hets.bcf"
