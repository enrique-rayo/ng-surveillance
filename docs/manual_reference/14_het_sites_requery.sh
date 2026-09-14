#!/usr/bin/env bash
set -euo pipefail

# Step 14 - Heterozygous-site QC: re-query each sample against the shared panel
# Run for every sample, after step 13.
#
# This produces the raw site-level table. Turning it into the final
# "Num. het. sites / Num. bins" summary needs an additional binning script
# (het_per_100-mer.py) that isn't included here.

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
RAW_VCF="/path/to/snippy/${SAMPLE}_snippy/${SAMPLE}.vcf.gz"
OUTDIR="/path/to/output/het_sites"
# =====================

bcftools query -R "${OUTDIR}/unique_sites.txt" \
  -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' "${RAW_VCF}" \
  -o "${OUTDIR}/${SAMPLE}_all_het_sites.txt"

echo "Done. Output: ${OUTDIR}/${SAMPLE}_all_het_sites.txt"
