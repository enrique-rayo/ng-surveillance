#!/usr/bin/env bash
set -euo pipefail

# Step 13 - Heterozygous-site QC: build the shared site panel
# Run once, after every sample has been through step 12.

# ==== EDIT THESE ====
OUTDIR="/path/to/output/het_sites"
# =====================

bcftools query -f '%CHROM\t%POS\n' "${OUTDIR}"/*_hets.bcf | sort -u > "${OUTDIR}/unique_sites.txt"

echo "Done. Panel: ${OUTDIR}/unique_sites.txt"
