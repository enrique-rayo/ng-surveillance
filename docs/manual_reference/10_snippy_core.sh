#!/usr/bin/env bash
set -euo pipefail

# Step 10 - Batch-wide core alignment (snippy-core)
# Run this once, after every sample has been through step 9.

# ==== EDIT THESE ====
BATCH="BATCH_NAME"
SNIPPY_DIR="/path/to/output/snippy"
REF_FASTA="/path/to/snippy_reference/FA1090.fasta"
MASK_BED="/path/to/snippy_reference/mask.bed"
# =====================

cd "${SNIPPY_DIR}"
snippy-core --mask "${MASK_BED}" --ref "${REF_FASTA}" \
  --prefix "${BATCH}_snippy_masked" ./*_snippy

echo "Done. Core alignment: ${SNIPPY_DIR}/${BATCH}_snippy_masked.full.aln"
