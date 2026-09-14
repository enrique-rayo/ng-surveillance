#!/usr/bin/env bash
set -euo pipefail

# Step 16 - Recombination filtering + phylogeny (Gubbins + IQ-TREE)
# Uses the batch-wide core alignment produced in step 10.

# ==== EDIT THESE ====
BATCH="BATCH_NAME"
FULL_ALN="/path/to/snippy/${BATCH}_snippy_masked.full.aln"
CORE_ALN="/path/to/snippy/${BATCH}_snippy_masked.aln"
OUTDIR="/path/to/output/phylogeny"
THREADS=16
# =====================

mkdir -p "${OUTDIR}"
FIXED_ALN="${OUTDIR}/${BATCH}_snippy_masked.full.Nfixed.aln"

# Snippy marks masked bases as 'X'; Gubbins/IQ-TREE expect 'N'
sed '/^>/!s/X/N/g' "${FULL_ALN}" > "${FIXED_ALN}"

run_gubbins.py --threads "${THREADS}" \
  --first-tree-builder fasttree --tree-builder iqtree \
  "${FIXED_ALN}"

GFF="${FIXED_ALN%.aln}.recombination_predictions.gff"

mask_gubbins_aln.py --aln "${CORE_ALN}" --gff "${GFF}" \
  --out "${OUTDIR}/${BATCH}_snippy_masked.gubbins.aln"

echo "Done. Final alignment: ${OUTDIR}/${BATCH}_snippy_masked.gubbins.aln"
