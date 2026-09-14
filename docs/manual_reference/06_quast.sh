#!/usr/bin/env bash
set -euo pipefail

# Step 6 - Assembly QC (QUAST) - run once per batch, across all assemblies together

# ==== EDIT THESE ====
ASSEMBLIES=(
  "/path/to/spades/SAMPLE1/contigs.fasta"
  "/path/to/spades/SAMPLE2/contigs.fasta"
)
OUTDIR="/path/to/output/QUAST/quast_results_summary"
THREADS=8
# =====================

quast.py "${ASSEMBLIES[@]}" -o "${OUTDIR}" --threads "${THREADS}"

echo "Done. Report: ${OUTDIR}/report.tsv"
