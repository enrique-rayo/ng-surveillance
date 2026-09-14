#!/usr/bin/env bash
set -euo pipefail

# Step 8 - Sequence typing (pyngoST: NG-STAR + MLST)

# ==== EDIT THESE ====
BATCH="BATCH_NAME"
ASSEMBLIES=(
  "/path/to/spades/SAMPLE1/contigs.fasta"
  "/path/to/spades/SAMPLE2/contigs.fasta"
)
PYNGOST_DB="/path/to/pyngost_allelesDB"
OUTDIR="/path/to/output/pyngoST"
# =====================

mkdir -p "${OUTDIR}"
python3 pyngoST.py -i "${ASSEMBLIES[@]}" -p "${PYNGOST_DB}" \
  -s NG-STAR,MLST -c -m -b -a -o "${OUTDIR}/${BATCH}_pyngoST_out.txt"

echo "Done. Output: ${OUTDIR}/${BATCH}_pyngoST_out.txt"
