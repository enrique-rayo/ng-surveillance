#!/usr/bin/env bash
set -euo pipefail

# Step 3 - Taxonomic classification (Kraken2)
# Run one sample at a time - this DB is large and memory-hungry.

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FWD="/path/to/trimmomatic/${SAMPLE}_forward_paired.fq.gz"
REV="/path/to/trimmomatic/${SAMPLE}_reverse_paired.fq.gz"
KRAKEN2_DB="/path/to/kraken2_db"
OUTDIR="/path/to/output/kraken2"
THREADS=8
# =====================

mkdir -p "${OUTDIR}"
kraken2 --db "${KRAKEN2_DB}" \
  --threads "${THREADS}" --paired --gzip-compressed --report-zero-counts \
  --output "${OUTDIR}/${SAMPLE}_kraken_results.txt" \
  --report "${OUTDIR}/${SAMPLE}_kraken_report.txt" \
  --classified-out "${OUTDIR}/${SAMPLE}_classified#.fq" \
  --unclassified-out "${OUTDIR}/${SAMPLE}_unclassified#.fq" \
  "${FWD}" "${REV}"

# taxid 0 = unclassified, 482 = genus Neisseria, 485 = species N. gonorrhoeae
REPORT="${OUTDIR}/${SAMPLE}_kraken_report.txt"
UNCLASS=$(awk -F'\t' '$5==0 {print $1}' "${REPORT}"); UNCLASS=${UNCLASS:-0}
GENUS=$(awk -F'\t' '$5==482 {print $1}' "${REPORT}"); GENUS=${GENUS:-0}
SPECIES=$(awk -F'\t' '$5==485 {print $1}' "${REPORT}"); SPECIES=${SPECIES:-0}
REST=$(echo "100 - ${UNCLASS} - ${SPECIES}" | bc)

SUMMARY="${OUTDIR}/${SAMPLE}_kraken_summary.tsv"
{
  echo -e "sample\tunclassified_pct\tgenus_neisseria_pct\tspecies_gonorrhoeae_pct\trest_pct"
  echo -e "${SAMPLE}\t${UNCLASS}\t${GENUS}\t${SPECIES}\t${REST}"
} > "${SUMMARY}"

echo "Done. Summary: ${SUMMARY}"
