#!/usr/bin/env bash
set -euo pipefail

# Step 11 - Assembly-mapped variant calling (Snippy, per sample)
# Different analysis from step 9: each sample is mapped against its own
# assembly instead of the reference genome, with different thresholds.
# Feeds into the heterozygous-site QC (steps 12-14).

# ==== EDIT THESE ====
SAMPLE="SAMPLE_ID"
FWD="/path/to/trimmomatic/${SAMPLE}_forward_paired.fq.gz"
REV="/path/to/trimmomatic/${SAMPLE}_reverse_paired.fq.gz"
ASSEMBLY="/path/to/spades/${SAMPLE}/contigs.fasta"
OUTDIR="/path/to/output/assembly_map"
THREADS=16
# =====================

snippy --outdir "${OUTDIR}/${SAMPLE}_snippy" --ref "${ASSEMBLY}" \
  --R1 "${FWD}" --R2 "${REV}" \
  --prefix "${SAMPLE}" --mincov 10 --minfrac 0.05 --unmapped --cpus "${THREADS}"

VCF="${OUTDIR}/${SAMPLE}_snippy/${SAMPLE}.vcf"
bgzip -f "${VCF}"
bcftools index -f "${VCF}.gz"

echo "Done. Output: ${OUTDIR}/${SAMPLE}_snippy"
