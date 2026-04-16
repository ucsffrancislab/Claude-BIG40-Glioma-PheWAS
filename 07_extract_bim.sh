#!/usr/bin/env bash
# =============================================================================
# 07a_extract_bim.sh
#
# Extracts a PLINK .bim file from imputed VCFs for PRS-CS.
# PRS-CS needs a BIM to know which SNPs are in the target genotype data.
# Only needs to be run once — any one cohort will do since imputed data
# shares the same variant set.
#
# Usage:
#   bash 07a_extract_bim.sh <VCF_DIR> <OUT_PREFIX>
#
#   VCF_DIR     Directory containing chr{1..22}.dose.vcf.gz
#   OUT_PREFIX  Output prefix for PLINK files (e.g., data/big40/target_bim/allchr)
#
# Example:
#   bash 07a_extract_bim.sh input/imputed-umich-cidr data/big40/target_bim/allchr
# =============================================================================

set -eu

VCF_DIR="${1:?Usage: 07a_extract_bim.sh <VCF_DIR> <OUT_PREFIX>}"
OUT_PREFIX="${2:?Usage: 07a_extract_bim.sh <VCF_DIR> <OUT_PREFIX>}"

OUT_DIR=$(dirname "$OUT_PREFIX")
mkdir -p "$OUT_DIR"

# Check if BIM already exists and is write-protected
if [ -f "${OUT_PREFIX}.bim" ] && [ ! -w "${OUT_PREFIX}.bim" ]; then
    echo "BIM already exists (write-protected): ${OUT_PREFIX}.bim"
    n=$(wc -l < "${OUT_PREFIX}.bim")
    echo "  ${n} variants"
    exit 0
fi

echo "============================================================"
echo "  Extracting BIM from VCFs"
echo "  VCF dir: ${VCF_DIR}"
echo "  Output:  ${OUT_PREFIX}.bim"
echo "============================================================"

# Process each chromosome: extract variant info only (no genotypes)
# Using plink2 --make-just-bim to avoid reading dosages
MERGE_LIST=$(mktemp)

for chr in $(seq 1 22); do
    vcf="${VCF_DIR}/chr${chr}.dose.vcf.gz"
    chr_prefix="${OUT_DIR}/tmp_chr${chr}"

    if [ ! -f "$vcf" ]; then
        echo "  WARNING: ${vcf} not found, skipping chr${chr}"
        continue
    fi

    echo "  chr${chr} ..."

    plink2 \
        --vcf "$vcf" dosage=DS \
        --make-just-bim \
        --out "$chr_prefix" \
        --threads 1 \
        --memory 2000 \
        > /dev/null 2>&1

    if [ -f "${chr_prefix}.bim" ]; then
        echo "${chr_prefix}" >> "$MERGE_LIST"
    fi
done

# Concatenate all chromosome BIMs into one
echo ""
echo "  Merging chromosomes ..."

cat /dev/null > "${OUT_PREFIX}.bim"
while IFS= read -r prefix; do
    cat "${prefix}.bim" >> "${OUT_PREFIX}.bim"
    rm -f "${prefix}.bim" "${prefix}.log"
done < "$MERGE_LIST"
rm -f "$MERGE_LIST"

n=$(wc -l < "${OUT_PREFIX}.bim")
echo "  Total variants: ${n}"

chmod a-w "${OUT_PREFIX}.bim"
echo "  Wrote: ${OUT_PREFIX}.bim"
echo "  Done."
