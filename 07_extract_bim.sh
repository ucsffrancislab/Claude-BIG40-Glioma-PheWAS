#!/usr/bin/env bash
# =============================================================================
# 07_extract_bim.sh
#
# Extracts PLINK .bim files from imputed VCFs for each cohort.
# PRS-CS needs per-cohort BIMs because different genotyping arrays and
# QC filters produce different variant sets even after imputation.
#
# Usage:
#   bash 07_extract_bim.sh <INPUT_DIR> <OUT_DIR>
#
#   INPUT_DIR  Directory containing cohort subdirs (imputed-umich-*)
#   OUT_DIR    Output directory for BIM files
#
# Example:
#   bash 07_extract_bim.sh input data/big40/target_bim
#
# Output:
#   data/big40/target_bim/imputed-umich-cidr.bim
#   data/big40/target_bim/imputed-umich-i370.bim
#   data/big40/target_bim/imputed-umich-onco.bim
#   data/big40/target_bim/imputed-umich-tcga.bim
# =============================================================================

set -eu

module load plink2

INPUT_DIR="${1:?Usage: 07_extract_bim.sh <INPUT_DIR> <OUT_DIR>}"
OUT_DIR="${2:?Usage: 07_extract_bim.sh <INPUT_DIR> <OUT_DIR>}"

mkdir -p "$OUT_DIR"

echo "============================================================"
echo "  Extracting BIM files from imputed VCFs"
echo "  Input:  ${INPUT_DIR}"
echo "  Output: ${OUT_DIR}"
echo "============================================================"
echo ""

for cohort_dir in "${INPUT_DIR}"/imputed-umich-*; do
    [ -d "$cohort_dir" ] || continue

    cohort=$(basename "$cohort_dir")
    out_prefix="${OUT_DIR}/${cohort}"

    # Write-protected = already complete
    if [ -f "${out_prefix}.bim" ] && [ ! -w "${out_prefix}.bim" ]; then
        n=$(wc -l < "${out_prefix}.bim")
        echo "SKIP  ${cohort}  (${n} variants, write-protected)"
        continue
    fi

    echo "${cohort} ..."

    # Concatenate per-chromosome BIMs
    cat /dev/null > "${out_prefix}.bim"

    for chr in $(seq 1 22); do
        vcf="${cohort_dir}/chr${chr}.dose.vcf.gz"

        if [ ! -f "$vcf" ]; then
            echo "  WARNING: ${vcf} not found, skipping chr${chr}"
            continue
        fi

        chr_prefix="${OUT_DIR}/tmp_${cohort}_chr${chr}"

        plink2 \
            --vcf "$vcf" dosage=DS \
            --make-just-bim \
            --out "$chr_prefix" \
            --threads 1 \
            --memory 2000 \
            > /dev/null 2>&1

        if [ -f "${chr_prefix}.bim" ]; then
            cat "${chr_prefix}.bim" >> "${out_prefix}.bim"
            rm -f "${chr_prefix}.bim" "${chr_prefix}.log"
        fi
    done

    n=$(wc -l < "${out_prefix}.bim")
    chmod a-w "${out_prefix}.bim"
    echo "  OK  ${cohort}  ${n} variants"
    echo ""

done

# ── Summary ──────────────────────────────────────────────────────────────────

echo "============================================================"
echo "  Summary"
echo "============================================================"

for bim in "${OUT_DIR}"/*.bim; do
    [ -f "$bim" ] || continue
    n=$(wc -l < "$bim")
    printf "  %-40s  %'d variants\n" "$(basename "$bim")" "$n"
done
