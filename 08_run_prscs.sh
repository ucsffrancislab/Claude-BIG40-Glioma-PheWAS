#!/usr/bin/env bash
# =============================================================================
# 07b_run_prscs.sbatch
#
# SLURM array job to run PRS-CS on all 3,935 BIG40 IDPs.
#
# Each array task processes one IDP across all 22 chromosomes.
# PRS-CS runs per-chromosome internally when --chrom is not specified,
# or we can loop over chromosomes within each task.
#
# Usage:
#   sbatch 07b_run_prscs.sbatch
#
# Before running:
#   1. Complete 06_format_for_prscs.sh (summary stats formatted)
#   2. Complete 07a_extract_bim.sh (target BIM extracted)
#   3. Adjust paths below to match your setup
# =============================================================================

#SBATCH --job-name=prscs
#SBATCH --array=1-3935
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/prscs/prscs_%a.out
#SBATCH --error=logs/prscs/prscs_%a.err

set -eu

# ── Configuration — EDIT THESE PATHS ─────────────────────────────────────────

PRSCS_PY="${HOME}/.local/PRScs/PRScs.py"
LD_REF="${HOME}/.local/ld_ref/ldblk_1kg_eur"
BIM_PREFIX="data/big40/target_bim/allchr"
SST_DIR="data/big40/prscs_input"
OUT_DIR="data/big40/prscs_output"
N_GWAS=33224

# ── Thread control (important for SLURM) ─────────────────────────────────────

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# ── Determine IDP for this array task ────────────────────────────────────────

IDP=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
SST_FILE="${SST_DIR}/${IDP}.txt"
IDP_OUT_DIR="${OUT_DIR}/${IDP}"

# Check input exists
if [ ! -f "$SST_FILE" ]; then
    echo "ERROR: Summary stats not found: ${SST_FILE}"
    echo "  IDP ${IDP} may not exist (gap in numbering). Exiting."
    exit 0
fi

# Check if already complete (all 22 chromosome output files exist and are write-protected)
n_complete=0
for chr in $(seq 1 22); do
    outfile="${IDP_OUT_DIR}/${IDP}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
    if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
        n_complete=$((n_complete + 1))
    fi
done

if [ "$n_complete" -eq 22 ]; then
    echo "SKIP  IDP ${IDP}  (all 22 chromosomes complete)"
    exit 0
fi

# ── Run PRS-CS ───────────────────────────────────────────────────────────────

mkdir -p "$IDP_OUT_DIR" "logs/prscs"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting PRS-CS for IDP ${IDP}"
echo "  Summary stats: ${SST_FILE}"
echo "  LD reference:  ${LD_REF}"
echo "  BIM prefix:    ${BIM_PREFIX}"
echo "  N GWAS:        ${N_GWAS}"
echo "  Output dir:    ${IDP_OUT_DIR}"

for chr in $(seq 1 22); do
    outfile="${IDP_OUT_DIR}/${IDP}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

    # Skip chromosomes already done
    if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
        echo "  chr${chr} ... SKIP"
        continue
    fi

    echo "  chr${chr} ... running"

    python3 "$PRSCS_PY" \
        --ref_dir="$LD_REF" \
        --bim_prefix="$BIM_PREFIX" \
        --sst_file="$SST_FILE" \
        --n_gwas="$N_GWAS" \
        --chrom="$chr" \
        --out_dir="$IDP_OUT_DIR" \
        --out_name="$IDP"

    # Write-protect completed output
    if [ -f "$outfile" ]; then
        chmod a-w "$outfile"
        echo "  chr${chr} ... OK"
    else
        echo "  chr${chr} ... FAIL (no output file)"
    fi
done

# ── Verify ───────────────────────────────────────────────────────────────────

n_done=$(find "$IDP_OUT_DIR" -name "${IDP}_pst_eff_*_chr*.txt" ! -writable 2>/dev/null | wc -l)
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] IDP ${IDP} complete: ${n_done}/22 chromosomes"
