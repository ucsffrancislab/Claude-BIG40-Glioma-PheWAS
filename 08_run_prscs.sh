#!/usr/bin/env bash
# =============================================================================
# 08_run_prscs.sh
#
# SLURM array job to run PRS-CS on all 3,935 BIG40 IDPs for one cohort.
#
# Run once per cohort:
#   sbatch 08_run_prscs.sh imputed-umich-cidr
#   sbatch 08_run_prscs.sh imputed-umich-i370
#   sbatch 08_run_prscs.sh imputed-umich-onco
#   sbatch 08_run_prscs.sh imputed-umich-tcga
#
# Each array task processes one IDP across all 22 chromosomes.
#
# Before running:
#   1. Complete 06_format_for_prscs.sh (summary stats formatted)
#   2. Complete 07_extract_bim.sh (per-cohort BIMs extracted)
#   3. Adjust paths below to match your setup
# =============================================================================

#SBATCH --job-name=prscs
#SBATCH --array=1-3935
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/prscs/%x_%a.out
#SBATCH --error=logs/prscs/%x_%a.err

set -eu

# ── Configuration — EDIT THESE PATHS ─────────────────────────────────────────

PRSCS_PY="${HOME}/.local/PRScs/PRScs.py"
LD_REF="${HOME}/.local/ld_ref/ldblk_1kg_eur"
BIM_DIR="data/big40/target_bim"
SST_DIR="data/big40/prscs_input"
OUT_BASE="data/big40/prscs_output"
N_GWAS=33224

# ── Cohort from command line argument ────────────────────────────────────────

COHORT="${1:?Usage: sbatch 08_run_prscs.sh <cohort_name>}"
BIM_PREFIX="${BIM_DIR}/${COHORT}"
OUT_DIR="${OUT_BASE}/${COHORT}"

if [ ! -f "${BIM_PREFIX}.bim" ]; then
    echo "ERROR: BIM not found: ${BIM_PREFIX}.bim"
    echo "  Run 07_extract_bim.sh first."
    exit 1
fi

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
    echo "SKIP  IDP ${IDP}  (no summary stats — gap in numbering)"
    exit 0
fi

# Check if already complete (all 22 chromosome output files write-protected)
n_complete=0
for chr in $(seq 1 22); do
    outfile="${IDP_OUT_DIR}/${IDP}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
    if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
        n_complete=$((n_complete + 1))
    fi
done

if [ "$n_complete" -eq 22 ]; then
    echo "SKIP  IDP ${IDP} ${COHORT}  (all 22 chromosomes complete)"
    exit 0
fi

# ── Run PRS-CS ───────────────────────────────────────────────────────────────

mkdir -p "$IDP_OUT_DIR" "logs/prscs"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] PRS-CS  IDP=${IDP}  cohort=${COHORT}"

for chr in $(seq 1 22); do
    outfile="${IDP_OUT_DIR}/${IDP}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

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

    if [ -f "$outfile" ]; then
        chmod a-w "$outfile"
        echo "  chr${chr} ... OK"
    else
        echo "  chr${chr} ... FAIL (no output)"
    fi
done

# ── Verify ───────────────────────────────────────────────────────────────────

n_done=$(find "$IDP_OUT_DIR" -name "${IDP}_pst_eff_*_chr*.txt" ! -writable 2>/dev/null | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] IDP ${IDP} ${COHORT}: ${n_done}/22 chromosomes"
