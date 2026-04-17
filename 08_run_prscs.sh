#!/usr/bin/env bash
# =============================================================================
# 08_run_prscs.sh
#
# SLURM array job: run PRS-CS on all 3,935 BIG40 IDPs for one cohort.
#
# Usage:
#   sbatch --export=LOG_DIR=/abs/path/logs/prscs \
#          08_run_prscs.sh <cohort> <out_base_dir>
#
#   LOG_DIR is OPTIONAL (default: <out_base_dir>/logs) but MUST be an
#   absolute path and MUST exist before sbatch is called — SLURM opens
#   log files before the script runs, so a missing dir = silent job
#   failure with no output.
#
#   <cohort> may be given with or without the 'imputed-umich-' prefix —
#   the prefix is stripped for output-path purposes.
#
# Examples:
#   sbatch 08_run_prscs.sh cidr /francislab/data1/working/BIG40/prscs_output
#   sbatch 08_run_prscs.sh imputed-umich-i370 /francislab/data1/working/BIG40/prscs_output
#
# Each array task = one IDP, looping chr1-22 internally.
# Outputs are write-protected on success so reruns skip completed chromosomes.
#
# Inputs (NOT VCFs — VCFs are used only in step 09 for scoring):
#   - Per-cohort BIM  : ${BIM_DIR}/<cohort>.bim  (w/ or w/o imputed-umich- prefix)
#   - Summary stats   : ${SST_DIR}/<IDP>.txt[.gz]
#   - LD reference    : ${LD_REF}/
# =============================================================================

#SBATCH --job-name=prscs
#SBATCH --array=1-3935
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
# NOTE: --output / --error intentionally omitted here — pass them on the
#       sbatch command line with an absolute LOG_DIR (see Usage above),
#       otherwise SLURM fails silently if the relative path can't be opened.

set -eu

# ── Configuration ────────────────────────────────────────────────────────────
PRSCS_PY="${HOME}/.local/PRScs/PRScs.py"
LD_REF="${HOME}/.local/ld_ref/ldblk_1kg_eur"
SST_DIR="/francislab/data1/refs/BIG40/prscs_input"
BIM_DIR="/francislab/data1/refs/BIG40/target_bim"
N_GWAS=33224

# ── CLI arguments ────────────────────────────────────────────────────────────
COHORT_RAW="${1:?Usage: sbatch 08_run_prscs.sh <cohort> <out_base_dir>}"
OUT_BASE="${2:?Usage: sbatch 08_run_prscs.sh <cohort> <out_base_dir>}"

# Strip 'imputed-umich-' prefix so it doesn't pollute output paths
COHORT="${COHORT_RAW#imputed-umich-}"

# Locate BIM: try clean name first, then with prefix
if   [ -f "${BIM_DIR}/${COHORT}.bim" ]; then
    BIM_PREFIX="${BIM_DIR}/${COHORT}"
elif [ -f "${BIM_DIR}/imputed-umich-${COHORT}.bim" ]; then
    BIM_PREFIX="${BIM_DIR}/imputed-umich-${COHORT}"
else
    echo "ERROR: BIM not found for cohort '${COHORT}'. Tried:"
    echo "  ${BIM_DIR}/${COHORT}.bim"
    echo "  ${BIM_DIR}/imputed-umich-${COHORT}.bim"
    exit 1
fi

OUT_DIR="${OUT_BASE}/${COHORT}"

if [ ! -d "$SST_DIR" ]; then
    echo "ERROR: Summary stats dir not found: $SST_DIR"
    exit 1
fi

# ── Thread control (important for SLURM) ─────────────────────────────────────
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# ── Determine IDP for this array task ────────────────────────────────────────
IDP=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
IDP_OUT_DIR="${OUT_DIR}/${IDP}"

# Locate summary stats (accept .txt or .txt.gz)
if   [ -f "${SST_DIR}/${IDP}.txt.gz" ]; then
    SST_SRC="${SST_DIR}/${IDP}.txt.gz"
    SST_IS_GZ=1
elif [ -f "${SST_DIR}/${IDP}.txt" ]; then
    SST_SRC="${SST_DIR}/${IDP}.txt"
    SST_IS_GZ=0
else
    echo "SKIP  IDP ${IDP}  (no summary stats — gap in numbering)"
    exit 0
fi

# Check if already complete (all 22 chromosome outputs write-protected)
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

# ── Decompress to node-local scratch if gzipped ──────────────────────────────
# PRS-CS opens the sst file with plain open() (twice) — no native gzip support.
# $TMPDIR is SLURM's per-job local scratch, auto-cleaned at job exit.
SCRATCH="${TMPDIR:-/tmp}"
if [ "$SST_IS_GZ" -eq 1 ]; then
    SST_FILE="${SCRATCH}/${IDP}.txt"
    zcat "$SST_SRC" > "$SST_FILE"
else
    SST_FILE="$SST_SRC"
fi

# ── Run PRS-CS ───────────────────────────────────────────────────────────────
mkdir -p "$IDP_OUT_DIR"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] PRS-CS  IDP=${IDP}  cohort=${COHORT}  out=${OUT_DIR}"

for chr in $(seq 1 22); do
    outfile="${IDP_OUT_DIR}/${IDP}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

    if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
        echo "  chr${chr} ... SKIP"
        continue
    fi

    echo "  chr${chr} ... running"

    # NOTE: PRS-CS's --out_dir is actually a filename PREFIX, not a directory.
    # It gets concatenated with '_pst_eff_a1_b0.5_phiauto_chr<N>.txt'.
    # So we pass <dir>/<IDP> to produce <dir>/<IDP>_pst_eff_..._chr<N>.txt.
    python3 "$PRSCS_PY" \
        --ref_dir="$LD_REF" \
        --bim_prefix="$BIM_PREFIX" \
        --sst_file="$SST_FILE" \
        --n_gwas="$N_GWAS" \
        --chrom="$chr" \
        --out_dir="${IDP_OUT_DIR}/${IDP}"

    if [ -f "$outfile" ]; then
        chmod a-w "$outfile"
        echo "  chr${chr} ... OK"
    else
        echo "  chr${chr} ... FAIL (no output)"
    fi
done

# ── Clean up scratch (belt-and-braces; $TMPDIR is auto-removed) ──────────────
if [ "$SST_IS_GZ" -eq 1 ] && [ -f "$SST_FILE" ]; then
    rm -f "$SST_FILE"
fi

# ── Verify ───────────────────────────────────────────────────────────────────
n_done=$(find "$IDP_OUT_DIR" -name "${IDP}_pst_eff_*_chr*.txt" ! -writable 2>/dev/null | wc -l)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] IDP ${IDP} ${COHORT}: ${n_done}/22 chromosomes complete"
