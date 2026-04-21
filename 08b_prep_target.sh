#!/usr/bin/env bash
# =============================================================================
# 08b_prep_target.sh  —  One-time: extract EUR-panel SNPs from imputed VCFs
#                         into small plink binary filesets for fast C+T scoring.
#
# Creates one merged .bed/.bim/.fam per cohort containing only the ~467k SNPs
# in the 1000G EUR LD reference panel. These are ~100x smaller than the full
# imputed VCFs and support random-access scoring in seconds.
#
# Usage:
#   sbatch 08b_prep_target.sh <input_dir> <out_base>
#
# Example:
#   sbatch 08b_prep_target.sh ${PWD}/input ${PWD}/ct_output
#
# Output:
#   <out_base>/target_bed/<cohort>.{bed,bim,fam}   (one per cohort)
#
# Safe to rerun — skips cohorts whose final .bed already exists.
# =============================================================================

#SBATCH --job-name=ct_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=14-0
#SBATCH --mail-type=FAIL
#SBATCH --export=None

set -eu

# ── CLI arguments ────────────────────────────────────────────────────────────
INPUT_DIR="${1:?Usage: sbatch 08b_prep_target.sh <input_dir> <out_base>}"
OUT_BASE="${2:?Usage: sbatch 08b_prep_target.sh <input_dir> <out_base>}"

# ── Configuration ────────────────────────────────────────────────────────────
EUR_REF="/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR"
COHORTS="cidr i370 onco tcga"
N_WORKERS="${SLURM_CPUS_PER_TASK:-32}"

# ── Load required modules ────────────────────────────────────────────────────
module load plink 2>/dev/null || true
module load plink2 2>/dev/null || true

# ── Setup ────────────────────────────────────────────────────────────────────
TARGET_DIR="${OUT_BASE}/target_bed"
SCRATCH="${TMPDIR:-/tmp}/ct_prep_$$"
mkdir -p "$TARGET_DIR" "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# ── Build chr:pos extract list from EUR panel (once) ─────────────────────────
EUR_CHRPOS="${OUT_BASE}/.eur_snps_chrpos.txt"
if [ ! -f "$EUR_CHRPOS" ]; then
    awk '{print $1":"$4}' "${EUR_REF}.bim" > "${EUR_CHRPOS}.tmp"
    mv "${EUR_CHRPOS}.tmp" "$EUR_CHRPOS"
    echo "[$(date '+%F %T')] EUR chr:pos list: $(wc -l < "$EUR_CHRPOS") SNPs"
fi

# ── Process each cohort ──────────────────────────────────────────────────────
for cohort in $COHORTS; do
    final="${TARGET_DIR}/${cohort}.bed"
    if [ -f "$final" ] && [ ! -w "$final" ]; then
        echo "[$(date '+%F %T')] SKIP $cohort (already done)"
        continue
    fi

    echo "[$(date '+%F %T')] START $cohort"
    t0=$(date +%s)

    chr_dir="${SCRATCH}/${cohort}"
    mkdir -p "$chr_dir"

    # ── Extract EUR SNPs from each chr VCF (parallel via xargs) ──────────
    extract_chr() {
        local chr="$1"
        local vcf="${INPUT_DIR}/imputed-umich-${cohort}/chr${chr}.dose.vcf.gz"
        [ ! -f "$vcf" ] && return 0

        plink2 --vcf "$vcf" dosage=DS \
               --set-all-var-ids @:# \
               --rm-dup force-first \
               --extract "$EUR_CHRPOS" \
               --make-bed \
               --out "${chr_dir}/chr${chr}" \
               --threads 1 \
               --memory 4000 \
               > /dev/null 2>&1

        if [ -f "${chr_dir}/chr${chr}.bed" ]; then
            local n
            n=$(wc -l < "${chr_dir}/chr${chr}.bim")
            echo "  chr${chr}: ${n} SNPs extracted"
        else
            echo "  chr${chr}: FAIL"
        fi
    }
    export -f extract_chr
    export INPUT_DIR cohort chr_dir EUR_CHRPOS

    seq 1 22 | xargs -P "$N_WORKERS" -I '{}' bash -c 'extract_chr "$1"' _ '{}'

    # ── Merge 22 chr files into one per cohort ───────────────────────────
    # plink1.9 --merge-list expects a file of prefixes (one per line)
    merge_list="${chr_dir}/merge_list.txt"
    : > "$merge_list"
    for chr in $(seq 1 22); do
        [ -f "${chr_dir}/chr${chr}.bed" ] && echo "${chr_dir}/chr${chr}" >> "$merge_list"
    done

    n_chr=$(wc -l < "$merge_list")
    echo "[$(date '+%F %T')] Merging ${n_chr} chr files for $cohort ..."

    if [ "$n_chr" -eq 0 ]; then
        echo "[$(date '+%F %T')] FAIL $cohort (no chr files to merge)"
        continue
    fi

    # plink1.9 merge-list needs at least 2 entries; if just 1, just copy
    if [ "$n_chr" -eq 1 ]; then
        prefix=$(head -1 "$merge_list")
        cp "${prefix}.bed" "${TARGET_DIR}/${cohort}.bed"
        cp "${prefix}.bim" "${TARGET_DIR}/${cohort}.bim"
        cp "${prefix}.fam" "${TARGET_DIR}/${cohort}.fam"
    else
        # First entry is the base, rest go in the merge list
        base=$(head -1 "$merge_list")
        tail -n +2 "$merge_list" > "${chr_dir}/merge_rest.txt"

        plink --bfile "$base" \
              --merge-list "${chr_dir}/merge_rest.txt" \
              --make-bed \
              --out "${TARGET_DIR}/${cohort}" \
              --threads "$N_WORKERS" \
              --memory 32000 \
              > /dev/null 2>&1
    fi

    if [ -f "${TARGET_DIR}/${cohort}.bed" ]; then
        chmod a-w "${TARGET_DIR}/${cohort}".{bed,bim,fam}
        n_snps=$(wc -l < "${TARGET_DIR}/${cohort}.bim")
        n_samples=$(wc -l < "${TARGET_DIR}/${cohort}.fam")
        elapsed=$(( $(date +%s) - t0 ))
        printf "[%s] DONE  %s (%d SNPs, %d samples, %dm%02ds)\n" \
               "$(date '+%F %T')" "$cohort" "$n_snps" "$n_samples" \
               $((elapsed / 60)) $((elapsed % 60))
    else
        echo "[$(date '+%F %T')] FAIL $cohort (merge failed)"
        # Check for error log
        [ -f "${TARGET_DIR}/${cohort}.log" ] && tail -5 "${TARGET_DIR}/${cohort}.log"
    fi

    # Cleanup chr intermediates
    rm -rf "$chr_dir"
done

echo ""
echo "[$(date '+%F %T')] === DONE ==="
for cohort in $COHORTS; do
    if [ -f "${TARGET_DIR}/${cohort}.bed" ]; then
        n=$(wc -l < "${TARGET_DIR}/${cohort}.bim")
        echo "  ${cohort}: ${n} SNPs"
    else
        echo "  ${cohort}: MISSING"
    fi
done
