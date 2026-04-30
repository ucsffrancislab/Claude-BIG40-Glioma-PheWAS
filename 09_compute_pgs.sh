#!/usr/bin/env bash
# =============================================================================
# 09_compute_pgs.sh  —  Score genotypes using PRS-CS posterior weights
#
# Reads PRS-CS output (per-IDP, per-chr weight files), concatenates into
# genome-wide weight files, and scores each cohort using plink2 --score
# against the pre-extracted target binaries (from 08b_prep_target.sh).
#
# Usage:
#   sbatch 09_compute_pgs.sh <prscs_output_dir> <target_bed_dir> <score_output_dir> [idp_list_file]
#
# Example:
#   sbatch 09_compute_pgs.sh ${PWD}/prscs_output ${PWD}/ct_output/target_bed ${PWD}/prscs_scores ${PWD}/prscs_idp_list.txt
#
# Output:
#   <score_output_dir>/<cohort>/<IDP>.sscore  (one per IDP per cohort)
# =============================================================================

#SBATCH --job-name=pgs_score
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=14-0
#SBATCH --mail-type=FAIL
#SBATCH --export=None

set -eu

# ── CLI arguments ────────────────────────────────────────────────────────────
PRSCS_DIR="${1:?Usage: 09_compute_pgs.sh <prscs_output_dir> <target_bed_dir> <score_output_dir> [idp_list_file]}"
TARGET_DIR="${2:?}"
SCORE_DIR="${3:?}"
IDP_LIST="${4:-}"

# ── Configuration ────────────────────────────────────────────────────────────
COHORTS="cidr i370 onco tcga"
N_WORKERS="${SLURM_CPUS_PER_TASK:-32}"

# ── Load modules ─────────────────────────────────────────────────────────────
module load plink2 2>/dev/null || true

# ── Setup ────────────────────────────────────────────────────────────────────
SCRATCH="${TMPDIR:-/tmp}/pgs_score_$$"
mkdir -p "$SCRATCH"
for cohort in $COHORTS; do
    mkdir -p "${SCORE_DIR}/${cohort}"
done
trap 'rm -rf "$SCRATCH"' EXIT

# ── Discover IDPs ────────────────────────────────────────────────────────────
if [ -n "$IDP_LIST" ] && [ -f "$IDP_LIST" ]; then
    echo "[$(date '+%F %T')] IDP list: $IDP_LIST ($(wc -l < "$IDP_LIST") IDPs)"
else
    # Auto-discover from prscs output (first cohort that has data)
    for cohort in $COHORTS; do
        if [ -d "${PRSCS_DIR}/${cohort}" ]; then
            IDP_LIST="${SCRATCH}/idp_list.txt"
            ls "${PRSCS_DIR}/${cohort}/" | sort > "$IDP_LIST"
            echo "[$(date '+%F %T')] Auto-discovered $(wc -l < "$IDP_LIST") IDPs from ${cohort}"
            break
        fi
    done
fi

# ── Validate target binaries ────────────────────────────────────────────────
for cohort in $COHORTS; do
    if [ ! -f "${TARGET_DIR}/${cohort}.bed" ]; then
        echo "ERROR: Target binary not found: ${TARGET_DIR}/${cohort}.bed" >&2
        exit 1
    fi
done

# ── Build rsID→chr:pos map from target BIM (once) ───────────────────────────
# PRS-CS outputs SNP IDs as rsIDs; target .bim has chr:pos IDs.
# We need to translate the weight file SNP column for --score matching.
RSID_MAP="${SCORE_DIR}/.rsid_to_chrpos.txt"
if [ ! -f "$RSID_MAP" ]; then
    # Use the PRS-CS LD reference snpinfo to build rsID → chr:pos
    SNPINFO="${HOME}/.local/ld_ref/ldblk_1kg_eur/snpinfo_1kg_hm3"
    awk '{OFS="\t"; print $2, $1":"$3}' "$SNPINFO" > "${RSID_MAP}.tmp"
    mv "${RSID_MAP}.tmp" "$RSID_MAP"
    echo "[$(date '+%F %T')] rsID→chr:pos map: $(wc -l < "$RSID_MAP") entries"
fi

# ── Worker function ──────────────────────────────────────────────────────────
score_idp() {
    local idp="$1"

    # Check if all cohorts already scored
    local all_done=1
    for cohort in $COHORTS; do
        local sf="${SCORE_DIR}/${cohort}/${idp}.sscore"
        if [ ! -f "$sf" ] || [ -w "$sf" ]; then
            all_done=0
            break
        fi
    done
    if [ "$all_done" -eq 1 ]; then
        return 0
    fi

    # Find the PRS-CS weight files (check first available cohort)
    local weight_dir=""
    for cohort in $COHORTS; do
        if [ -d "${PRSCS_DIR}/${cohort}/${idp}" ]; then
            weight_dir="${PRSCS_DIR}/${cohort}/${idp}"
            break
        fi
    done
    if [ -z "$weight_dir" ]; then
        echo "[$(date '+%F %T')] MISS $idp (no PRS-CS output)"
        return 0
    fi

    # Check all 22 chr files exist
    local n_chr=0
    for chr in $(seq 1 22); do
        [ -f "${weight_dir}/${idp}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt" ] && n_chr=$((n_chr + 1))
    done
    if [ "$n_chr" -lt 22 ]; then
        echo "[$(date '+%F %T')] SKIP $idp (only ${n_chr}/22 chr weight files)"
        return 0
    fi

    # Concatenate all chr weight files and translate rsID → chr:pos
    # PRS-CS output format: CHR  SNP  BP  A1  A2  POSTERIOR_EFFECT (space-delimited)
    local weight_file="${SCRATCH}/${idp}_weights.txt"
    cat "${weight_dir}/${idp}_pst_eff_a1_b0.5_phiauto_chr"*.txt | \
        awk 'NR==FNR{map[$1]=$2; next} ($2 in map){$2=map[$2]; print}' \
            "$RSID_MAP" - | \
        awk '{OFS="\t"; print $2, $4, $6}' > "$weight_file"
    # Output: SNP(chr:pos)  A1  POSTERIOR_EFFECT

    local n_snps
    n_snps=$(wc -l < "$weight_file")

    # Score each cohort
    for cohort in $COHORTS; do
        local sf="${SCORE_DIR}/${cohort}/${idp}.sscore"
        if [ -f "$sf" ] && [ ! -w "$sf" ]; then
            continue
        fi
        [ -f "$sf" ] && rm -f "$sf"

        plink2 --bfile "${TARGET_DIR}/${cohort}" \
               --score "$weight_file" 1 2 3 header-read cols=+scoresums,-scoreavgs \
               --out "${SCRATCH}/${idp}_${cohort}" \
               --threads 1 \
               --memory 4000 \
               > /dev/null 2>&1

        if [ -f "${SCRATCH}/${idp}_${cohort}.sscore" ]; then
            mv "${SCRATCH}/${idp}_${cohort}.sscore" "$sf"
            chmod a-w "$sf"
        else
            echo "[$(date '+%F %T')] SCORE ${idp}/${cohort} ... FAIL"
        fi
        rm -f "${SCRATCH}/${idp}_${cohort}".log
    done

    rm -f "$weight_file"
    echo "[$(date '+%F %T')] DONE  $idp (${n_snps} SNPs, 4 cohorts)"
}
export -f score_idp

export PRSCS_DIR TARGET_DIR SCORE_DIR SCRATCH COHORTS RSID_MAP N_WORKERS

# ── Dispatch ─────────────────────────────────────────────────────────────────
n_idps=$(wc -l < "$IDP_LIST")
echo "[$(date '+%F %T')] Scoring ${n_idps} IDPs × ${#COHORTS} cohorts, ${N_WORKERS} workers"

cat "$IDP_LIST" | xargs -P "$N_WORKERS" -I '{}' bash -c 'score_idp "$1"' _ '{}'

# ── Final tally ──────────────────────────────────────────────────────────────
echo ""
echo "[$(date '+%F %T')] === FINAL TALLY ==="
for cohort in $COHORTS; do
    n=$(find "${SCORE_DIR}/${cohort}" -name '*.sscore' ! -writable 2>/dev/null | wc -l)
    echo "  ${cohort}: ${n} IDPs scored"
done
