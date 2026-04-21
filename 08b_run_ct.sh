#!/usr/bin/env bash
# =============================================================================
# 08b_run_ct.sh  —  Clump + Threshold PGS for all BIG40 IDPs
#
# Fast alternative to PRS-CS: LD-clump summary stats using 1000G EUR,
# then score all 4 cohorts. Minutes per IDP instead of hours.
#
# Usage:
#   sbatch 08b_run_ct.sh <input_dir> <out_base> [start_idp] [end_idp]
#
# Examples:
#   sbatch 08b_run_ct.sh /path/to/input /francislab/data1/working/BIG40/ct_output
#   sbatch 08b_run_ct.sh /path/to/input /path/to/ct_output 1 2126
#
# <input_dir> contains the target genotype plink filesets:
#   imputed-umich-{cidr,i370,onco,tcga}.{bed,bim,fam}
#
# Crash / timeout recovery:
#   - Score files are chmod a-w on success
#   - Workers skip IDPs whose 4 cohort scores are all write-protected
#   - Partial IDPs: only missing cohort scores are recomputed
#   - Resubmit same command to resume
# =============================================================================

#SBATCH --job-name=ct_pgs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=14-0
#SBATCH --mail-type=FAIL
#SBATCH --export=None

set -eu

# ── CLI arguments ────────────────────────────────────────────────────────────
INPUT_DIR="${1:?Usage: sbatch 08b_run_ct.sh <input_dir> <out_base> [start_idp] [end_idp]}"
OUT_BASE="${2:?Usage: sbatch 08b_run_ct.sh <input_dir> <out_base> [start_idp] [end_idp]}"

# ── Configuration ────────────────────────────────────────────────────────────
EUR_REF="/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR"
SST_DIR="/francislab/data1/refs/BIG40/prscs_input"
COHORTS="cidr i370 onco tcga"
N_GWAS=33224
START_IDP="${3:-${START_IDP:-1}}"
END_IDP="${4:-${END_IDP:-3935}}"
N_WORKERS="${SLURM_CPUS_PER_TASK:-32}"

# C+T parameters
CLUMP_P1="${CLUMP_P1:-1}"        # p-value threshold (1 = keep all SNPs)
CLUMP_R2="${CLUMP_R2:-0.1}"      # LD r² clumping threshold
CLUMP_KB="${CLUMP_KB:-250}"      # clumping window in kb

# ── Load required modules ────────────────────────────────────────────────────
module load plink 2>/dev/null || true

# ── Validate inputs ──────────────────────────────────────────────────────────
if [ ! -f "${EUR_REF}.bed" ]; then
    echo "ERROR: EUR reference not found: ${EUR_REF}.bed" >&2
    exit 1
fi
for cohort in $COHORTS; do
    if [ ! -f "${INPUT_DIR}/imputed-umich-${cohort}/chr1.dose.vcf.gz" ]; then
        echo "ERROR: VCF not found: ${INPUT_DIR}/imputed-umich-${cohort}/chr1.dose.vcf.gz" >&2
        exit 1
    fi
done

# ── Directory setup ──────────────────────────────────────────────────────────
CLUMP_DIR="${OUT_BASE}/clumped"
SCORE_DIR="${OUT_BASE}/scores"
LOG_DIR="${OUT_BASE}/logs/${SLURM_JOB_ID:-local}"
SCRATCH="${TMPDIR:-/tmp}/ct_$$"
for cohort in $COHORTS; do
    mkdir -p "${SCORE_DIR}/${cohort}"
done
mkdir -p "$CLUMP_DIR" "$LOG_DIR" "$SCRATCH"

# ── Pre-filter: extract EUR SNP list + rsID→chr:pos map (once) ────────────────
EUR_SNPS="${OUT_BASE}/.eur_snps.txt"
RSID_MAP="${OUT_BASE}/.rsid_to_chrpos.txt"
if [ ! -f "$EUR_SNPS" ]; then
    awk '{print $2}' "${EUR_REF}.bim" > "${EUR_SNPS}.tmp"
    mv "${EUR_SNPS}.tmp" "$EUR_SNPS"
    echo "[$(date '+%F %T')] EUR SNP list: $(wc -l < "$EUR_SNPS") SNPs"
fi
if [ ! -f "$RSID_MAP" ]; then
    awk '{OFS="\t"; print $2, $1":"$4}' "${EUR_REF}.bim" > "${RSID_MAP}.tmp"
    mv "${RSID_MAP}.tmp" "$RSID_MAP"
    echo "[$(date '+%F %T')] rsID-to-chr:pos map: $(wc -l < "$RSID_MAP") entries"
fi

# ── Worker function ──────────────────────────────────────────────────────────
process_idp() {
    local idp="$1"
    local log_file="${LOG_DIR}/${idp}.log"

    {
        # ── Check if fully done (all 4 cohort scores write-protected) ────
        local n_done=0
        for cohort in $COHORTS; do
            local scorefile="${SCORE_DIR}/${cohort}/${idp}.profile"
            if [ -f "$scorefile" ] && [ ! -w "$scorefile" ]; then
                n_done=$((n_done + 1))
            fi
        done
        if [ "$n_done" -eq 4 ]; then
            echo "[$(date '+%F %T')] SKIP ${idp} (all 4 cohorts complete)"
            exit 0
        fi

        # ── Locate and decompress summary stats ─────────────────────────
        local sst_file
        if [ -f "${SST_DIR}/${idp}.txt.gz" ]; then
            sst_file="${SCRATCH}/${idp}.txt"
            zcat "${SST_DIR}/${idp}.txt.gz" > "$sst_file"
        elif [ -f "${SST_DIR}/${idp}.txt" ]; then
            sst_file="${SST_DIR}/${idp}.txt"
        else
            echo "[$(date '+%F %T')] MISS ${idp} (no sumstats)"
            exit 0
        fi

        local t0=$(date +%s)

        # ── Clump (once per IDP, reuse across cohorts) ───────────────────
        local clump_snps="${CLUMP_DIR}/${idp}.snps"
        if [ -f "$clump_snps" ] && [ ! -w "$clump_snps" ]; then
            echo "[$(date '+%F %T')] CLUMP ${idp} ... cached"
        else
            # Pre-filter sumstats to EUR SNPs (15M → ~1.6M lines; huge speedup)
            local sst_filt="${SCRATCH}/${idp}_filt.txt"
            awk 'NR==FNR{a[$1]; next} FNR==1 || ($1 in a)' \
                "$EUR_SNPS" "$sst_file" > "$sst_filt"

            plink --bfile "$EUR_REF" \
                  --clump "$sst_filt" \
                  --clump-snp-field SNP \
                  --clump-field P \
                  --clump-p1 "$CLUMP_P1" \
                  --clump-r2 "$CLUMP_R2" \
                  --clump-kb "$CLUMP_KB" \
                  --out "${SCRATCH}/${idp}" \
                  --threads 1 \
                  --memory 4000 \
                  > /dev/null 2>&1

            if [ -f "${SCRATCH}/${idp}.clumped" ]; then
                awk 'NR>1 && $3!="" {print $3}' "${SCRATCH}/${idp}.clumped" \
                    > "${clump_snps}.tmp"
                mv "${clump_snps}.tmp" "$clump_snps"
                chmod a-w "$clump_snps"
                local n_clumped=$(wc -l < "$clump_snps")
                echo "[$(date '+%F %T')] CLUMP ${idp} ... ${n_clumped} SNPs"
            else
                echo "[$(date '+%F %T')] CLUMP ${idp} ... FAIL (no output)"
                rm -f "$sst_filt" "${SCRATCH}/${idp}".*
                [ "$sst_file" = "${SCRATCH}/${idp}.txt" ] && rm -f "$sst_file"
                exit 1
            fi
            rm -f "$sst_filt" "${SCRATCH}/${idp}".{clumped,log,nosex}
        fi

        # ── Build chr:pos score file from clumped rsIDs ────────────────────
        local chrpos_snps="${SCRATCH}/${idp}_chrpos.snps"
        local chrpos_sst="${SCRATCH}/${idp}_chrpos.txt"

        # Translate clumped rsIDs to chr:pos
        awk 'NR==FNR{map[$1]=$2; next} ($1 in map){print map[$1]}' \
            "$RSID_MAP" "$clump_snps" > "$chrpos_snps"

        # Create chr:pos version of sumstats (header preserved, rsIDs replaced)
        awk -F'\t' -v OFS='\t' \
            'NR==FNR{map[$1]=$2; next} FNR==1{print; next} ($1 in map){$1=map[$1]; print}' \
            "$RSID_MAP" "$sst_file" > "$chrpos_sst"

        local n_chrpos
        n_chrpos=$(wc -l < "$chrpos_snps")
        echo "[$(date '+%F %T')] CHRPOS ${idp} ... ${n_chrpos} SNPs translated"

        # ── Score each cohort (per-chr VCFs, then sum) ───────────────────
        for cohort in $COHORTS; do
            local scorefile="${SCORE_DIR}/${cohort}/${idp}.profile"
            if [ -f "$scorefile" ] && [ ! -w "$scorefile" ]; then
                continue    # already done
            fi
            [ -f "$scorefile" ] && rm -f "$scorefile"

            local chr_prefix="${SCRATCH}/${idp}_${cohort}"
            local any_scored=0

            for chr in $(seq 1 22); do
                local vcf="${INPUT_DIR}/imputed-umich-${cohort}/chr${chr}.dose.vcf.gz"
                [ ! -f "$vcf" ] && continue

                plink --vcf "$vcf" \
                      --extract "$chrpos_snps" \
                      --score "$chrpos_sst" 1 2 4 header sum \
                      --out "${chr_prefix}_chr${chr}" \
                      --threads 1 \
                      --memory 4000 \
                      > /dev/null 2>&1

                [ -f "${chr_prefix}_chr${chr}.profile" ] && any_scored=1
            done

            if [ "$any_scored" -eq 1 ]; then
                # Sum SCORESUM across chromosomes per individual
                # .profile cols: FID IID PHENO CNT CNT2 SCORESUM
                awk 'FNR==1{next}
                     {key=$1"\t"$2"\t"$3; cnt[key]+=$5; score[key]+=$6}
                     END{print "FID\tIID\tPHENO\tCNT\tCNT2\tSCORESUM";
                         for(k in cnt) print k"\t"cnt[k]"\t"cnt[k]"\t"score[k]}' \
                    "${chr_prefix}"_chr*.profile > "$scorefile"
                chmod a-w "$scorefile"
            else
                echo "[$(date '+%F %T')] SCORE ${idp}/${cohort} ... FAIL (no chr scored)"
            fi

            # Cleanup per-chr intermediates
            rm -f "${chr_prefix}"_chr*.{profile,log,nosex,nopred}
        done

        # Cleanup chr:pos intermediates
        rm -f "$chrpos_snps" "$chrpos_sst"

        # ── Cleanup scratch ──────────────────────────────────────────────
        [ "$sst_file" = "${SCRATCH}/${idp}.txt" ] && rm -f "$sst_file"

        local elapsed=$(( $(date +%s) - t0 ))
        printf "[%s] DONE  %s (%dm%02ds)\n" \
               "$(date '+%F %T')" "$idp" \
               $((elapsed / 60)) $((elapsed % 60))

    } >>"$log_file" 2>&1
}
export -f process_idp

# Export everything workers need
export COHORTS SST_DIR INPUT_DIR CLUMP_DIR SCORE_DIR LOG_DIR SCRATCH
export EUR_REF EUR_SNPS RSID_MAP N_GWAS
export CLUMP_P1 CLUMP_R2 CLUMP_KB

# ── Master heartbeat ─────────────────────────────────────────────────────────
(
    t_start=$(date +%s)
    while true; do
        sleep 300
        elapsed=$(( $(date +%s) - t_start ))
        h=$(( elapsed / 3600 )); m=$(( (elapsed % 3600) / 60 ))
        n_done=0
        for cohort in $COHORTS; do
            cohort_done=$(find "${SCORE_DIR}/${cohort}" -name '*.profile' ! -writable 2>/dev/null | wc -l)
            n_done=$((n_done + cohort_done))
        done
        avg=$((n_done / 4))
        printf "[HEARTBEAT %s] elapsed %dh%02dm  avg_idps_done=%d/3935\n" \
               "$(date '+%H:%M:%S')" "$h" "$m" "$avg"
    done
) &
HEARTBEAT_PID=$!
trap 'kill $HEARTBEAT_PID 2>/dev/null || true; rm -rf "$SCRATCH"; wait 2>/dev/null || true' EXIT

# ── Dispatch ─────────────────────────────────────────────────────────────────
echo "[$(date '+%F %T')] C+T PGS  idps=${START_IDP}..${END_IDP}  workers=${N_WORKERS}"
echo "[$(date '+%F %T')]   eur_ref=${EUR_REF}"
echo "[$(date '+%F %T')]   clump: p1=${CLUMP_P1} r2=${CLUMP_R2} kb=${CLUMP_KB}"
echo "[$(date '+%F %T')]   out=${OUT_BASE}"

seq -f "%04.0f" "$START_IDP" "$END_IDP" |
    xargs -P "$N_WORKERS" -I '{}' bash -c 'process_idp "$1"' _ '{}'

# ── Final tally ──────────────────────────────────────────────────────────────
echo ""
echo "[$(date '+%F %T')] === FINAL TALLY ==="
for cohort in $COHORTS; do
    n=$(find "${SCORE_DIR}/${cohort}" -name '*.profile' ! -writable 2>/dev/null | wc -l)
    echo "  ${cohort}: ${n} IDPs scored"
done
