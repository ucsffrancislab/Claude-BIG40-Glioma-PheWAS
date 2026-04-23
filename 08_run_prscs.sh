#!/usr/bin/env bash
# =============================================================================
# 08_run_prscs.sh  (v8 — single big job, xargs parallelism, full resume)
#
# One SLURM job per cohort. Internally runs N_WORKERS concurrent PRS-CS
# invocations via xargs -P. Each invocation processes one IDP, all 22
# chromosomes in a single python call (amortizes sumstats/BIM parsing).
#
# Usage:
#   sbatch 08_run_prscs.sh <out_base> <cohort> [start_idp] [end_idp]
#
# Examples:
#   sbatch 08_run_prscs.sh /francislab/data1/working/BIG40/prscs_output cidr
#   sbatch 08_run_prscs.sh /francislab/data1/working/BIG40/prscs_output i370 1 2000
#
# Crash / timeout recovery:
#   - Output files chmod a-w on successful completion
#   - Worker skips any IDP whose 22 outputs are all write-protected
#   - Partial IDPs: worker runs only the missing chromosomes
#   - flock-based atomic claim prevents concurrent workers (or jobs) from
#     double-running the same IDP; locks auto-release on process death (even
#     SIGKILL) because the kernel drops fd-backed locks
#
# Resubmit to resume: same command line, no flags — it just picks up.
# =============================================================================

#SBATCH --job-name=prscs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=14-0
#SBATCH --mail-type=FAIL
#SBATCH --export=None

set -eu

# ── CLI arguments ────────────────────────────────────────────────────────────
OUT_BASE="${1:?Usage: sbatch 08_run_prscs.sh <out_base> <cohort> [start_idp] [end_idp]  OR  <out_base> <cohort> <idp_list_file>}"
COHORT_RAW="${2:?}"
# Accept either (start end) range or a file of IDP numbers
if [ -f "${3:-}" ]; then
    IDP_LIST_FILE="$3"
    echo "[$(date '+%F %T')] IDP list file: $IDP_LIST_FILE ($(wc -l < "$IDP_LIST_FILE") IDPs)"
else
    START_IDP="${3:-1}"
    END_IDP="${4:-3935}"
    IDP_LIST_FILE=""
fi

COHORT="${COHORT_RAW#imputed-umich-}"

# ── Configuration ────────────────────────────────────────────────────────────
PRSCS_PY="${HOME}/.local/PRScs/PRScs.py"
LD_REF="${HOME}/.local/ld_ref/ldblk_1kg_eur"
SST_DIR="/francislab/data1/refs/BIG40/prscs_input"
BIM_DIR="/francislab/data1/refs/BIG40/target_bim"
N_GWAS=33224

# Number of parallel workers (edit here to change; --export=None blocks env overrides)
N_WORKERS="${SLURM_CPUS_PER_TASK:-32}"

# MCMC iterations (default 500/250 for screening; use 1000/500 for final)
N_ITER="${N_ITER:-500}"
N_BURNIN="${N_BURNIN:-250}"

# Locate BIM (with or without imputed-umich- prefix)
if   [ -f "${BIM_DIR}/${COHORT}.bim" ];                 then BIM_PREFIX="${BIM_DIR}/${COHORT}"
elif [ -f "${BIM_DIR}/imputed-umich-${COHORT}.bim" ];   then BIM_PREFIX="${BIM_DIR}/imputed-umich-${COHORT}"
else
    echo "ERROR: BIM not found for cohort '${COHORT}' in $BIM_DIR" >&2
    exit 1
fi

# ── Pre-flight: verify PRS-CS NumPy 2.x patch is applied ─────────────────────
MCMC_PY="${PRSCS_PY%/*}/mcmc_gtb.py"
missing=0
for marker in \
    "float(np.sum(beta*beta_mrg))" \
    "float(2.0*delta[jj,0])" \
    "float(np.sum(delta))" \
    "beta_est.flatten()" \
    "psi_est.flatten()"
do
    if ! grep -qF -- "$marker" "$MCMC_PY"; then
        echo "ERROR: PRS-CS NumPy 2.x patch missing marker: $marker" >&2
        missing=$((missing + 1))
    fi
done
if [ "$missing" -gt 0 ]; then
    echo "ERROR: PRS-CS at $MCMC_PY is not patched (${missing}/5 markers missing)." >&2
    echo "       Run: bash patch_prscs_numpy2.sh" >&2
    exit 1
fi

# ── Paths / setup ────────────────────────────────────────────────────────────
COHORT_OUT="${OUT_BASE}/${COHORT}"
LOCK_DIR="${OUT_BASE}/.locks/${COHORT}"
WORKER_LOG_DIR="${OUT_BASE}/logs/workers/${COHORT}/${SLURM_JOB_ID:-local}"
mkdir -p "$COHORT_OUT" "$LOCK_DIR" "$WORKER_LOG_DIR"

# ── Threading control (each worker = 1 thread; we parallelize across workers) ─
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONUNBUFFERED=1

# Export everything worker subshells need
export COHORT BIM_PREFIX PRSCS_PY LD_REF SST_DIR COHORT_OUT LOCK_DIR
export WORKER_LOG_DIR N_GWAS N_ITER N_BURNIN

# ── Worker function ──────────────────────────────────────────────────────────
# One call = one IDP. All 22 chrs in a single python invocation (only missing
# chrs if resuming). flock-claims the IDP so multiple workers/jobs can't double-
# run it. Writes to a per-IDP log file.
process_idp() {
    local idp="$1"
    local idp_out_dir="${COHORT_OUT}/${idp}"
    local log_file="${WORKER_LOG_DIR}/${idp}.log"
    local lock_file="${LOCK_DIR}/${idp}.lock"

    mkdir -p "$idp_out_dir"

    # All output/logging inside the flock guard so concurrent attempts are silent
    (
        # fd 9 → lock file; non-blocking try-lock, exit silently if held
        exec 9>"$lock_file"
        flock -n 9 || exit 0

        {
            # ── Determine which chrs still need work ─────────────────────────
            local needed=()
            for chr in $(seq 1 22); do
                local outfile="${idp_out_dir}/${idp}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
                if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
                    continue                       # finished, protected
                fi
                if [ -f "$outfile" ] && [ -w "$outfile" ]; then
                    rm -f "$outfile"               # stale partial from a crash
                fi
                needed+=("$chr")
            done

            if [ ${#needed[@]} -eq 0 ]; then
                echo "[$(date '+%F %T')] SKIP ${COHORT}/${idp} (all 22 chrs complete)"
                exit 0
            fi

            # ── Locate sumstats ──────────────────────────────────────────────
            local sst_src sst_file
            if   [ -f "${SST_DIR}/${idp}.txt.gz" ]; then
                sst_src="${SST_DIR}/${idp}.txt.gz"
                sst_file="${TMPDIR:-/tmp}/${idp}_$$.txt"
                zcat "$sst_src" > "$sst_file"
            elif [ -f "${SST_DIR}/${idp}.txt" ]; then
                sst_src=""
                sst_file="${SST_DIR}/${idp}.txt"
            else
                echo "[$(date '+%F %T')] MISS ${COHORT}/${idp} (no sumstats — numbering gap)"
                exit 0
            fi

            # ── Run PRS-CS (single invocation for all needed chrs) ───────────
            local chrom_arg
            chrom_arg=$(IFS=,; echo "${needed[*]}")
            local n_needed=${#needed[@]}

            echo "[$(date '+%F %T')] START ${COHORT}/${idp} (${n_needed} chrs: ${chrom_arg})"
            local t0=$(date +%s)
            local prscs_rc=0

            python3 -u "$PRSCS_PY" \
                --ref_dir="$LD_REF" \
                --bim_prefix="$BIM_PREFIX" \
                --sst_file="$sst_file" \
                --n_gwas="$N_GWAS" \
                --n_iter="$N_ITER" \
                --n_burnin="$N_BURNIN" \
                --chrom="$chrom_arg" \
                --out_dir="${idp_out_dir}/${idp}" \
                || prscs_rc=$?

            local t_elapsed=$(( $(date +%s) - t0 ))

            # Cleanup scratch sumstats (if we decompressed)
            [ -n "$sst_src" ] && rm -f "$sst_file" 2>/dev/null || true

            # ── Write-protect each successful output ─────────────────────────
            local n_ok=0
            for chr in "${needed[@]}"; do
                local outfile="${idp_out_dir}/${idp}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
                if [ -f "$outfile" ]; then
                    chmod a-w "$outfile"
                    n_ok=$((n_ok + 1))
                fi
            done

            if [ "$n_ok" -eq "$n_needed" ]; then
                printf "[%s] DONE  %s/%s (%d chrs, %dh%02dm)\n" \
                       "$(date '+%F %T')" "$COHORT" "$idp" "$n_ok" \
                       $((t_elapsed / 3600)) $(((t_elapsed % 3600) / 60))
            else
                printf "[%s] FAIL  %s/%s (%d/%d chrs, rc=%d)\n" \
                       "$(date '+%F %T')" "$COHORT" "$idp" "$n_ok" \
                       "$n_needed" "$prscs_rc"
                exit 1
            fi
        } >>"$log_file" 2>&1
    )
    # flock auto-releases on subshell exit (even SIGKILL)
}
export -f process_idp

# ── Master heartbeat (every 5 min: overall progress) ─────────────────────────
(
    t_start=$(date +%s)
    while true; do
        sleep 300
        elapsed=$(( $(date +%s) - t_start ))
        h=$(( elapsed / 3600 )); m=$(( (elapsed % 3600) / 60 ))
        # Count IDPs with chr22 output write-protected (i.e. fully done)
        n_done=$(find "$COHORT_OUT" -name '*_pst_eff_*_chr22.txt' ! -writable 2>/dev/null | wc -l)
        # Active workers = flock-held lock files
        n_active=$(find "$LOCK_DIR" -name '*.lock' 2>/dev/null | wc -l)
        printf "[HEARTBEAT %s] elapsed %dh%02dm  done=%d  locks_held=%d\n" \
               "$(date '+%H:%M:%S')" "$h" "$m" "$n_done" "$n_active"
    done
) &
HEARTBEAT_PID=$!
trap 'kill $HEARTBEAT_PID 2>/dev/null || true; wait 2>/dev/null || true' EXIT

# ── Dispatch IDPs to worker pool ─────────────────────────────────────────────
if [ -n "$IDP_LIST_FILE" ]; then
    n_idps=$(wc -l < "$IDP_LIST_FILE")
    echo "[$(date '+%F %T')] START cohort=${COHORT}  idps=${n_idps} from list  workers=${N_WORKERS}"
else
    n_idps=$(( END_IDP - START_IDP + 1 ))
    echo "[$(date '+%F %T')] START cohort=${COHORT}  idps=${START_IDP}..${END_IDP}  workers=${N_WORKERS}"
fi
echo "[$(date '+%F %T')]   out=${COHORT_OUT}"
echo "[$(date '+%F %T')]   worker_logs=${WORKER_LOG_DIR}"

if [ -n "$IDP_LIST_FILE" ]; then
    cat "$IDP_LIST_FILE" | xargs -P "$N_WORKERS" -I '{}' bash -c 'process_idp "$1"' _ '{}'
else
    seq -f "%04.0f" "$START_IDP" "$END_IDP" |
        xargs -P "$N_WORKERS" -I '{}' bash -c 'process_idp "$1"' _ '{}'
fi

# ── Final verify ─────────────────────────────────────────────────────────────
n_complete=$(find "$COHORT_OUT" -name '*_pst_eff_*_chr22.txt' ! -writable 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] FINISH cohort=${COHORT}  fully_complete=${n_complete}/${n_idps}"
