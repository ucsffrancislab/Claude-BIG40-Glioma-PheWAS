#!/usr/bin/env bash
# =============================================================================
# 01_download_big40.sh
#
# Downloads BIG40 brain imaging GWAS summary statistics (stats33k release)
# for all 3,935 IDPs, plus variant annotation.
#
# Crash-safe design:
#   - Before each download, checks if the file exists and is write-protected.
#   - Write-protected  →  already complete, skip.
#   - Missing or writable  →  download, then write-protect on success.
#   - Rerun this script as many times as needed after network interruptions.
#
# Usage:
#   bash 01_download_big40.sh [DATA_DIR] [NJOBS]
#
#   DATA_DIR  Base directory for all downloads  (default: data/big40)
#   NJOBS     Parallel download workers         (default: 8)
#
# Estimated sizes:
#   GWAS summary stats  ~40 GB compressed  (~1.4 TB uncompressed)
#   Variant annotation   ~270 MB compressed
#
# Reference:
#   Smith et al. (2021) Nature Neuroscience 24(5):737-745
#   https://open.win.ox.ac.uk/ukbiobank/big40/
# =============================================================================

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

# The stats33k files live on open.win; the variant table lives on open.oxcin.
STATS_URL="https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"
VARIANT_URL="https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/variants.txt.gz"

DATA_DIR="${1:-data/big40}"
NJOBS="${2:-8}"
STATS_DIR="${DATA_DIR}/stats33k"
LOG_FILE="${DATA_DIR}/download.log"

# IDP files are 4-digit zero-padded: 0001.txt.gz .. 3935.txt.gz
IDP_FIRST=1
IDP_LAST=3935

# ── Helper functions ─────────────────────────────────────────────────────────

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

# safe_download <url> <output_path>
#   Returns 0 on success or skip, nonzero on failure.
safe_download() {
    local url="$1"
    local out="$2"

    if [[ -f "$out" && ! -w "$out" ]]; then
        return 0
    fi

    if curl -fSL \
            --retry 3 \
            --retry-delay 10 \
            --connect-timeout 30 \
            --max-time 3600 \
            -o "$out" \
            "$url" 2>>"$LOG_FILE"; then
        chmod a-w "$out"
        return 0
    else
        local rc=$?
        rm -f "$out"
        return $rc
    fi
}

# ── Setup ────────────────────────────────────────────────────────────────────

mkdir -p "$STATS_DIR"
log "Download session started.  DATA_DIR=${DATA_DIR}  NJOBS=${NJOBS}"

# ── Phase 1: Variant annotation ─────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 1 / 2 :  Variant annotation  (variants.txt.gz)"
echo "============================================================"

VARIANT_FILE="${DATA_DIR}/variants.txt.gz"

if [[ -f "$VARIANT_FILE" && ! -w "$VARIANT_FILE" ]]; then
    log "variants.txt.gz already complete (write-protected)."
else
    log "Downloading variants.txt.gz (~270 MB) ..."
    if safe_download "$VARIANT_URL" "$VARIANT_FILE"; then
        log "  OK   variants.txt.gz  ($(du -h "$VARIANT_FILE" | cut -f1))"
    else
        log "  FAIL variants.txt.gz — will retry on next run"
    fi
fi

# ── Phase 2: GWAS summary statistics (parallel) ─────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 2 / 2 :  GWAS summary statistics  (${NJOBS} workers)"
echo "============================================================"

total=$(( IDP_LAST - IDP_FIRST + 1 ))
log "Downloading ${total} IDP files with ${NJOBS} parallel workers ..."

# Write a small worker script that xargs will call.
# Each invocation handles one IDP number.
WORKER=$(mktemp)
cat > "$WORKER" <<'WORKER_EOF'
#!/usr/bin/env bash
set -euo pipefail
idp="$1"
stats_url="$2"
stats_dir="$3"
log_file="$4"

outfile="${stats_dir}/${idp}.txt.gz"

# Write-protected = already complete
if [[ -f "$outfile" && ! -w "$outfile" ]]; then
    echo "SKIP  ${idp}"
    exit 0
fi

if curl -fSL \
        --retry 3 \
        --retry-delay 10 \
        --connect-timeout 30 \
        --max-time 3600 \
        -o "$outfile" \
        "${stats_url}/${idp}.txt.gz" 2>>"$log_file"; then
    chmod a-w "$outfile"
    echo "OK    ${idp}  $(du -h "$outfile" | cut -f1)"
else
    rc=$?
    rm -f "$outfile"
    if [[ $rc -eq 22 ]]; then
        echo "404   ${idp}"
    else
        echo "FAIL  ${idp}"
    fi
fi
WORKER_EOF
chmod +x "$WORKER"

# Generate 4-digit zero-padded IDP numbers, pipe to xargs for parallel exec
seq -f "%04g" "$IDP_FIRST" "$IDP_LAST" \
    | xargs -P "$NJOBS" -I {} \
        bash "$WORKER" {} "$STATS_URL" "$STATS_DIR" "$LOG_FILE"

rm -f "$WORKER"

# ── Summary (count results from filesystem) ──────────────────────────────────

echo ""
echo "============================================================"
echo "  Download Summary"
echo "============================================================"

# Count write-protected files in stats_dir = successfully downloaded
completed=$(find "$STATS_DIR" -name '*.txt.gz' ! -writable 2>/dev/null | wc -l)

log "Total IDPs requested : ${total}"
log "Files complete        : ${completed}  (write-protected in ${STATS_DIR})"
log "Remaining             : $(( total - completed ))"
echo ""

if [[ "$completed" -lt "$total" ]]; then
    log "Some files still missing (gaps in IDP numbering are normal)."
    log "Rerun to retry any network failures."
fi
