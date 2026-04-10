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
#   bash 01_download_big40.sh [DATA_DIR]
#
#   DATA_DIR  Base directory for all downloads (default: data/big40)
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
STATS_DIR="${DATA_DIR}/stats33k"
LOG_FILE="${DATA_DIR}/download.log"

# IDP files are 4-digit zero-padded: 0001.txt.gz .. 3935.txt.gz
# Some numbers in this range may not exist on the server (gaps) — that's normal.
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

    # Write-protected = already complete
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
log "Download session started.  DATA_DIR=${DATA_DIR}"

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
    log "  URL: ${VARIANT_URL}"
    if safe_download "$VARIANT_URL" "$VARIANT_FILE"; then
        log "  OK   variants.txt.gz  ($(du -h "$VARIANT_FILE" | cut -f1))"
    else
        log "  FAIL variants.txt.gz — will retry on next run"
    fi
fi

# ── Phase 2: GWAS summary statistics ────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 2 / 2 :  GWAS summary statistics  (${IDP_FIRST}..${IDP_LAST})"
echo "============================================================"

total=$(( IDP_LAST - IDP_FIRST + 1 ))
i=0
downloaded=0
skipped=0
not_found=0
failed=0

for (( n=IDP_FIRST; n<=IDP_LAST; n++ )); do
    i=$((i + 1))
    idp=$(printf "%04d" "$n")
    outfile="${STATS_DIR}/${idp}.txt.gz"

    # Already complete
    if [[ -f "$outfile" && ! -w "$outfile" ]]; then
        skipped=$((skipped + 1))
        printf "\r  [%4d / %d]  IDP %s  SKIP" "$i" "$total" "$idp"
        continue
    fi

    printf "\r  [%4d / %d]  IDP %s  downloading ...          " "$i" "$total" "$idp"

    if safe_download "${STATS_URL}/${idp}.txt.gz" "$outfile"; then
        downloaded=$((downloaded + 1))
        printf "\r  [%4d / %d]  IDP %s  OK   (%s)              \n" \
            "$i" "$total" "$idp" "$(du -h "$outfile" | cut -f1)"
    else
        rc=$?
        if [[ $rc -eq 22 ]]; then
            not_found=$((not_found + 1))
            printf "\r  [%4d / %d]  IDP %s  NOT FOUND (gap)       \n" \
                "$i" "$total" "$idp"
        else
            failed=$((failed + 1))
            printf "\r  [%4d / %d]  IDP %s  FAIL (curl rc=%d)     \n" \
                "$i" "$total" "$idp" "$rc"
        fi
    fi
done

echo ""

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Download Summary"
echo "============================================================"
log "Total IDPs  : ${total}"
log "Downloaded  : ${downloaded}  (this run)"
log "Skipped     : ${skipped}  (already complete)"
log "Not found   : ${not_found}  (gaps in IDP numbering)"
log "Failed      : ${failed}  (network errors — rerun to retry)"
echo ""

if [[ "$failed" -gt 0 ]]; then
    log "Some downloads failed due to network errors. Rerun this script to retry."
    exit 1
else
    log "All downloads complete."
    exit 0
fi
