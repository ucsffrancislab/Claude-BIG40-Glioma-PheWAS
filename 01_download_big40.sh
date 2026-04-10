#!/usr/bin/env bash
# =============================================================================
# 01_download_big40.sh
#
# Downloads BIG40 reference files: variant annotation table.
#
# The stats33k GWAS files are commented out below — the Oxford server is
# extremely slow (~500 KB/s, no parallel connections).  Use 02_download_ebi.sh
# to grab the 22k discovery stats from EBI instead, and optionally meta-analyse
# with the 11k replication stats later.
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
# Reference:
#   Smith et al. (2021) Nature Neuroscience 24(5):737-745
#   https://open.win.ox.ac.uk/ukbiobank/big40/
# =============================================================================

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

VARIANT_URL="https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/variants.txt.gz"

DATA_DIR="${1:-data/big40}"
LOG_FILE="${DATA_DIR}/download.log"

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

mkdir -p "$DATA_DIR"
log "Download session started.  DATA_DIR=${DATA_DIR}"

# ── Variant annotation ──────────────────────────────────────────────────────
#
# Shared SNP table for all IDPs: chr rsid pos a1 a2 af info
# This is needed regardless of whether we use stats33k, stats (22k), or repro (11k).

echo ""
echo "============================================================"
echo "  Variant annotation  (variants.txt.gz, ~270 MB)"
echo "============================================================"

VARIANT_FILE="${DATA_DIR}/variants.txt.gz"

if [[ -f "$VARIANT_FILE" && ! -w "$VARIANT_FILE" ]]; then
    log "variants.txt.gz already complete (write-protected)."
else
    log "Downloading variants.txt.gz ..."
    log "  URL: ${VARIANT_URL}"
    if safe_download "$VARIANT_URL" "$VARIANT_FILE"; then
        log "  OK   variants.txt.gz  ($(du -h "$VARIANT_FILE" | cut -f1))"
    else
        log "  FAIL variants.txt.gz — will retry on next run"
    fi
fi

# ── Stats33k GWAS summary statistics (DISABLED) ─────────────────────────────
#
# The Oxford server throttles downloads to ~500 KB/s per connection and blocks
# parallel connections.  At ~10 min per file × 3,935 files ≈ 1 month.
#
# Instead, use 02_download_ebi.sh to get the 22k discovery stats from the EBI
# GWAS Catalog FTP.  If you need full 33k power, you can also download the 11k
# replication stats from release2/repro/ and meta-analyse with the 22k.
#
# To re-enable: uncomment the block below and run this script.
#
# STATS_URL="https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"
# STATS_DIR="${DATA_DIR}/stats33k"
# mkdir -p "$STATS_DIR"
#
# for (( n=1; n<=3935; n++ )); do
#     idp=$(printf "%04d" "$n")
#     outfile="${STATS_DIR}/${idp}.txt.gz"
#     if [[ -f "$outfile" && ! -w "$outfile" ]]; then
#         continue
#     fi
#     printf "  [%04d / 3935]  IDP %s  downloading ...\n" "$n" "$idp"
#     if safe_download "${STATS_URL}/${idp}.txt.gz" "$outfile"; then
#         printf "  [%04d / 3935]  IDP %s  OK  (%s)\n" "$n" "$idp" "$(du -h "$outfile" | cut -f1)"
#     else
#         printf "  [%04d / 3935]  IDP %s  FAIL\n" "$n" "$idp"
#     fi
# done

log "Done."
