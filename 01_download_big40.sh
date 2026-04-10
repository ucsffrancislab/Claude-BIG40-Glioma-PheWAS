#!/usr/bin/env bash
# =============================================================================
# 01_download_big40.sh
#
# Downloads BIG40 brain imaging GWAS summary statistics (stats33k release)
# for all 3,935 IDPs, plus variant annotation and IDP metadata.
#
# Crash-safe design:
#   - Before each download, checks if the file exists and is write-protected.
#   - Write-protected  →  skip  (already downloaded successfully).
#   - Missing or writable  →  download, then write-protect on success.
#   - Rerun this script as many times as needed to complete all downloads.
#
# Usage:
#   bash 01_download_big40.sh [DATA_DIR]
#
#   DATA_DIR  Base directory for all downloads (default: data/big40)
#
# Estimated sizes:
#   GWAS summary stats  ~1.4 TB uncompressed (~40 GB compressed)
#   SNP annotation       ~270 MB compressed
#
# Reference:
#   Smith et al. (2021) Nature Neuroscience 24(5):737-745
#   https://open.win.ox.ac.uk/ukbiobank/big40/
# =============================================================================

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

BASE_URL="https://open.win.ox.ac.uk/ukbiobank/big40"
STATS_URL="${BASE_URL}/release2/stats33k"

DATA_DIR="${1:-data/big40}"
STATS_DIR="${DATA_DIR}/stats33k"
META_DIR="${DATA_DIR}/metadata"
LOG_FILE="${DATA_DIR}/download.log"

# ── Helper functions ─────────────────────────────────────────────────────────

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

# safe_download <url> <output_path>
#   Downloads a file with the write-protection crash-safe pattern.
#   Returns 0 on success or skip, 1 on failure.
safe_download() {
    local url="$1"
    local out="$2"
    local fname
    fname="$(basename "$out")"

    # Already downloaded: write-protected files are complete
    if [[ -f "$out" && ! -w "$out" ]]; then
        return 0
    fi

    # Download (or re-download if a previous attempt left a partial file)
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
        rm -f "$out"          # remove partial file
        return 1
    fi
}

# ── Create directory structure ───────────────────────────────────────────────

mkdir -p "$STATS_DIR" "$META_DIR"
log "Download session started.  DATA_DIR=${DATA_DIR}"

# ── Phase 1: Metadata ───────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 1 / 3 :  Metadata"
echo "============================================================"

# IDP descriptions (HTML table mapping IDP numbers → phenotype names)
log "Downloading IDP metadata table ..."
safe_download "${BASE_URL}/IDPs.html" "${META_DIR}/IDPs.html" \
    && log "  OK   IDPs.html" \
    || log "  FAIL IDPs.html"

# SNP annotation file: allele frequencies, INFO scores (~270 MB compressed)
log "Downloading SNP annotation file (~270 MB compressed) ..."
safe_download "${STATS_URL}/snp_stats33k.txt.gz" "${META_DIR}/snp_stats33k.txt.gz" \
    && log "  OK   snp_stats33k.txt.gz" \
    || log "  FAIL snp_stats33k.txt.gz"

# ── Phase 2: Build IDP number list ──────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 2 / 3 :  Build IDP list"
echo "============================================================"

IDP_LIST="${META_DIR}/idp_numbers.txt"

if [[ -f "$IDP_LIST" && ! -w "$IDP_LIST" ]]; then
    log "IDP list already exists (write-protected). Using existing list."
else
    log "Generating IDP number list from server directory listing ..."

    # Strategy: fetch the stats33k directory index and extract filenames
    TMPFILE=$(mktemp)

    if curl -fsSL "${STATS_URL}/" -o "$TMPFILE" 2>>"$LOG_FILE"; then
        grep -oE '[0-9]+\.txt\.gz' "$TMPFILE" \
            | sed 's/\.txt\.gz//' \
            | sort -un \
            > "$IDP_LIST"
    fi
    rm -f "$TMPFILE"

    n_idps=0
    [[ -f "$IDP_LIST" ]] && n_idps=$(wc -l < "$IDP_LIST")

    # Fallback: parse IDP numbers from the metadata HTML table
    if [[ "$n_idps" -lt 100 && -f "${META_DIR}/IDPs.html" ]]; then
        log "  Directory listing yielded ${n_idps} entries; trying IDPs.html ..."
        grep -oE 'IDP_[0-9]+' "${META_DIR}/IDPs.html" \
            | sed 's/IDP_//' \
            | sort -un \
            > "$IDP_LIST"
        n_idps=$(wc -l < "$IDP_LIST")
    fi

    # Last resort: use the expected numeric range
    if [[ "$n_idps" -lt 100 ]]; then
        log "  Auto-detection found ${n_idps} IDPs; falling back to range 0..3934."
        log "  Edit ${IDP_LIST} manually if your IDP numbering differs."
        seq 0 3934 > "$IDP_LIST"
        n_idps=3935
    fi

    chmod a-w "$IDP_LIST"
    log "IDP list written: ${IDP_LIST}  (${n_idps} IDPs)"
fi

# ── Phase 3: Download GWAS summary statistics ───────────────────────────────

echo ""
echo "============================================================"
echo "  Phase 3 / 3 :  GWAS summary statistics"
echo "============================================================"

total=$(wc -l < "$IDP_LIST")
i=0
downloaded=0
skipped=0
failed=0

while IFS= read -r idp_num; do
    i=$((i + 1))

    outfile="${STATS_DIR}/${idp_num}.txt.gz"

    # Decide action BEFORE calling safe_download (cleaner counting)
    if [[ -f "$outfile" && ! -w "$outfile" ]]; then
        skipped=$((skipped + 1))
        printf "\r  [%4d / %d]  IDP %-6s  SKIP" "$i" "$total" "$idp_num"
        continue
    fi

    printf "\r  [%4d / %d]  IDP %-6s  downloading ..." "$i" "$total" "$idp_num"

    if safe_download "${STATS_URL}/${idp_num}.txt.gz" "$outfile"; then
        downloaded=$((downloaded + 1))
        printf "\r  [%4d / %d]  IDP %-6s  OK   (%s)\n" \
            "$i" "$total" "$idp_num" "$(du -h "$outfile" | cut -f1)"
    else
        failed=$((failed + 1))
        printf "\r  [%4d / %d]  IDP %-6s  FAIL\n" "$i" "$total" "$idp_num"
    fi

done < "$IDP_LIST"

# Clear the last \r line if we ended on a skip
echo ""

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Download Summary"
echo "============================================================"
log "Total IDPs : ${total}"
log "Downloaded : ${downloaded}  (this run)"
log "Skipped    : ${skipped}  (already complete)"
log "Failed     : ${failed}"
echo ""

if [[ "$failed" -gt 0 ]]; then
    log "Some downloads failed. Rerun this script to retry."
    exit 1
else
    log "All downloads complete."
    exit 0
fi
