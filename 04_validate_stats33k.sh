#!/usr/bin/env bash
# =============================================================================
# 04_validate_stats33k.sh
#
# Fast validation of downloaded stats33k GWAS summary statistics.
#
# Checks per file (all fast — reads only the first ~1000 lines):
#   1. Not empty
#   2. Column count matches the first file (consistent format)
#   3. Numeric spot-check at row 1000
#
# Skips full gzip integrity check (gzip -t) since the write-protection
# download pattern already ensures complete downloads.
#
# Usage:
#   bash 04_validate_stats33k.sh [DATA_DIR]
#
#   DATA_DIR  Base directory (default: data/big40)
# =============================================================================

set -eu

DATA_DIR="${1:-data/big40}"
STATS_DIR="${DATA_DIR}/stats33k"
REPORT="${DATA_DIR}/validation_stats33k.txt"

# ── Setup ────────────────────────────────────────────────────────────────────

echo "" > "$REPORT"
log() { echo "$*" | tee -a "$REPORT"; }

log "================================================================"
log "  stats33k Validation Report"
log "  $(date '+%Y-%m-%d %H:%M:%S')"
log "  Directory: ${STATS_DIR}"
log "================================================================"
log ""

n_ok=0
n_fail=0
total=0

ref_header=""
ref_ncols=""
ref_fname=""

# ── Validate ─────────────────────────────────────────────────────────────────

for f in "${STATS_DIR}"/*.txt.gz; do
    [ -f "$f" ] || continue
    total=$((total + 1))
    fname=$(basename "$f")
    idp="${fname%.txt.gz}"
    status="OK"
    note=""

    # Read first 1000 lines once, reuse for all checks
    chunk=$(zcat "$f" | head -1000 || true)

    # 1. Empty check
    if [ -z "$chunk" ]; then
        status="FAIL"
        note="empty file"
        n_fail=$((n_fail + 1))
        log "$(printf "  %-6s %-16s %s" "$status" "$fname" "$note")"
        continue
    fi

    header=$(echo "$chunk" | head -1)
    ncols=$(echo "$header" | awk -F'\t' '{print NF}')

    # 2. Column count — compare against first file
    if [ -z "$ref_header" ]; then
        ref_header="$header"
        ref_ncols="$ncols"
        ref_fname="$fname"
    elif [ "$ncols" != "$ref_ncols" ]; then
        status="WARN"
        note="cols=${ncols} (expected ${ref_ncols})"
        n_fail=$((n_fail + 1))
    fi

    # 3. Numeric spot-check at row 1000
    spot=$(echo "$chunk" | tail -1)
    if [ -n "$spot" ]; then
        bad=$(echo "$spot" | tr '\t' '\n' \
            | grep -cvE '^-?[0-9]|^NA$|^NaN$|^na$|^nan$|^\.$|^$' || true)
        if [ "$bad" -gt 0 ]; then
            status="WARN"
            note="${note:+${note}; }${bad} non-numeric field(s)"
            n_fail=$((n_fail + 1))
        fi
    fi

    if [ "$status" = "OK" ]; then
        n_ok=$((n_ok + 1))
    fi

    printf "  %-6s %-16s %d cols  %s\n" "$status" "$fname" "$ncols" "$note"

done

# ── Summary ──────────────────────────────────────────────────────────────────

log ""
log "================================================================"
log "  Summary"
log "================================================================"
log "  Total files : ${total}"
log "  Passed      : ${n_ok}"
log "  Warnings    : ${n_fail}"
log ""

if [ -n "$ref_header" ]; then
    log "  Reference header (${ref_fname}, ${ref_ncols} columns):"
    log "  $(echo "$ref_header" | tr '\t' '\n' | cat -n || true)"
fi

log ""
log "  Report saved to: ${REPORT}"
