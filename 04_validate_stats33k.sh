#!/usr/bin/env bash
# =============================================================================
# 04_validate_stats33k.sh
#
# Fast validation of downloaded stats33k GWAS summary statistics.
#
# Checks per file:
#   1. gzip integrity  (gzip -t — reads full file but no pipe overhead)
#   2. Column count    (header should have 8 fields: beta se neglog10p per SNP,
#                       but the variant info is in variants.txt.gz so each
#                       per-IDP file may just have 3 columns)
#   3. Numeric spot-check (row 1000 — are the values actually numbers?)
#
# Deliberately skips full row count (zcat | wc -l) since that was taking
# minutes per file for a marginal check.
#
# Usage:
#   bash 04_validate_stats33k.sh [DATA_DIR]
#
#   DATA_DIR  Base directory (default: data/big40)
# =============================================================================

set -euo pipefail

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
n_gzip_fail=0
n_header_fail=0
n_numeric_fail=0
n_empty=0
total=0

# Grab the header from the first valid file to use as reference
ref_header=""

# ── Validate ─────────────────────────────────────────────────────────────────

for f in "${STATS_DIR}"/*.txt.gz; do
    [ -f "$f" ] || continue
    total=$((total + 1))
    fname=$(basename "$f")
    idp="${fname%.txt.gz}"
    errors=""

    # 1. gzip integrity
    if ! gzip -t "$f" 2>/dev/null; then
        n_gzip_fail=$((n_gzip_fail + 1))
        log "FAIL  ${fname}  gzip corrupt"
        continue
    fi

    # 2. Header / column check
    header=$(zcat "$f" | head -1 || true)

    if [ -z "$header" ]; then
        n_empty=$((n_empty + 1))
        log "FAIL  ${fname}  empty file"
        continue
    fi

    ncols=$(echo "$header" | awk -F'\t' '{print NF}')

    # Store first header as reference; flag if others differ
    if [ -z "$ref_header" ]; then
        ref_header="$header"
        ref_ncols="$ncols"
        ref_fname="$fname"
    elif [ "$ncols" != "$ref_ncols" ]; then
        errors="${errors}  WARN column count ${ncols} differs from ${ref_fname} (${ref_ncols})\n"
        n_header_fail=$((n_header_fail + 1))
    fi

    # 3. Numeric spot-check: row 1000
    spot=$(zcat "$f" | head -1000 | tail -1 || true)
    if [ -n "$spot" ]; then
        # Check that all fields look numeric (allow scientific notation, NA, -, .)
        bad_fields=$(echo "$spot" | tr '\t' '\n' \
            | grep -cvE '^-?[0-9]|^NA$|^NaN$|^na$|^nan$|^\.$|^$' || true)
        if [ "$bad_fields" -gt 0 ]; then
            errors="${errors}  WARN ${bad_fields} non-numeric field(s) at row 1000\n"
            n_numeric_fail=$((n_numeric_fail + 1))
        fi
    fi

    # Report
    if [ -z "$errors" ]; then
        n_ok=$((n_ok + 1))
        # Print progress every 100 files, or for failures
        if (( total % 100 == 0 )); then
            printf "  [%4d]  OK so far: %d\n" "$total" "$n_ok"
        fi
    else
        printf "  WARN  %-20s  %d cols\n" "$fname" "$ncols"
        printf "$errors" | tee -a "$REPORT"
    fi
done

# ── Summary ──────────────────────────────────────────────────────────────────

log ""
log "================================================================"
log "  Summary"
log "================================================================"
log "  Total files      : ${total}"
log "  Passed           : ${n_ok}"
log "  gzip failures    : ${n_gzip_fail}"
log "  Empty files      : ${n_empty}"
log "  Column mismatches: ${n_header_fail}"
log "  Numeric warns    : ${n_numeric_fail}"
log ""

# Show reference header
if [ -n "$ref_header" ]; then
    log "  Reference header (${ref_fname}):"
    log "  Columns: ${ref_ncols}"
    log "  $(echo "$ref_header" | tr '\t' '\n' | cat -n || true)"
fi

log ""
log "  Report saved to: ${REPORT}"
