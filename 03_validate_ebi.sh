#!/usr/bin/env bash
# =============================================================================
# 03_validate_ebi.sh
#
# Validates downloaded EBI GWAS summary statistics files (stats22k).
#
# Checks per file (fast — reads only the first ~1000 lines):
#   1. Not empty
#   2. Column headers match expected format (harmonised vs raw)
#   3. Numeric spot-check at row 1000
#
# Usage:
#   bash 03_validate_ebi.sh [DATA_DIR]
#
#   DATA_DIR  Base directory (default: data/big40)
# =============================================================================

set -eu

DATA_DIR="${1:-data/big40}"
STATS_DIR="${DATA_DIR}/stats22k"
REPORT="${DATA_DIR}/validation_ebi.txt"

# ── Setup ────────────────────────────────────────────────────────────────────

log() { echo "$*" | tee -a "$REPORT"; }

echo "" > "$REPORT"
log "================================================================"
log "  EBI Download Validation Report"
log "  $(date '+%Y-%m-%d %H:%M:%S')"
log "  Directory: ${STATS_DIR}"
log "================================================================"
log ""

n_harmonised=0
n_raw=0
n_col_fail=0
n_numeric_fail=0
n_empty=0
n_ok=0
total=0

# ── Validate ─────────────────────────────────────────────────────────────────

for f in "${STATS_DIR}"/*.tsv.gz; do
    [ -f "$f" ] || continue
    total=$((total + 1))

    fname=$(basename "$f")
    idp="${fname%%.*}"
    status="OK"
    note=""

    # Determine file type
    if [[ "$fname" == *.h.tsv.gz ]]; then
        ftype="harmonised"
        n_harmonised=$((n_harmonised + 1))
    elif [[ "$fname" == *.raw_hg19.tsv.gz ]]; then
        ftype="raw_hg19"
        n_raw=$((n_raw + 1))
    else
        ftype="unknown"
    fi

    # Read first 1000 lines once
    chunk=$(zcat "$f" 2>/dev/null | head -1000 || true)

    # 1. Empty check
    if [ -z "$chunk" ]; then
        n_empty=$((n_empty + 1))
        log "$(printf "  FAIL  %-30s  %-12s  empty file" "$fname" "$ftype")"
        continue
    fi

    header=$(echo "$chunk" | head -1)
    ncols=$(echo "$header" | awk -F'\t' '{print NF}')

    # 2. Column check — look for expected columns by file type
    missing=""
    if [ "$ftype" = "harmonised" ]; then
        for col in variant_id chromosome base_pair_location effect_allele other_allele beta standard_error p_value; do
            if ! echo "$header" | grep -qi "$col"; then
                missing="${missing} ${col}"
            fi
        done
    elif [ "$ftype" = "raw_hg19" ]; then
        for col in beta se; do
            if ! echo "$header" | grep -qi "$col"; then
                missing="${missing} ${col}"
            fi
        done
    fi

    if [ -n "$missing" ]; then
        status="WARN"
        note="missing:${missing}"
        n_col_fail=$((n_col_fail + 1))
    fi

    # 3. Numeric spot-check at row 1000
    spot=$(echo "$chunk" | tail -1)
    if [ -n "$spot" ]; then
        beta_col=$(echo "$header" | tr '\t' '\n' | grep -ni "^beta$" | head -1 | cut -d: -f1 || true)
        if [ -n "$beta_col" ]; then
            beta_val=$(echo "$spot" | cut -f"$beta_col")
            if ! echo "$beta_val" | grep -qE '^-?[0-9]|^NA$|^na$|^NaN$|^$'; then
                status="WARN"
                note="${note:+${note}; }beta not numeric: '${beta_val}'"
                n_numeric_fail=$((n_numeric_fail + 1))
            fi
        fi
    fi

    if [ "$status" = "OK" ]; then
        n_ok=$((n_ok + 1))
    fi

    printf "  %-6s %-30s  %-12s  %d cols  %s\n" "$status" "$fname" "$ftype" "$ncols" "$note"

done

# ── Summary ──────────────────────────────────────────────────────────────────

log ""
log "================================================================"
log "  Summary"
log "================================================================"
log "  Total files     : ${total}"
log "  Harmonised      : ${n_harmonised}"
log "  Raw (hg19)      : ${n_raw}"
log "  Passed          : ${n_ok}"
log ""
log "  Empty files     : ${n_empty}"
log "  Column issues   : ${n_col_fail}"
log "  Numeric warns   : ${n_numeric_fail}"
log ""

# Dump sample headers
sample_h=$(find "$STATS_DIR" -name '*.h.tsv.gz' 2>/dev/null | head -1)
if [ -n "$sample_h" ]; then
    log "  Harmonised header ($(basename "$sample_h")):"
    zcat "$sample_h" 2>/dev/null | head -1 | tr '\t' '\n' | cat -n | while read line; do
        log "    $line"
    done
fi

sample_raw=$(find "$STATS_DIR" -name '*.raw_hg19.tsv.gz' 2>/dev/null | head -1)
if [ -n "$sample_raw" ]; then
    log ""
    log "  Raw header ($(basename "$sample_raw")):"
    zcat "$sample_raw" 2>/dev/null | head -1 | tr '\t' '\n' | cat -n | while read line; do
        log "    $line"
    done
fi

log ""
log "  Report saved to: ${REPORT}"
