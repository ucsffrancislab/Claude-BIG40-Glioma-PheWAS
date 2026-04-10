#!/usr/bin/env bash
# =============================================================================
# 03_validate_ebi.sh
#
# Validates downloaded EBI GWAS summary statistics files.
#
# Checks:
#   1. gzip integrity (can we decompress it?)
#   2. Column headers match expected format (harmonised vs raw)
#   3. Row count is in the expected ballpark (~10-17M SNPs)
#   4. Numeric columns are actually numeric (spot-check)
#   5. Summary report of all files
#
# Usage:
#   bash 03_validate_ebi.sh [DATA_DIR]
#
#   DATA_DIR  Base directory (default: data/big40)
# =============================================================================

set -euo pipefail

DATA_DIR="${1:-data/big40}"
STATS_DIR="${DATA_DIR}/stats22k"
REPORT="${DATA_DIR}/validation_report.txt"

# Expected columns in harmonised files (EBI GWAS-SSF format)
# These are the typical harmonised column names; we check a key subset.
EXPECTED_H_COLS="variant_id|rsid|chromosome|base_pair_location|effect_allele|other_allele|beta|standard_error|p_value"

# Expected columns in BIG40 raw files
EXPECTED_RAW_COLS="chr|rsid|pos|a1|a2|beta|se"

# ── Functions ────────────────────────────────────────────────────────────────

log() {
    echo "$*" | tee -a "$REPORT"
}

# ── Setup ────────────────────────────────────────────────────────────────────

echo "" > "$REPORT"
log "================================================================"
log "  EBI Download Validation Report"
log "  $(date '+%Y-%m-%d %H:%M:%S')"
log "  Directory: ${STATS_DIR}"
log "================================================================"
log ""

n_harmonised=0
n_raw=0
n_gzip_fail=0
n_col_fail=0
n_rowcount_warn=0
n_numeric_fail=0
n_ok=0
total=0

# ── Validate each file ──────────────────────────────────────────────────────

for f in "${STATS_DIR}"/*.tsv.gz; do
    [ -f "$f" ] || continue
    total=$((total + 1))

    fname=$(basename "$f")
    idp="${fname%%.*}"
    errors=""

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

    # 1. gzip integrity
    if ! gzip -t "$f" 2>/dev/null; then
        errors="${errors}  FAIL gzip integrity\n"
        n_gzip_fail=$((n_gzip_fail + 1))
        log "FAIL  ${fname}  gzip corrupt"
        continue
    fi

    # 2. Column headers
    header=$(zcat "$f" | head -1)

    if [ "$ftype" = "harmonised" ]; then
        missing_cols=""
        for col in variant_id chromosome base_pair_location effect_allele other_allele beta standard_error p_value; do
            if ! echo "$header" | grep -qi "$col"; then
                missing_cols="${missing_cols} ${col}"
            fi
        done
        if [ -n "$missing_cols" ]; then
            errors="${errors}  WARN missing harmonised columns:${missing_cols}\n"
            n_col_fail=$((n_col_fail + 1))
        fi
    elif [ "$ftype" = "raw_hg19" ]; then
        missing_cols=""
        # Raw BIG40 files may have: chr, rsid, pos, a1, a2, beta, se, neglog10p
        # or the EBI non-harmonised format — check what we actually got
        for col in beta se; do
            if ! echo "$header" | grep -qi "$col"; then
                missing_cols="${missing_cols} ${col}"
            fi
        done
        if [ -n "$missing_cols" ]; then
            errors="${errors}  WARN missing raw columns:${missing_cols}\n"
            n_col_fail=$((n_col_fail + 1))
        fi
    fi

    # 3. Row count (sample first 10M lines, extrapolate if needed)
    #    BIG40 GWAS has ~10-17M SNPs.  Flag if < 1M or > 25M.
    rowcount=$(zcat "$f" | wc -l)
    rowcount=$((rowcount - 1))  # subtract header

    if [ "$rowcount" -lt 1000000 ]; then
        errors="${errors}  WARN low row count: ${rowcount} (expected ~10-17M)\n"
        n_rowcount_warn=$((n_rowcount_warn + 1))
    elif [ "$rowcount" -gt 25000000 ]; then
        errors="${errors}  WARN high row count: ${rowcount} (expected ~10-17M)\n"
        n_rowcount_warn=$((n_rowcount_warn + 1))
    fi

    # 4. Numeric spot-check: grab line 1000, check that beta-like column is numeric
    spot=$(zcat "$f" | sed -n '1000p')
    if [ -n "$spot" ]; then
        # Try to extract a field that should be numeric (beta or similar)
        # Use the header to find the beta column index
        beta_col=$(echo "$header" | tr '\t' '\n' | grep -ni "^beta$" | head -1 | cut -d: -f1)
        if [ -n "$beta_col" ]; then
            beta_val=$(echo "$spot" | cut -f"$beta_col")
            # Check if it looks numeric (allow scientific notation, NA)
            if ! echo "$beta_val" | grep -qE '^-?[0-9]|^NA$|^na$|^NaN$|^$'; then
                errors="${errors}  WARN beta column not numeric at row 1000: '${beta_val}'\n"
                n_numeric_fail=$((n_numeric_fail + 1))
            fi
        fi
    fi

    # Report
    if [ -z "$errors" ]; then
        n_ok=$((n_ok + 1))
        printf "  OK    %-30s  %-12s  %'d rows\n" "$fname" "$ftype" "$rowcount"
    else
        printf "  WARN  %-30s  %-12s  %'d rows\n" "$fname" "$ftype" "$rowcount"
        printf "$errors" | tee -a "$REPORT"
    fi

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
log "  gzip failures   : ${n_gzip_fail}"
log "  Column issues   : ${n_col_fail}"
log "  Row count warns : ${n_rowcount_warn}"
log "  Numeric warns   : ${n_numeric_fail}"
log ""
log "  Report saved to : ${REPORT}"

# Also dump the header of one harmonised and one raw file for reference
log ""
log "================================================================"
log "  Sample Headers"
log "================================================================"

sample_h=$(find "$STATS_DIR" -name '*.h.tsv.gz' | head -1)
if [ -n "$sample_h" ]; then
    log ""
    log "  Harmonised ($(basename "$sample_h")):"
    log "  $(zcat "$sample_h" | head -1 | tr '\t' '\n' | cat -n)"
fi

sample_raw=$(find "$STATS_DIR" -name '*.raw_hg19.tsv.gz' | head -1)
if [ -n "$sample_raw" ]; then
    log ""
    log "  Raw hg19 ($(basename "$sample_raw")):"
    log "  $(zcat "$sample_raw" | head -1 | tr '\t' '\n' | cat -n)"
fi
