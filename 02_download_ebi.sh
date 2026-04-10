#!/usr/bin/env bash
# =============================================================================
# 02_download_ebi.sh
#
# Downloads BIG40 22k discovery GWAS summary statistics from the EBI GWAS
# Catalog FTP server.  These are the same 3,935 brain imaging IDPs as
# stats33k but from the discovery cohort only (N ≤ 22,138).
#
# EBI accessions: GCST90002426 through GCST90006360
# Source: Smith et al. (2021) Nature Neuroscience 24(5):737-745
#
# The EBI FTP is a proper data server and should be much faster than the
# Oxford BIG40 web server used for stats33k.
#
# Crash-safe: same write-protection pattern as 01_download_big40.sh.
#
# Usage:
#   bash 02_download_ebi.sh [DATA_DIR]
#
#   DATA_DIR  Base directory for downloads (default: data/big40)
# =============================================================================

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

EBI_BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

# BIG40 22k discovery: GCST90002426 (IDP 0001) .. GCST90006360 (IDP 3935)
# Accession = 90002425 + IDP number
ACCESSION_OFFSET=90002425
IDP_FIRST=1
IDP_LAST=3935

DATA_DIR="${1:-data/big40}"
STATS_DIR="${DATA_DIR}/stats22k"
LOG_FILE="${DATA_DIR}/download_ebi.log"

# ── Helper functions ─────────────────────────────────────────────────────────

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

# safe_download <url> <output_path>
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

# ebi_bucket <accession_number>
#   Returns the EBI bucket directory name for a given accession number.
#   EBI groups accessions in ranges of 1000:
#     GCST90002001-GCST90003000, GCST90003001-GCST90004000, etc.
ebi_bucket() {
    local acc_num=$1
    local bucket_start=$(( ((acc_num - 1) / 1000) * 1000 + 1 ))
    local bucket_end=$(( bucket_start + 999 ))
    echo "GCST${bucket_start}-GCST${bucket_end}"
}

# ── Setup ────────────────────────────────────────────────────────────────────

mkdir -p "$STATS_DIR"
log "EBI download session started.  DATA_DIR=${DATA_DIR}"

# ── Download ─────────────────────────────────────────────────────────────────
#
# EBI GWAS Catalog FTP directory structure:
#   summary_statistics/
#     GCST90002001-GCST90003000/       <-- bucket (groups of 1000)
#       GCST90002426/                  <-- accession directory
#         harmonised/
#           GCST90002426.h.tsv.gz      <-- harmonised summary stats
#
# We download the harmonised file and rename it to {IDP}.h.tsv.gz for clarity.

echo ""
echo "============================================================"
echo "  EBI GWAS Catalog — BIG40 22k discovery (3,935 IDPs)"
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
    acc_num=$(( ACCESSION_OFFSET + n ))
    accession="GCST${acc_num}"
    bucket=$(ebi_bucket "$acc_num")
    outfile="${STATS_DIR}/${idp}.h.tsv.gz"

    # Already complete
    if [[ -f "$outfile" && ! -w "$outfile" ]]; then
        skipped=$((skipped + 1))
        printf "\r  [%4d / %d]  IDP %s  SKIP" "$i" "$total" "$idp"
        continue
    fi

    printf "\r  [%4d / %d]  IDP %s (%s)  downloading ...          " \
        "$i" "$total" "$idp" "$accession"

    url="${EBI_BASE}/${bucket}/${accession}/harmonised/${accession}.h.tsv.gz"

    if safe_download "$url" "$outfile"; then
        downloaded=$((downloaded + 1))
        printf "\r  [%4d / %d]  IDP %s (%s)  OK   (%s)              \n" \
            "$i" "$total" "$idp" "$accession" "$(du -h "$outfile" | cut -f1)"
    else
        rc=$?
        if [[ $rc -eq 22 ]]; then
            not_found=$((not_found + 1))
            printf "\r  [%4d / %d]  IDP %s (%s)  NOT FOUND             \n" \
                "$i" "$total" "$idp" "$accession"
        else
            failed=$((failed + 1))
            printf "\r  [%4d / %d]  IDP %s (%s)  FAIL (curl rc=%d)    \n" \
                "$i" "$total" "$idp" "$accession" "$rc"
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
log "Not found   : ${not_found}  (accession may not exist or different path)"
log "Failed      : ${failed}  (network errors — rerun to retry)"
echo ""

if [[ "$failed" -gt 0 ]]; then
    log "Some downloads failed. Rerun this script to retry."
    exit 1
else
    log "All downloads complete."
    exit 0
fi
