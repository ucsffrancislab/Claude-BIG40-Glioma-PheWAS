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
# Most accessions have harmonised files (~900 MB gzipped).  A small number
# (~4%) only have the raw TSV (~1.3 GB uncompressed); for those the script
# falls back to downloading the raw file and gzipping it locally.
#
# Crash-safe: same write-protection pattern as 01_download_big40.sh.
#
# Usage:
#   bash 02_download_ebi.sh [DATA_DIR] [NJOBS]
#
#   DATA_DIR  Base directory for downloads  (default: data/big40)
#   NJOBS     Parallel download workers     (default: 4)
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
NJOBS="${2:-4}"
STATS_DIR="${DATA_DIR}/stats22k"
LOG_FILE="${DATA_DIR}/download_ebi.log"

# ── Setup ────────────────────────────────────────────────────────────────────

mkdir -p "$STATS_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] EBI download session started.  DATA_DIR=${DATA_DIR}  NJOBS=${NJOBS}" \
    | tee -a "$LOG_FILE"

# ── Worker script ────────────────────────────────────────────────────────────

WORKER=$(mktemp)
cat > "$WORKER" << 'WORKER_EOF'
#!/usr/bin/env bash
set -euo pipefail

n="$1"
ebi_base="$2"
acc_offset="$3"
stats_dir="$4"
log_file="$5"

idp=$(printf "%04d" "$n")
acc_num=$(( acc_offset + n ))
accession="GCST${acc_num}"

# EBI bucket: groups of 1000
bucket_start=$(( ((acc_num - 1) / 1000) * 1000 + 1 ))
bucket_end=$(( bucket_start + 999 ))
bucket="GCST${bucket_start}-GCST${bucket_end}"

acc_dir="${ebi_base}/${bucket}/${accession}"

# Check for either harmonised or raw output
outfile_h="${stats_dir}/${idp}.h.tsv.gz"
outfile_raw="${stats_dir}/${idp}.raw_hg19.tsv.gz"

# Already complete (either format)
if [ -f "$outfile_h" ] && [ ! -w "$outfile_h" ]; then
    echo "SKIP  ${idp}"
    exit 0
fi
if [ -f "$outfile_raw" ] && [ ! -w "$outfile_raw" ]; then
    echo "SKIP  ${idp}"
    exit 0
fi

# Try 1: harmonised file (most accessions have this)
url_h="${acc_dir}/harmonised/${accession}.h.tsv.gz"

if curl -fSL \
        --retry 3 \
        --retry-delay 10 \
        --connect-timeout 30 \
        --max-time 3600 \
        -o "$outfile_h" \
        "$url_h" 2>>"$log_file"; then
    chmod a-w "$outfile_h"
    echo "OK    ${idp} (${accession})  $(du -h "$outfile_h" | cut -f1)"
    exit 0
fi
rm -f "$outfile_h"

# Try 2: raw TSV — download and gzip on the fly
url_raw="${acc_dir}/${accession}_buildGRCh37.tsv"

if curl -fSL \
        --retry 3 \
        --retry-delay 10 \
        --connect-timeout 30 \
        --max-time 7200 \
        "$url_raw" 2>>"$log_file" \
    | gzip > "$outfile_raw"; then

    # Check that we actually got data (gzip of empty = 20 bytes)
    if [ -s "$outfile_raw" ] && [ "$(stat -c%s "$outfile_raw" 2>/dev/null || stat -f%z "$outfile_raw")" -gt 100 ]; then
        chmod a-w "$outfile_raw"
        echo "OK    ${idp} (${accession})  $(du -h "$outfile_raw" | cut -f1)  [raw→gzip]"
        exit 0
    fi
fi
rm -f "$outfile_raw"

echo "FAIL  ${idp} (${accession})"
WORKER_EOF
chmod +x "$WORKER"

# ── Download ─────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  EBI GWAS Catalog — BIG40 22k discovery (3,935 IDPs)"
echo "  Workers: ${NJOBS}"
echo "============================================================"

seq "$IDP_FIRST" "$IDP_LAST" \
    | xargs -P "$NJOBS" -I {} \
        bash "$WORKER" {} "$EBI_BASE" "$ACCESSION_OFFSET" "$STATS_DIR" "$LOG_FILE"

rm -f "$WORKER"

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Download Summary"
echo "============================================================"

completed=$(find "$STATS_DIR" -name '*.tsv.gz' ! -writable 2>/dev/null | wc -l)
total=$(( IDP_LAST - IDP_FIRST + 1 ))

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Files complete: ${completed} / ${total}" \
    | tee -a "$LOG_FILE"

if [ "$completed" -lt "$total" ]; then
    echo "Some files still missing. Rerun to retry." | tee -a "$LOG_FILE"
fi
