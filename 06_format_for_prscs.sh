#!/usr/bin/env bash
# =============================================================================
# 06_format_for_prscs.sh
#
# Converts BIG40 stats33k GWAS files into PRS-CS input format.
#
# BIG40 stats33k format (space-delimited):
#   chr rsid pos a1 a2 beta se pval(-log10)
#
#   a1 = reference allele
#   a2 = alternative allele = EFFECT allele
#
# PRS-CS expected format (tab-delimited, with header):
#   SNP  A1  A2  BETA  P
#
#   A1 = effect allele (BIG40 a2)
#   A2 = other allele  (BIG40 a1)
#   P  = p-value       (convert from -log10)
#
# Crash-safe: write-protection pattern. Rerunnable.
#
# Usage:
#   bash 06_format_for_prscs.sh [DATA_DIR] [NJOBS]
#
#   DATA_DIR  Base directory  (default: data/big40)
#   NJOBS     Parallel workers (default: 4)
# =============================================================================

set -eu

DATA_DIR="${1:-data/big40}"
NJOBS="${2:-4}"
IN_DIR="${DATA_DIR}/stats33k"
OUT_DIR="${DATA_DIR}/prscs_input"
LOG_FILE="${DATA_DIR}/format_prscs.log"

mkdir -p "$OUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Format conversion started.  NJOBS=${NJOBS}" \
    | tee -a "$LOG_FILE"

# ── Worker ───────────────────────────────────────────────────────────────────

WORKER=$(mktemp)
cat > "$WORKER" << 'WORKER_EOF'
#!/usr/bin/env bash
set -eu

infile="$1"
out_dir="$2"

fname=$(basename "$infile" .txt.gz)
outfile="${out_dir}/${fname}.txt"

# Write-protected = already complete
if [ -f "$outfile" ] && [ ! -w "$outfile" ]; then
    echo "SKIP  ${fname}"
    exit 0
fi

# Convert:
#   Input:  chr rsid pos a1 a2 beta se neglog10p  (space-delimited, no header)
#   Output: SNP A1   A2  BETA P                   (tab-delimited, with header)
#
#   SNP  = $2  (rsid)
#   A1   = $5  (a2, effect allele)
#   A2   = $4  (a1, other allele)
#   BETA = $6  (beta)
#   P    = 10^(-$8)  (convert neglog10p to p-value)
#
# Notes:
#   - Skip any rows where rsid is "." or empty
#   - Cap very small p-values: if neglog10p > 300, set P = 1e-300

{
    printf "SNP\tA1\tA2\tBETA\tP\n"
    zcat "$infile" | awk '
    {
        rsid = $2
        a1   = $5     # effect allele (BIG40 a2)
        a2   = $4     # other allele  (BIG40 a1)
        beta = $6
        nlp  = $8     # -log10(p)

        # Skip missing rsids
        if (rsid == "." || rsid == "" || rsid == "rsid") next

        # Convert -log10(p) to p-value, cap at 1e-300
        if (nlp + 0 > 300) {
            p = 1e-300
        } else if (nlp + 0 <= 0) {
            p = 1
        } else {
            p = 10 ^ (-nlp)
        }

        printf "%s\t%s\t%s\t%s\t%.6e\n", rsid, a1, a2, beta, p
    }
    '
} > "$outfile"

# Sanity check: file should have at least 1M lines
nlines=$(wc -l < "$outfile")
if [ "$nlines" -lt 1000000 ]; then
    echo "WARN  ${fname}  only ${nlines} lines"
    rm -f "$outfile"
    exit 1
fi

chmod a-w "$outfile"
echo "OK    ${fname}  ${nlines} lines"
WORKER_EOF
chmod +x "$WORKER"

# ── Run ──────────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Converting stats33k → PRS-CS format  (${NJOBS} workers)"
echo "============================================================"

ls "${IN_DIR}"/*.txt.gz \
    | xargs -P "$NJOBS" -I {} \
        bash "$WORKER" {} "$OUT_DIR"

rm -f "$WORKER"

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Summary"
echo "============================================================"

completed=$(find "$OUT_DIR" -name '*.txt' ! -writable 2>/dev/null | wc -l)
total=$(ls "${IN_DIR}"/*.txt.gz 2>/dev/null | wc -l)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converted: ${completed} / ${total}" \
    | tee -a "$LOG_FILE"

# Spot check: show header + first 3 lines of one output file
sample=$(find "$OUT_DIR" -name '*.txt' ! -writable 2>/dev/null | head -1)
if [ -n "$sample" ]; then
    echo ""
    echo "  Sample output ($(basename "$sample")):"
    head -4 "$sample" | column -t
fi
