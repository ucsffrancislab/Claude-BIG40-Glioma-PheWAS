#!/usr/bin/env bash
# =============================================================================
# patch_prscs_numpy2.sh
#
# Patch a local PRS-CS install for NumPy 2.x compatibility.
#
# Two fixes in mcmc_gtb.py:
#   (1) gigrnd call site: index 2-D arrays as [jj,0] and cast to float
#       (fixes "float(b)" error in gigrnd.py line 38)
#   (2) sigma: cast to float after computation, because the preceding
#       err computation uses Python's built-in sum() on 2-D arrays which
#       returns 1-D results, making sigma a (1,) array that propagates
#       through n*beta[jj,0]**2/sigma.
#
# Idempotent. Backs up the original as mcmc_gtb.py.bak (first run only;
# re-running does not overwrite the .bak).
#
# Usage:
#   bash patch_prscs_numpy2.sh [PRSCS_DIR]
#   default PRSCS_DIR: /c4/home/gwendt/.local/PRScs
# =============================================================================
set -eu

PRSCS_DIR="${1:-/c4/home/gwendt/.local/PRScs}"
FILE="$PRSCS_DIR/mcmc_gtb.py"

if [ ! -f "$FILE" ]; then
    echo "ERROR: $FILE not found"
    exit 1
fi

# Backup once (preserve the truly-original)
if [ ! -f "${FILE}.bak" ]; then
    cp "$FILE" "${FILE}.bak"
    echo "backup created: ${FILE}.bak"
fi

# --- Fix 1: gigrnd call site ---
OLD1='psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)'
NEW1='psi[jj] = gigrnd.gigrnd(a-0.5, float(2.0*delta[jj,0]), float(n*beta[jj,0]**2/sigma))'

# --- Fix 2: sigma scalar cast ---
OLD2='sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)'
NEW2='sigma = float(1.0/random.gamma((n+p)/2.0, 1.0/err))'

python3 - "$FILE" <<'PYEOF'
import sys
fp = sys.argv[1]
s = open(fp).read()

subs = [
    ("psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)",
     "psi[jj] = gigrnd.gigrnd(a-0.5, float(2.0*delta[jj,0]), float(n*beta[jj,0]**2/sigma))"),
    ("sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)",
     "sigma = float(1.0/random.gamma((n+p)/2.0, 1.0/err))"),
]

changes = 0
already = 0
for old, new in subs:
    if new in s:
        already += 1
        continue
    if old not in s:
        print(f"ERROR: expected pattern not found: {old!r}", file=sys.stderr)
        sys.exit(1)
    s = s.replace(old, new)
    changes += 1

open(fp, "w").write(s)
print(f"applied {changes} new change(s); {already} already in place")
PYEOF

echo "patched: $FILE"
