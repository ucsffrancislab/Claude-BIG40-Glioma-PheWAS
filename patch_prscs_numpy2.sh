#!/usr/bin/env bash
# =============================================================================
# patch_prscs_numpy2.sh
#
# Patch a local PRS-CS install for NumPy 2.x compatibility.
#
# Fixes TypeError in gigrnd.py ("only 0-dimensional arrays can be converted
# to Python scalars") by rewriting the call site in mcmc_gtb.py to pass
# scalars instead of 1-element arrays.
#
# Idempotent — safe to re-run. Backs up the original as mcmc_gtb.py.bak.
#
# Usage:
#   bash patch_prscs_numpy2.sh [PRSCS_DIR]      # default: ~/.local/PRScs
# =============================================================================
set -eu

PRSCS_DIR="${1:-/sessions/operon-k-5bbbcdfd-548a-a657/.local/PRScs}"
FILE="$PRSCS_DIR/mcmc_gtb.py"

if [ ! -f "$FILE" ]; then
    echo "ERROR: $FILE not found"
    exit 1
fi

# The original line (as shipped by upstream):
OLD='psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)'
# Patched: index 2-D arrays with [jj,0] so we pass scalars to gigrnd
NEW='psi[jj] = gigrnd.gigrnd(a-0.5, float(2.0*delta[jj,0]), float(n*beta[jj,0]**2/sigma))'

if grep -qF "$NEW" "$FILE"; then
    echo "already patched: $FILE"
    exit 0
fi

if ! grep -qF "$OLD" "$FILE"; then
    echo "ERROR: expected line not found in $FILE"
    echo "       (upstream may have changed; please inspect manually)"
    exit 1
fi

cp "$FILE" "${FILE}.bak"
python3 - "$FILE" "$OLD" "$NEW" <<PYEOF
import sys
fp, old, new = sys.argv[1], sys.argv[2], sys.argv[3]
s = open(fp).read()
assert old in s
open(fp, "w").write(s.replace(old, new))
print(f"patched: {fp}")
print(f"backup:  {fp}.bak")
PYEOF
