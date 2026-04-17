#!/usr/bin/env bash
# =============================================================================
# patch_prscs_numpy2.sh
#
# Patch a local PRS-CS install for NumPy 2.x compatibility.
#
# Root cause: PRS-CS uses Python's built-in sum() on 2-D numpy arrays,
# which iterates along axis 0 and returns 1-D arrays instead of scalars.
# This propagates through err -> sigma -> gigrnd call site. In NumPy 2.x,
# float(1-D-array) raises TypeError; older NumPy silently unboxed.
#
# Fixes three lines in mcmc_gtb.py:
#   L66: err = max(...)              -- replace sum() with np.sum(), cast to float
#   L72: psi[jj] = gigrnd.gigrnd(...) -- index 2-D arrays as [jj,0], cast to float
#   L77: phi = random.gamma(...)     -- replace sum() with np.sum()
#
# Idempotent and self-healing:
#   - On first run, backs up original as mcmc_gtb.py.bak (immutable thereafter).
#   - On every run, restores from .bak before patching, so partial prior
#     patches are cleanly replaced with the full current fix.
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

# Create immutable backup on first run
if [ ! -f "${FILE}.bak" ]; then
    cp "$FILE" "${FILE}.bak"
    echo "backup created: ${FILE}.bak"
fi

# Always restore from backup so we start from a known clean state
cp "${FILE}.bak" "$FILE"

python3 - "$FILE" <<'PYEOF'
import sys
fp = sys.argv[1]
s = open(fp).read()

subs = [
    # Fix 1: err line -- use np.sum, cast to scalar
    (
        "        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))",
        "        err = float(max(n/2.0*(1.0-2.0*float(np.sum(beta*beta_mrg))+float(np.sum(quad))), n/2.0*float(np.sum(beta**2/psi))))",
    ),
    # Fix 2: gigrnd call site -- 2-D index and float cast
    (
        "            psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)",
        "            psi[jj] = gigrnd.gigrnd(a-0.5, float(2.0*delta[jj,0]), float(n*beta[jj,0]**2/sigma))",
    ),
    # Fix 3: phi update -- np.sum instead of Python sum
    (
        "            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))",
        "            phi = random.gamma(p*b+0.5, 1.0/(float(np.sum(delta))+w))",
    ),
]

for old, new in subs:
    if old not in s:
        print(f"ERROR: pattern not found (unexpected upstream change?):", file=sys.stderr)
        print(f"  looking for: {old!r}", file=sys.stderr)
        sys.exit(1)
    s = s.replace(old, new, 1)

open(fp, "w").write(s)
print(f"applied 3 fixes to {fp}")
PYEOF

echo "patched: $FILE"
echo "backup:  ${FILE}.bak"
