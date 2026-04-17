#!/usr/bin/env bash
# =============================================================================
# patch_prscs_numpy2.sh
#
# Patch a local PRS-CS install for NumPy 2.x compatibility.
#
# Root cause: PRS-CS uses 2-D column-vector arrays (shape (p,1)) and relies on
# Python's built-in sum() + implicit scalar conversion to unbox 1-element
# ndarrays. NumPy 2.x made float(1-D-array) and "%d"/"%e" % 1-D-array into
# hard errors. We fix the four places where this leaks out.
#
# Fixes in mcmc_gtb.py:
#   L66: err = max(...)            -- sum() -> np.sum(), cast to float
#   L72: psi[jj] = gigrnd(...)     -- index 2-D arrays as [jj,0], cast to float
#   L77: phi = random.gamma(...)   -- sum(delta) -> float(np.sum(delta))
#   L110: write loop for beta_est  -- zip over beta_est.flatten() so each
#                                     element yielded is a numpy scalar
#   L120: write loop for psi_est   -- zip over psi_est.flatten() (symmetry;
#                                     matters only if --write_psi=TRUE)
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
    # Fix 1: err line
    (
        "        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))",
        "        err = float(max(n/2.0*(1.0-2.0*float(np.sum(beta*beta_mrg))+float(np.sum(quad))), n/2.0*float(np.sum(beta**2/psi))))",
    ),
    # Fix 2: gigrnd call site
    (
        "            psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)",
        "            psi[jj] = gigrnd.gigrnd(a-0.5, float(2.0*delta[jj,0]), float(n*beta[jj,0]**2/sigma))",
    ),
    # Fix 3: phi update
    (
        "            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))",
        "            phi = random.gamma(p*b+0.5, 1.0/(float(np.sum(delta))+w))",
    ),
    # Fix 4: write beta_est (default output path)
    (
        "            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_est):",
        "            for snp, bp, a1, a2, beta in zip(sst_dict['SNP'], sst_dict['BP'], sst_dict['A1'], sst_dict['A2'], beta_est.flatten()):",
    ),
    # Fix 5: write psi_est (only if write_psi=TRUE)
    (
        "            for snp, psi in zip(sst_dict['SNP'], psi_est):",
        "            for snp, psi in zip(sst_dict['SNP'], psi_est.flatten()):",
    ),
]

for old, new in subs:
    if old not in s:
        print(f"ERROR: pattern not found (unexpected upstream change?):", file=sys.stderr)
        print(f"  looking for: {old!r}", file=sys.stderr)
        sys.exit(1)
    s = s.replace(old, new, 1)

open(fp, "w").write(s)
print(f"applied 5 fixes to {fp}")
PYEOF

echo "patched: $FILE"
echo "backup:  ${FILE}.bak"
