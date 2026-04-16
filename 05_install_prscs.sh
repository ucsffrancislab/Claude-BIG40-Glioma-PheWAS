#!/usr/bin/env bash
# =============================================================================
# 05_install_prscs.sh
#
# Installs PRS-CS and downloads the 1000 Genomes EUR LD reference panel.
#
# PRS-CS: Polygenic prediction via Bayesian regression and continuous
#         shrinkage priors (Ge et al., 2019, Nature Communications)
#         https://github.com/getian107/PRScs
#
# Requirements:
#   - Python 3 with scipy and h5py
#   - ~4.2 GB disk for the LD reference panel
#
# Usage:
#   bash 05_install_prscs.sh [INSTALL_DIR]
#
#   INSTALL_DIR  Where to install PRS-CS (default: ~/software/PRScs)
# =============================================================================

set -eu

INSTALL_DIR="${1:-${HOME}/software/PRScs}"
LD_DIR="${INSTALL_DIR}/ld_ref"

echo "============================================================"
echo "  PRS-CS Installation"
echo "  Install dir: ${INSTALL_DIR}"
echo "============================================================"
echo ""

# ── Step 1: Clone PRS-CS ────────────────────────────────────────────────────

if [ -d "${INSTALL_DIR}/PRScs" ]; then
    echo "PRS-CS already cloned at ${INSTALL_DIR}/PRScs"
else
    echo "Cloning PRS-CS from GitHub ..."
    mkdir -p "${INSTALL_DIR}"
    cd "${INSTALL_DIR}"
    git clone https://github.com/getian107/PRScs.git
    echo "  Done."
fi

# ── Step 2: Check Python dependencies ───────────────────────────────────────

echo ""
echo "Checking Python dependencies ..."

missing=""
python3 -c "import scipy" 2>/dev/null || missing="${missing} scipy"
python3 -c "import h5py"  2>/dev/null || missing="${missing} h5py"

if [ -n "$missing" ]; then
    echo ""
    echo "  Missing Python packages:${missing}"
    echo ""
    echo "  Install them with one of:"
    echo "    pip install${missing}"
    echo "    conda install${missing}"
    echo ""
    echo "  Or if you use a module system:"
    echo "    module load python/3.x"
    echo "    pip install --user${missing}"
    echo ""
    echo "  After installing, rerun this script."
    exit 1
else
    echo "  scipy ... OK"
    echo "  h5py  ... OK"
fi

# ── Step 3: Download 1000 Genomes EUR LD reference panel ────────────────────
#
# The LD reference panels are HDF5 files, one per chromosome.
# 1000G EUR panel is ~4.2 GB compressed.
#
# Download URL from the PRS-CS GitHub:
#   https://www.dropbox.com/s/p9aqanhxhxcruzz/ldblk_1kg_eur.tar.gz
#
# Alternative: UK Biobank EUR panel (~7 GB)
#   https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz

echo ""
echo "============================================================"
echo "  LD Reference Panel (1000 Genomes EUR)"
echo "============================================================"

LD_TARBALL="${INSTALL_DIR}/ldblk_1kg_eur.tar.gz"
LD_EXTRACTED="${LD_DIR}/ldblk_1kg_eur"

if [ -d "$LD_EXTRACTED" ] && [ -f "${LD_EXTRACTED}/snpinfo_1kg_hm3" ]; then
    echo "LD reference already extracted at ${LD_EXTRACTED}"
else
    mkdir -p "${LD_DIR}"

    if [ -f "$LD_TARBALL" ]; then
        echo "Tarball already downloaded."
    else
        echo "Downloading 1000G EUR LD panel (~4.2 GB) ..."
        echo "  This will take a while."
        wget -O "$LD_TARBALL"             "https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=1"
        echo "  Download complete."
    fi

    echo "Extracting ..."
    tar -xzf "$LD_TARBALL" -C "${LD_DIR}"
    echo "  Extracted to ${LD_DIR}"

    # Verify
    n_chr=$(ls "${LD_EXTRACTED}"/ldblk_1kg_chr*.hdf5 2>/dev/null | wc -l)
    echo "  Found ${n_chr} chromosome LD files."
    if [ "$n_chr" -lt 22 ]; then
        echo "  WARNING: Expected 22 chromosome files!"
    fi
    if [ -f "${LD_EXTRACTED}/snpinfo_1kg_hm3" ]; then
        echo "  SNP info file: snpinfo_1kg_hm3 ... OK"
    else
        echo "  WARNING: snpinfo_1kg_hm3 not found!"
    fi
fi

# ── Step 4: Verify installation ─────────────────────────────────────────────

echo ""
echo "============================================================"
echo "  Verification"
echo "============================================================"

PRSCS_PY="${INSTALL_DIR}/PRScs/PRScs.py"

if [ -f "$PRSCS_PY" ]; then
    echo "  PRS-CS script : ${PRSCS_PY}"
else
    echo "  ERROR: PRScs.py not found at ${PRSCS_PY}"
    exit 1
fi

echo "  LD ref dir    : ${LD_EXTRACTED}"
echo ""

# Quick test: print help
python3 "${PRSCS_PY}" --help 2>&1 | head -5 || true

echo ""
echo "============================================================"
echo "  Installation complete!"
echo "============================================================"
echo ""
echo "  To use PRS-CS in your scripts:"
echo ""
echo "    PRSCS=${PRSCS_PY}"
echo "    LD_REF=${LD_EXTRACTED}"
echo ""
echo "    python3 \${PRSCS} \\"
echo "      --ref_dir=\${LD_REF} \\"
echo "      --bim_prefix=<plink_bim_prefix> \\"
echo "      --sst_file=<summary_stats.txt> \\"
echo "      --n_gwas=33224 \\"
echo "      --out_dir=<output_dir> \\"
echo "      --chrom=<1-22>"
echo ""
echo "  IMPORTANT for SLURM: limit threads to avoid hogging the node:"
echo ""
echo "    export MKL_NUM_THREADS=1"
echo "    export NUMEXPR_NUM_THREADS=1"
echo "    export OMP_NUM_THREADS=1"
echo ""
