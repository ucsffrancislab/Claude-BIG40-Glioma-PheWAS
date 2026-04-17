#!/usr/bin/env python3
"""
07c_remap_bim_rsids.py

Rewrite a PLINK .bim file's variant-ID column (col 2) from chr:pos-style
IDs (as emitted by Michigan Imputation Server) to rsIDs, using PRS-CS's
HM3 SNP-info file as the (chr, bp) -> rsID lookup.

Needed because PRS-CS matches SNPs by rsID (its LD reference is HM3-keyed),
but imputed-VCF BIMs typically have variant IDs like '1:15585'.

Allele match: we require {bim_a1, bim_a2} == {hm3_a1, hm3_a2} (set equality
-- covers the two possible strand/REF-ALT orientations; PRS-CS handles the
flip internally). Variants with no match OR an allele mismatch keep their
original ID; PRS-CS simply won't use them (they're not in its LD ref anyway).

Usage:
    python 07c_remap_bim_rsids.py <input.bim> <output.bim> <snpinfo_1kg_hm3>
"""
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <input.bim> <output.bim> <snpinfo_1kg_hm3>")
    inp, outp, snpinfo = map(Path, sys.argv[1:4])

    print(f"[1/3] loading HM3 SNP info: {snpinfo}", flush=True)
    lookup = {}  # (chr_str, bp_str) -> (rsid, a1, a2)
    with open(snpinfo) as f:
        header = f.readline().split()
        ci, si, bi, a1i, a2i = (header.index(k) for k in ("CHR", "SNP", "BP", "A1", "A2"))
        for line in f:
            ll = line.split()
            lookup[(ll[ci], ll[bi])] = (ll[si], ll[a1i].upper(), ll[a2i].upper())
    print(f"       {len(lookup):,} HM3 variants loaded", flush=True)

    print(f"[2/3] remapping {inp} -> {outp}", flush=True)
    total = matched = allele_mismatch = 0
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(inp) as fi, open(outp, "w") as fo:
        for line in fi:
            # PLINK .bim: chr  snp  cm  bp  a1  a2   (tab OR whitespace)
            ll = line.split()
            if len(ll) < 6:
                fo.write(line)
                continue
            chrom, bp = ll[0], ll[3]
            a1, a2 = ll[4].upper(), ll[5].upper()
            total += 1
            hit = lookup.get((chrom, bp))
            if hit is not None:
                rsid, ha1, ha2 = hit
                if {a1, a2} == {ha1, ha2}:
                    ll[1] = rsid
                    matched += 1
                else:
                    allele_mismatch += 1
            fo.write("\t".join(ll) + "\n")

    print("[3/3] done", flush=True)
    print(f"       total BIM variants:       {total:,}", flush=True)
    pct = 100.0 * matched / total if total else 0.0
    print(f"       matched to HM3 rsID:      {matched:,}  ({pct:.2f}%)", flush=True)
    print(f"       (chr,pos) hit, allele differs: {allele_mismatch:,}", flush=True)


if __name__ == "__main__":
    main()
