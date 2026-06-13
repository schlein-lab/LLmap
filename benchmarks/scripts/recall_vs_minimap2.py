#!/usr/bin/env python3
"""Recall + junction concordance: llmap vs minimap2 on the SAME spliced reads.

Measures the recall gap the Prio-1 fix targets (llmap spliced 5 vs minimap2 18)
so the fix can be quantified instead of guessed. Fair by construction — it
compares on the INTERSECTION of read IDs both tools processed (no easier-subset
artefact), per the shared-reads recommendation.

Usage:
  recall_vs_minimap2.py <llmap.sam> <minimap2.sam> [--tol BP] [--llmap-sec SEC] [--mm2-sec SEC]

Reports (to stdout, parseable):
  mapped(llmap/mm2), spliced(llmap/mm2),
  SPLICE-RECALL = |reads spliced by BOTH| / |reads spliced by mm2|   ← the headline
  junction concordance (Jaccard, donor/acceptor ±TOL on shared spliced reads)
  speed (reads/sec) if --*-sec timings given (the second benchmark axis)
"""
import re
import sys


def parse_args(argv):
    pos, opt = [], {"tol": 10, "llmap_sec": None, "mm2_sec": None}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--tol":
            opt["tol"] = int(argv[i + 1]); i += 2
        elif a == "--llmap-sec":
            opt["llmap_sec"] = float(argv[i + 1]); i += 2
        elif a == "--mm2-sec":
            opt["mm2_sec"] = float(argv[i + 1]); i += 2
        else:
            pos.append(a); i += 1
    return pos, opt


def n_ops_positions(pos1, cigar):
    """Genomic (donor, acceptor) of each N op, walking the CIGAR from POS (1-based)."""
    juncs = []
    ref = pos1
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op == "N":
            juncs.append((ref, ref + n))
            ref += n
        elif op in "MD=X":
            ref += n
    return juncs


def parse(sam):
    """read_id -> {mapped, spliced, juncs:set, primary_only}."""
    reads = {}
    with open(sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.split("\t")
            if len(c) < 6:
                continue
            name, flag, pos, cig = c[0], int(c[1]), int(c[3]), c[5]
            if flag & 0x100:  # skip secondary (count primary + supplementary)
                continue
            r = reads.setdefault(name, {"mapped": False, "spliced": False, "juncs": set()})
            if not (flag & 0x4):
                r["mapped"] = True
            if "N" in cig and pos > 0:
                r["spliced"] = True
                for d, a in n_ops_positions(pos, cig):
                    r["juncs"].add((d, a))
    return reads


def jaccard(a, b, tol):
    if not a and not b:
        return 1.0
    matched = sum(1 for x in a if any(abs(x[0] - y[0]) <= tol and abs(x[1] - y[1]) <= tol for y in b))
    union = len(a) + len(b) - matched
    return matched / union if union else 0.0


def main(llmap_sam, mm2_sam, opt):
    L, M = parse(llmap_sam), parse(mm2_sam)
    shared = set(L) & set(M)  # fair: only reads both tools processed

    l_map = sum(1 for r in L.values() if r["mapped"])
    m_map = sum(1 for r in M.values() if r["mapped"])
    l_spl = sum(1 for r in L.values() if r["spliced"])
    m_spl = sum(1 for r in M.values() if r["spliced"])

    mm2_spliced = {r for r in shared if M[r]["spliced"]}
    both_spliced = {r for r in mm2_spliced if L.get(r, {}).get("spliced")}
    recall = len(both_spliced) / len(mm2_spliced) if mm2_spliced else 0.0

    jac = [jaccard(L[r]["juncs"], M[r]["juncs"], opt["tol"]) for r in both_spliced]
    mean_jac = sum(jac) / len(jac) if jac else 0.0

    print(f"reads_shared={len(shared)}")
    print(f"mapped llmap={l_map} mm2={m_map}")
    print(f"spliced llmap={l_spl} mm2={m_spl}")
    print(f"SPLICE_RECALL={recall:.4f}  (both_spliced={len(both_spliced)} / mm2_spliced={len(mm2_spliced)})")
    print(f"junction_jaccard_mean={mean_jac:.4f} (tol={opt['tol']}bp, over shared spliced reads)")
    if opt["llmap_sec"] and opt["mm2_sec"]:
        ln, mn = len(L) or 1, len(M) or 1
        print(f"speed reads/sec llmap={ln/opt['llmap_sec']:.1f} mm2={mn/opt['mm2_sec']:.1f} "
              f"(slowdown={opt['llmap_sec']*mn/(opt['mm2_sec']*ln):.0f}x)")
    return 0


if __name__ == "__main__":
    pos, opt = parse_args(sys.argv[1:])
    if len(pos) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(pos[0], pos[1], opt))
