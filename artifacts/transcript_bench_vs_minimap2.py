#!/usr/bin/env python3
"""Benchmark llmap --mode transcript vs minimap2 -ax splice on real FLNC reads.

Compares the two spliced aligners on the SAME reads + genome reference and
reports, per tool: mapping rate, #spliced reads (>=1 N op), total introns, and
the cross-tool junction concordance (Jaccard of the genomic intron set) +
per-read primary-position agreement. minimap2 -ax splice is the field-standard
reference; this measures how close llmap's spliced mapping is to it.

Usage: transcript_bench_vs_minimap2.py <llmap.sam> <minimap2.sam>
"""
import re, sys
from collections import defaultdict

def parse(path):
    reads = {}            # qname -> (rname, pos, cigar, flag) for primary
    introns = set()       # (rname, intron_start, intron_end) genomic, 0-based
    n_reads = n_mapped = n_spliced = 0
    seen = set()
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            qname, flag, rname, pos, cigar = c[0], int(c[1]), c[2], int(c[3]), c[5]
            if flag & 0x100 or flag & 0x800:   # skip secondary/supplementary
                continue
            if qname not in seen:
                seen.add(qname); n_reads += 1
            if flag & 0x4:                      # unmapped
                continue
            n_mapped += 1
            ops = [(int(n), o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
            has_n = any(o == "N" for _, o in ops)
            if has_n:
                n_spliced += 1
            # walk ref coords to collect intron (N) spans
            ref = pos - 1
            for ln, o in ops:
                if o in "MD=X":
                    ref += ln
                elif o == "N":
                    introns.add((rname, ref, ref + ln))
                    ref += ln
            reads[qname] = (rname, pos)
    return dict(n_reads=n_reads, n_mapped=n_mapped, n_spliced=n_spliced,
                introns=introns, reads=reads)

def jaccard(a, b, tol=10):
    # match introns within tol bp on both ends
    b_by_chrom = defaultdict(list)
    for r, s, e in b:
        b_by_chrom[r].append((s, e))
    inter = 0
    for r, s, e in a:
        if any(abs(s - bs) <= tol and abs(e - be) <= tol
               for bs, be in b_by_chrom.get(r, [])):
            inter += 1
    union = len(a) + len(b) - inter
    return inter, union, (inter / union if union else 1.0)

def main(llmap_sam, mm2_sam):
    L = parse(llmap_sam); M = parse(mm2_sam)
    print(f"{'metric':<26}{'llmap':>12}{'minimap2':>12}")
    print(f"{'reads':<26}{L['n_reads']:>12}{M['n_reads']:>12}")
    print(f"{'mapped':<26}{L['n_mapped']:>12}{M['n_mapped']:>12}")
    print(f"{'  mapping rate':<26}{L['n_mapped']/max(1,L['n_reads']):>12.3f}"
          f"{M['n_mapped']/max(1,M['n_reads']):>12.3f}")
    print(f"{'spliced reads (>=1 N)':<26}{L['n_spliced']:>12}{M['n_spliced']:>12}")
    print(f"{'distinct introns':<26}{len(L['introns']):>12}{len(M['introns']):>12}")
    inter, union, jac = jaccard(L['introns'], M['introns'])
    print(f"\nIntron concordance (±10bp): {inter} shared / {union} union "
          f"= Jaccard {jac:.3f}")
    # per-read primary placement agreement (same chrom, POS within 20bp)
    common = set(L['reads']) & set(M['reads'])
    agree = sum(1 for q in common
                if L['reads'][q][0] == M['reads'][q][0]
                and abs(L['reads'][q][1] - M['reads'][q][1]) <= 20)
    print(f"Primary placement agreement: {agree}/{len(common)} reads "
          f"(same chrom, POS ±20bp)")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
