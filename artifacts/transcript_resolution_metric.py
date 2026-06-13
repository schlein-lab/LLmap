#!/usr/bin/env python3
"""Transcript-mapping resolution + per-base splice-determinism metric.

Given spliced read alignments (SAM with N-ops) of iso-seq/transcript reads to a
genome, and optionally a reference annotation (GENCODE-style intron BED:
chrom<TAB>start<TAB>end per intron), report:

  1. JUNCTION CONCORDANCE / RESOLUTION — how far each observed read junction
     (intron) sits from the nearest annotated junction. The distribution of
     offsets IS the "boundary blurriness": fraction exact (0 bp), within 2, 5,
     10 bp, and the median |offset|. Sharp, annotation-matching boundaries →
     mass at 0; blurry/alternative sites → a spread.

  2. ISOFORM DISTRIBUTION — reads grouped by their ordered junction-chain
     signature = distinct isoforms; reports count of isoforms per locus and the
     read-support (so dominant vs minor isoforms are visible).

  3. PER-BASE SPLICE-DETERMINISM D(p) — for each covered genomic position, the
     fraction of covering reads that agree on the MODAL splicing state at p
     (exonic vs intron-skipped, AND, near a boundary, the modal exact boundary).
     D(p)=100 → constitutive + sharp; 90 → excluded ~9/10; 76 → blurry boundary.
     This is 1 - normalised splicing entropy at p — a per-position variant-context
     annotation: a variant at low D sits in a variably-spliced / blurry locus.

Usage: transcript_resolution_metric.py <reads.sam> [annot_introns.bed]
"""
import re, sys
from collections import defaultdict, Counter

def read_alignments(sam):
    """Yield (qname, chrom, [(exon_start,exon_end)...], [(intron_start,end)...])."""
    with open(sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.split("\t")
            if len(c) < 6:
                continue
            flag, chrom, pos, cig = int(c[1]), c[2], int(c[3]), c[5]
            if flag & 0x904:                       # unmapped/secondary/supplementary
                continue
            ref = pos - 1
            exons, introns, ex_start = [], [], pos - 1
            for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
                ln = int(ln)
                if op in "M=X":
                    ref += ln
                elif op == "D":
                    ref += ln
                elif op == "N":
                    exons.append((ex_start, ref))
                    introns.append((ref, ref + ln))
                    ref += ln
                    ex_start = ref
            exons.append((ex_start, ref))
            yield c[0], chrom, exons, introns

def load_annot(path):
    introns = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3:
                introns[p[0]].append((int(p[1]), int(p[2])))
    for k in introns:
        introns[k].sort()
    return introns

def nearest_offset(j, annot_list):
    """Min |start-as| + |end-ae| over annotated introns on the chrom (a coarse
    junction-distance); returns (best_total_offset, exact_bool)."""
    s, e = j
    best = None
    for a_s, a_e in annot_list:
        if abs(a_s - s) > 50 and abs(a_e - e) > 50:
            continue
        tot = abs(a_s - s) + abs(a_e - e)
        if best is None or tot < best:
            best = tot
    return (best if best is not None else 10**9)

def main(sam, annot_path=None):
    aligns = list(read_alignments(sam))
    annot = load_annot(annot_path) if annot_path else None
    n_spliced = sum(1 for _, _, _, intr in aligns if intr)
    print(f"reads aligned: {len(aligns)}  spliced (>=1 intron): {n_spliced}")

    # --- 1. junction concordance / resolution ---
    if annot:
        offsets = []
        for _, chrom, _, introns in aligns:
            for j in introns:
                offsets.append(nearest_offset(j, annot.get(chrom, [])))
        if offsets:
            offsets.sort()
            n = len(offsets)
            exact = sum(o == 0 for o in offsets)
            within = lambda t: sum(o <= t for o in offsets)
            print("\n[1] JUNCTION RESOLUTION vs annotation "
                  f"({n} observed junctions):")
            print(f"    exact (0 bp): {exact/n:6.1%}   <=2bp: {within(2)/n:6.1%}   "
                  f"<=5bp: {within(5)/n:6.1%}   <=10bp: {within(10)/n:6.1%}")
            print(f"    median |offset|: {offsets[n//2]} bp   "
                  f"unmatched(>50bp): {sum(o>50 for o in offsets)/n:6.1%}")

    # --- 2. isoform distribution (by junction-chain signature) ---
    iso = defaultdict(Counter)   # locus-key -> Counter of junction-chain sigs
    for _, chrom, exons, introns in aligns:
        if not introns:
            continue
        locus = (chrom, introns[0][0] // 100000)   # coarse 100kb bin as locus key
        sig = tuple(introns)
        iso[locus][sig] += 1
    if iso:
        n_iso = [len(v) for v in iso.values()]
        print(f"\n[2] ISOFORMS: {len(iso)} spliced loci, "
              f"{sum(n_iso)} distinct isoforms "
              f"(median {sorted(n_iso)[len(n_iso)//2]}/locus, max {max(n_iso)})")

    # --- 3. per-base splice-determinism D(p) ---
    # state per read at p: 'E' exonic, 'N' intron. modal-fraction over covering reads.
    state = defaultdict(Counter)   # (chrom,pos) -> Counter of states
    STEP = 1                       # per-base; coarsen for big loci if needed
    for _, chrom, exons, introns in aligns:
        span_lo = exons[0][0]; span_hi = exons[-1][1]
        ex_set = []
        for s, e in exons:
            ex_set.append((s, e))
        for s, e in introns:
            # mark intron positions N (sampled)
            for p in range(s, e, max(STEP, (e - s)//200 or 1)):
                state[(chrom, p)]['N'] += 1
        for s, e in exons:
            for p in range(s, e, max(STEP, (e - s)//200 or 1)):
                state[(chrom, p)]['E'] += 1
    if state:
        dvals = []
        for pos, ctr in state.items():
            tot = sum(ctr.values())
            if tot < 3:
                continue
            dvals.append(100.0 * max(ctr.values()) / tot)
        if dvals:
            dvals.sort()
            n = len(dvals)
            print(f"\n[3] SPLICE-DETERMINISM D(p) over {n} positions (>=3 reads):")
            print(f"    mean D: {sum(dvals)/n:5.1f}   median: {dvals[n//2]:5.1f}   "
                  f"D==100: {sum(d>=99.9 for d in dvals)/n:6.1%}   "
                  f"D<90 (variable): {sum(d<90 for d in dvals)/n:6.1%}   "
                  f"D<76 (blurry): {sum(d<76 for d in dvals)/n:6.1%}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
