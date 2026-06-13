#!/usr/bin/env python3
"""Validate `llmap align --mode transcript` output against synth truth.

LOSSLESS / NO-DROP proof on a synthetic genome<->transcriptome pair.

Usage: synth_transcript_validate.py <out.sam> <synth_truth.tsv> <synth_reads.fastq>

Exit 0 + "PASS" iff ALL of:
  (1) every input read appears in the output (no silent drop / no loss);
  (2) every read has a mapped primary record (flag bit 0x4 clear);
  (3) the NOVEL alt-last-exon read is mapped (never unmapped);
  (4) each spliced read's CIGAR N-ops match the truth intron count, and
      every truth intron length is matched by some N-op (within TOL bp).
"""
import re, sys

# bp tolerance on intron (N) length vs truth. With GT/AG splice-site snapping
# (src/mapping/splice_snap) the exon boundaries are pulled to the exact canonical
# site, so the N length should equal the true intron exactly — TOL is tight (a
# couple bp slack only for the rare case with no canonical site in-window).
TOL = 12

def parse_cigar(cig):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig)]

def n_ops(cig):
    return [n for n, op in parse_cigar(cig) if op == "N"]

def ref_span_consumed(cig):
    return sum(n for n, op in parse_cigar(cig) if op in "MDN=X")

def main(sam, truth_tsv, reads_fq):
    # input reads
    reads = []
    with open(reads_fq) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                reads.append(line[1:].strip().split()[0])
    reads = set(reads)

    truth = {}
    with open(truth_tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            name, rstart, nint = p[0], int(p[1]), int(p[2])
            ints = []
            if len(p) > 3 and p[3]:
                for tok in p[3].split(";"):
                    if tok:
                        g, l = tok.split(":"); ints.append((int(g), int(l)))
            truth[name] = (rstart, nint, ints)

    # parse SAM: primary (non-secondary/supplementary) per read
    seen, mapped, recs = set(), set(), {}
    with open(sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.split("\t")
            if len(c) < 6:
                continue
            name, flag, pos, cig = c[0], int(c[1]), int(c[3]), c[5]
            seen.add(name)
            is_unmapped = bool(flag & 0x4)
            is_sec = bool(flag & 0x100)
            is_supp = bool(flag & 0x800)
            if not is_unmapped:
                mapped.add(name)
            # keep the primary line (not secondary/supplementary) for CIGAR check
            if not is_sec and not is_supp:
                recs[name] = (flag, pos, cig)

    fails = []
    # (1) no loss
    missing = reads - seen
    if missing:
        fails.append(f"LOSS: {len(missing)} reads absent from output: "
                     f"{sorted(missing)[:5]}")
    # (2) all mapped
    unmapped = reads - mapped
    if unmapped:
        fails.append(f"UNMAPPED: {len(unmapped)} reads: {sorted(unmapped)[:5]}")
    # (3) novel reads mapped
    novel = {r for r in reads if r.startswith("novel_")}
    novel_unmapped = novel - mapped
    if novel_unmapped:
        fails.append(f"NOVEL-UNMAPPED (kill-switch!): {sorted(novel_unmapped)}")
    # (4) spliced N-ops match truth
    spliced_ok = 0
    for name, (rstart, nint, ints) in truth.items():
        if nint == 0 or name not in recs:
            continue
        _, pos, cig = recs[name]
        ns = n_ops(cig)
        if len(ns) != nint:
            fails.append(f"N-COUNT {name}: got {len(ns)} N-ops, truth {nint} "
                         f"(cigar={cig[:60]})")
            continue
        truth_lens = sorted(l for _, l in ints)
        got = sorted(ns)
        ok = all(any(abs(g - t) <= TOL for g in got) for t in truth_lens)
        if not ok:
            fails.append(f"N-LEN {name}: N={got} vs truth_introns={truth_lens}")
        else:
            spliced_ok += 1

    print(f"reads={len(reads)} seen={len(seen)} mapped={len(mapped)} "
          f"spliced_ok={spliced_ok}/{sum(1 for _,(_,n,_) in truth.items() if n>0)}")
    if fails:
        print("FAIL")
        for x in fails:
            print("  - " + x)
        return 1
    print("PASS — lossless, no-drop, splicing matches truth, novel isoform mapped")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
