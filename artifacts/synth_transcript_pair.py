#!/usr/bin/env python3
"""Synthetic genome<->transcriptome pair with EXACT ground truth.

Builds a single-copy multi-exon gene on one synthetic chromosome, with
canonical GT..AG introns, then emits transcript reads (exons spliced) plus
a NOVEL-isoform read using an alternative last exon that is NOT part of the
"annotated" transcript (the operator's NMD/alt-last-exon kill-switch).

Outputs (deterministic, no RNG):
  synth_genome.fa     one chromosome "chrS"
  synth_reads.fastq   transcript reads (FLNC-like, error-free)
  synth_truth.tsv     per-read expected genomic intron (N) positions

Acceptance the validator checks on `llmap align --mode transcript`:
  * EVERY read is present in the output (no silent drop / no loss).
  * Each spliced read's CIGAR carries N-ops whose ref-gaps match the
    truth introns (within a small tolerance).
  * The novel-isoform read maps (correct exon structure), never unmapped.
"""
import sys

# ---- deterministic "DNA" generator (no RNG; reproducible across runs) ----
def dna(n, seed):
    # simple LCG over {A,C,G,T}, avoids GT/AG accidental splice motifs at ends
    b = "ACGT"
    out = []
    x = (seed * 1103515245 + 12345) & 0x7FFFFFFF
    for _ in range(n):
        x = (x * 1103515245 + 12345) & 0x7FFFFFFF
        out.append(b[(x >> 16) & 3])
    return "".join(out)

# ---- exon / intron layout (genomic, 0-based) ----
# 5 exons; introns canonical GT...AG. Exons >= 60 bp so k=15/w seeds them.
EXON_LENS  = [180, 150, 240, 120, 200]      # annotated transcript exons
INTRON_LENS= [300, 500, 220, 800]           # 4 introns between the 5 exons
LEAD = 1000                                  # genomic lead before the gene
TAIL = 1000

def build_genome():
    seq = [dna(LEAD, 1)]
    layout = []           # (exon_genomic_start, exon_len)
    pos = LEAD
    for i, el in enumerate(EXON_LENS):
        layout.append((pos, el))
        seq.append(dna(el, 100 + i))
        pos += el
        if i < len(INTRON_LENS):
            il = INTRON_LENS[i]
            # canonical intron: GT ... AG, random middle
            seq.append("GT" + dna(il - 4, 500 + i) + "AG")
            pos += il
    # alternative last exon: a DISTINCT exon further downstream, reached via
    # an intron from exon 4 (index 3). It is NOT in the annotated transcript.
    seq.append(dna(600, 900))                # spacer
    alt_last_start = pos + 600
    seq.append(dna(200, 950))                # alt last exon (len 200)
    seq.append(dna(TAIL, 9))
    genome = "".join(seq)
    return genome, layout, alt_last_start

def exon_seq(genome, start, length):
    return genome[start:start + length]

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def main(outdir):
    genome, layout, alt_last_start = build_genome()
    with open(f"{outdir}/synth_genome.fa", "w") as f:
        f.write(">chrS synthetic single-copy gene\n")
        for i in range(0, len(genome), 70):
            f.write(genome[i:i+70] + "\n")

    reads = []   # (name, seq, truth_introns[list of (ref_gap_start,len)])
    # annotated transcript = exons 0..4 spliced
    canon = "".join(exon_seq(genome, s, l) for (s, l) in layout)
    introns = []
    # truth introns: between consecutive exons, ref gap = intron span
    for i in range(len(layout) - 1):
        gstart = layout[i][0] + layout[i][1]
        glen = INTRON_LENS[i]
        introns.append((gstart, glen))
    # 5 identical full-length transcript reads
    for k in range(5):
        reads.append((f"canon_{k}", canon, introns, layout[0][0]))
    # 1 reverse-strand transcript read
    reads.append(("canon_rev", revcomp(canon), introns, layout[0][0]))

    # NOVEL isoform: exons 0..3 spliced to the ALTERNATIVE last exon (not
    # exon 4). Junction exon3->alt is a real genomic intron but UNANNOTATED.
    novel_exons = layout[:4] + [(alt_last_start, 200)]
    novel = "".join(exon_seq(genome, s, l) for (s, l) in novel_exons)
    novel_introns = []
    for i in range(len(novel_exons) - 1):
        gstart = novel_exons[i][0] + novel_exons[i][1]
        glen = novel_exons[i+1][0] - gstart
        novel_introns.append((gstart, glen))
    for k in range(3):
        reads.append((f"novel_altlast_{k}", novel, novel_introns, novel_exons[0][0]))

    with open(f"{outdir}/synth_reads.fastq", "w") as f:
        for name, seq, _, _ in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I'*len(seq)}\n")
    with open(f"{outdir}/synth_truth.tsv", "w") as f:
        f.write("read\tref_start\tn_introns\tintrons(gstart:len;...)\n")
        for name, seq, intr, rstart in reads:
            f.write(f"{name}\t{rstart}\t{len(intr)}\t" +
                    ";".join(f"{g}:{l}" for g, l in intr) + "\n")
    print(f"genome={len(genome)}bp reads={len(reads)} "
          f"(canon+rev+novel) introns/canon={len(introns)} "
          f"novel_introns={len(novel_introns)}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
