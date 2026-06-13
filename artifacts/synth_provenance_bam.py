#!/usr/bin/env python3
"""Synthetic provenance SAM with EXACT ground truth.

Validates the `llmap provenance-spectrum` BAM-derive path + Σ-invariant
deterministically (no RNG), the way synth_transcript_pair validated splicing.
Each read is planted with a known provenance class via the cheap BAM-derivable
signals the derive path reads:

  host : flag 0, MAPQ 60, plain CIGAR, no tags         -> host
  dup  : flag has 0x400 (PCR/optical duplicate)         -> dup
  chim : soft-clipped primary + SA:Z tag (split read)   -> chim

Catalog-dependent classes (para/numt/pseudo/exo/mei) are NOT exercised here —
they need PSV/NUMT/TE catalogs (Block 3) and would simply resolve to host from a
plain BAM. So the expected spectrum is over {host, dup, chim} and MUST sum to N.

Outputs:
  synth_provenance.sam      the tagged reads
  synth_provenance_truth.tsv  class -> expected count

Acceptance the validator checks on `llmap provenance-spectrum`:
  * Σ(all partition classes) == N_reads          (lossless, no read dropped)
  * count(host)==N_HOST, count(dup)==N_DUP, count(chim)==N_CHIM
"""
import sys

N_HOST = 70
N_DUP = 15
N_CHIM = 15
REF = "chrS"
REF_LEN = 100000

def header():
    return [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:{REF}\tLN:{REF_LEN}",
        "@PG\tID:synth_provenance\tPN:synth_provenance\tVN:1.0",
    ]

def host_read(i):
    # clean, unique, high MAPQ — host
    return f"host_{i}\t0\t{REF}\t{1000 + i}\t60\t100M\t*\t0\t0\t*\t*"

def dup_read(i):
    # flag 0x400 (1024) = PCR/optical duplicate
    return f"dup_{i}\t1024\t{REF}\t{2000 + i}\t60\t100M\t*\t0\t0\t*\t*"

def chim_read(i):
    # soft-clipped primary + SA tag = split/chimeric read
    pos = 3000 + i
    sa = f"SA:Z:{REF},{pos + 5000},+,40M60S,60,0;"
    return f"chim_{i}\t0\t{REF}\t{pos}\t60\t60M40S\t*\t0\t0\t*\t*\t{sa}"

def main(out_sam, out_truth):
    lines = header()
    for i in range(N_HOST):
        lines.append(host_read(i))
    for i in range(N_DUP):
        lines.append(dup_read(i))
    for i in range(N_CHIM):
        lines.append(chim_read(i))
    with open(out_sam, "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(out_truth, "w") as f:
        f.write("class\texpected_count\n")
        f.write(f"host\t{N_HOST}\n")
        f.write(f"dup\t{N_DUP}\n")
        f.write(f"chim\t{N_CHIM}\n")
        f.write(f"_total\t{N_HOST + N_DUP + N_CHIM}\n")
    print(f"wrote {out_sam} ({N_HOST + N_DUP + N_CHIM} reads) + {out_truth}")

if __name__ == "__main__":
    out_sam = sys.argv[1] if len(sys.argv) > 1 else "synth_provenance.sam"
    out_truth = sys.argv[2] if len(sys.argv) > 2 else "synth_provenance_truth.tsv"
    main(out_sam, out_truth)
