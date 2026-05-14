#!/usr/bin/env python3
"""Pairwise concordance analysis between two BAMs of the same input.

For every read present in both BAMs, classify the agreement:

    same_position       both mapped, same chrom, |pos_a - pos_b| <= tolerance
    same_chrom_diff_pos same chrom, distance > tolerance
    diff_chrom          different chrom
    oneA                mapped by A only
    oneB                mapped by B only
    both_unmapped       both unmapped

Writes:
    concordance.json    counts + percentages
    divergence_examples.tsv  up to --max-examples rows of disagreements
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import pysam  # type: ignore
except ImportError:
    print("ERROR: pysam not installed", file=sys.stderr)
    raise


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-a", required=True)
    ap.add_argument("--bam-b", required=True)
    ap.add_argument("--name-a", required=True)
    ap.add_argument("--name-b", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tolerance-bp", type=int, default=10)
    ap.add_argument("--max-examples", type=int, default=100)
    return ap.parse_args()


def index_primaries(bam_path: str) -> dict[str, tuple[str | None, int]]:
    out: dict[str, tuple[str | None, int]] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary:
                continue
            if r.is_unmapped:
                out[r.query_name] = (None, -1)
            else:
                chrom = bam.get_reference_name(r.reference_id)
                out[r.query_name] = (chrom, r.reference_start)
    return out


def classify(a: tuple[str | None, int], b: tuple[str | None, int], tol: int) -> str:
    chrom_a, pos_a = a
    chrom_b, pos_b = b
    if chrom_a is None and chrom_b is None:
        return "both_unmapped"
    if chrom_a is None:
        return "oneB"
    if chrom_b is None:
        return "oneA"
    if chrom_a != chrom_b:
        return "diff_chrom"
    if abs(pos_a - pos_b) <= tol:
        return "same_position"
    return "same_chrom_diff_pos"


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    a = index_primaries(args.bam_a)
    b = index_primaries(args.bam_b)

    counts = {
        "same_position": 0,
        "same_chrom_diff_pos": 0,
        "diff_chrom": 0,
        "oneA": 0,
        "oneB": 0,
        "both_unmapped": 0,
    }
    examples: list[tuple[str, str, tuple, tuple]] = []

    all_reads = set(a.keys()) | set(b.keys())
    for r in all_reads:
        ra = a.get(r, (None, -1))
        rb = b.get(r, (None, -1))
        cls = classify(ra, rb, args.tolerance_bp)
        counts[cls] += 1
        if cls != "same_position" and cls != "both_unmapped" and len(examples) < args.max_examples:
            examples.append((r, cls, ra, rb))

    total = sum(counts.values()) or 1
    summary = {
        "name_a": args.name_a,
        "name_b": args.name_b,
        "tolerance_bp": args.tolerance_bp,
        "counts": counts,
        "percentages": {k: round(100.0 * v / total, 3) for k, v in counts.items()},
        "total_reads": total,
    }
    (out_dir / "concordance.json").write_text(json.dumps(summary, indent=2))

    with open(out_dir / "divergence_examples.tsv", "w") as f:
        f.write("read_id\tclass\t%s_chrom\t%s_pos\t%s_chrom\t%s_pos\n" %
                (args.name_a, args.name_a, args.name_b, args.name_b))
        for r, cls, ra, rb in examples:
            f.write(f"{r}\t{cls}\t{ra[0] or '*'}\t{ra[1]}\t{rb[0] or '*'}\t{rb[1]}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
