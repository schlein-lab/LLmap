#!/usr/bin/env python3
"""Compute per-run mapping metrics from a BAM file.

Inputs:
    --bam        sorted+indexed BAM
    --truth      optional TSV of ground truth: read_id<TAB>chrom<TAB>pos
    --out-dir    where to write JSON outputs (mapping_summary.json,
                 mapq_histogram.json, ground_truth.json)
    --total      total reads in input (denominator for drop_rate)

Outputs:
    mapping_summary.json
    mapq_histogram.json
    ground_truth.json (only when --truth given)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import pysam  # type: ignore
except ImportError as e:
    print("ERROR: pysam not installed. Run: pip install pysam", file=sys.stderr)
    raise


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--truth", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--total", type=int, required=True,
                    help="total reads in input file (denominator for drop_rate)")
    ap.add_argument("--tolerance-bp", type=int, default=10,
                    help="positional tolerance for ground-truth match")
    return ap.parse_args()


def summarize(bam_path: str) -> dict:
    primary = 0
    secondary = 0
    supplementary = 0
    unmapped = 0
    mapq_hist = [0] * 256
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped:
                unmapped += 1
                continue
            if r.is_secondary:
                secondary += 1
                continue
            if r.is_supplementary:
                supplementary += 1
                continue
            primary += 1
            mq = r.mapping_quality if r.mapping_quality is not None else 0
            mq = max(0, min(255, int(mq)))
            mapq_hist[mq] += 1
    return {
        "primary_mapped": primary,
        "secondary_count": secondary,
        "supplementary_count": supplementary,
        "unmapped": unmapped,
        "mapq_hist": mapq_hist,
    }


def load_truth(path: str) -> dict[str, tuple[str, int]]:
    truth: dict[str, tuple[str, int]] = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            read_id, chrom, pos = parts[0], parts[1], int(parts[2])
            truth[read_id] = (chrom, pos)
    return truth


def evaluate_truth(bam_path: str, truth: dict[str, tuple[str, int]],
                   tolerance_bp: int) -> dict:
    tp = 0
    fp = 0
    fn = 0
    seen: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary:
                continue
            name = r.query_name
            seen.add(name)
            if name not in truth:
                continue
            true_chrom, true_pos = truth[name]
            if r.is_unmapped:
                fn += 1
                continue
            called_chrom = bam.get_reference_name(r.reference_id)
            called_pos = r.reference_start
            if called_chrom == true_chrom and abs(called_pos - true_pos) <= tolerance_bp:
                tp += 1
            else:
                fp += 1
    for name in truth:
        if name not in seen:
            fn += 1
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "true_positive": tp,
        "false_positive": fp,
        "false_negative": fn,
        "recall": recall,
        "precision": precision,
        "f1": f1,
        "tolerance_bp": tolerance_bp,
    }


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    s = summarize(args.bam)
    mapped_any = s["primary_mapped"] + s["secondary_count"] + s["supplementary_count"]
    mapping_summary = {
        "total_input_reads": args.total,
        "primary_mapped":    s["primary_mapped"],
        "secondary_count":   s["secondary_count"],
        "supplementary_count": s["supplementary_count"],
        "unmapped":          s["unmapped"],
        "mapping_rate":      s["primary_mapped"] / args.total if args.total else 0.0,
        "drop_rate":         1.0 - (mapped_any / args.total) if args.total else 0.0,
    }
    (out_dir / "mapping_summary.json").write_text(json.dumps(mapping_summary, indent=2))

    hist = {str(i): s["mapq_hist"][i] for i in range(61) if s["mapq_hist"][i]}
    (out_dir / "mapq_histogram.json").write_text(json.dumps(hist, indent=2))

    if args.truth:
        truth = load_truth(args.truth)
        gt = evaluate_truth(args.bam, truth, args.tolerance_bp)
        (out_dir / "ground_truth.json").write_text(json.dumps(gt, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
