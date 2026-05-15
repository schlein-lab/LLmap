#!/usr/bin/env python3
"""
analyze_bench.py — per-read accuracy + per-region-class accuracy for one
tier directory laid out by run_species_bench.sh:

  <tier_dir>/
    truth.tsv
    <mapper>/rep0/alignments.bam
    ...

Outputs:
  <tier_dir>/<mapper>/rep0/accuracy.json
  <tier_dir>/comparison.tsv     # one row per read, each mapper's call
  <tier_dir>/summary.json

BUGFIX 2026-05-15: load_truth() previously stored (source_start+source_end)//2 as
truth["pos"], while score_one() treated truth["pos"] as the *start* of the source
interval. The effective truth interval became [midpoint, end] instead of
[start, end], so reads aligned to their true source start were scored as off by
~half the read length (multi-kb for HiFi tiers). Result: F1@1kb=0.000 across
all non-human organisms with long-read tiers. truth["pos"] now holds the actual
source_start; truth["end"] is source_end. BAM coords were correct all along —
this was a truth-loading bug, not a mapper bug and not a coord-space mismatch
between the fake reference build and gen_synth_reads (those are consistent).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    print("ERROR: pysam not installed (pip install pysam)", file=sys.stderr)
    sys.exit(2)


# gen_synth_reads.py produces these classes (and SV_*/chimera variants).
REGION_CLASSES = [
    "unique", "SD", "centromere", "telomere", "low_complexity", "synth_stress",
    "SV_DEL", "SV_DUP", "SV_INV", "chimera",
    # legacy lowercase fallbacks
    "sd",
]

WINDOWS = {"within_100bp": 100, "within_1kb": 1_000}


def load_truth(path: Path) -> Dict[str, Dict]:
    """read_id -> {chrom, pos, strand, region_class}

    Accepts two layouts:
      (legacy) read_id  chrom  pos                    [strand] [region_class]
      (new)    read_id  source_contig  source_start  source_end  true_strand  read_class  subclass
    The header may start with '#' or be plain text; we sniff by column count.
    """
    out: Dict[str, Dict] = {}
    with open(path) as f:
        header_seen = False
        for ln in f:
            if not ln.strip():
                continue
            if ln.startswith("#") or (not header_seen and ln.lower().startswith("read_id")):
                header_seen = True
                continue
            header_seen = True
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rid = parts[0]
            chrom = parts[1]
            if len(parts) >= 6:
                # new schema: start, end, strand, class
                # IMPORTANT: source_start/source_end are fake-reference coords
                # (the same space the mapper sees), so we store source_start as
                # pos. score_one() treats truth["pos"]..truth["end"] as the
                # closed interval and computes distance to the nearest edge.
                try:
                    start = int(parts[2]); end = int(parts[3])
                except ValueError:
                    continue
                pos = start
                strand = parts[4]
                cls = parts[5]
                end_val = end
            else:
                try:
                    pos = int(parts[2])
                except ValueError:
                    continue
                strand = parts[3] if len(parts) > 3 else "+"
                cls = parts[4] if len(parts) > 4 else "unique"
                end_val = pos
            out[rid] = {"chrom": chrom, "pos": pos, "end": end_val,
                        "strand": strand, "region_class": cls}
    return out


def parse_bam(bam_path: Path) -> Dict[str, Dict]:
    """
    Return read_id -> {chrom, pos, strand, mapq, is_mapped, is_secondary}
    Picks the primary alignment per read.
    """
    calls: Dict[str, Dict] = {}
    if not bam_path.is_file():
        return calls
    bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
    try:
        for rec in bam.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            rid = rec.query_name
            if rid is None:
                continue
            # Some pipelines append /1 /2; strip mate suffix for comparison
            if rid.endswith("/1") or rid.endswith("/2"):
                rid = rid[:-2]
            if rec.is_unmapped:
                # only record if we haven't already
                if rid not in calls:
                    calls[rid] = {
                        "chrom": None, "pos": None, "strand": None,
                        "mapq": 0, "is_mapped": False,
                    }
                continue
            calls[rid] = {
                "chrom": rec.reference_name,
                "pos": int(rec.reference_start),
                "strand": "-" if rec.is_reverse else "+",
                "mapq": int(rec.mapping_quality),
                "is_mapped": True,
            }
    finally:
        bam.close()
    return calls


def score_one(truth: Dict, call: Optional[Dict]) -> Dict:
    if not call or not call.get("is_mapped"):
        return {
            "mapped": False, "correct_chrom": False,
            "within_100bp": False, "within_1kb": False,
            "correct_strand": False, "mapq": 0,
            "abs_dist": None,
        }
    correct_chrom = (call["chrom"] == truth["chrom"])
    abs_dist = None
    if correct_chrom:
        truth_start = int(truth["pos"])
        truth_end = int(truth.get("end", truth_start))
        called = int(call["pos"])
        # distance to nearest edge of the truth interval (0 if inside)
        if called < truth_start:
            abs_dist = truth_start - called
        elif called > truth_end:
            abs_dist = called - truth_end
        else:
            abs_dist = 0
    return {
        "mapped": True,
        "correct_chrom": correct_chrom,
        "within_100bp": correct_chrom and abs_dist is not None and abs_dist <= 100,
        "within_1kb":   correct_chrom and abs_dist is not None and abs_dist <= 1_000,
        "correct_strand": call.get("strand") == truth.get("strand"),
        "mapq": int(call.get("mapq", 0)),
        "abs_dist": abs_dist,
    }


def aggregate(rows: List[Dict]) -> Dict:
    """Compute precision/recall/F1 with within_1kb as the 'correct' threshold."""
    n_total = len(rows)
    n_mapped = sum(1 for r in rows if r["mapped"])
    n_correct_100 = sum(1 for r in rows if r["within_100bp"])
    n_correct_1k  = sum(1 for r in rows if r["within_1kb"])
    n_correct_strand = sum(1 for r in rows if r["within_1kb"] and r["correct_strand"])

    precision_1k = (n_correct_1k / n_mapped) if n_mapped else 0.0
    recall_1k    = (n_correct_1k / n_total)  if n_total else 0.0
    f1_1k = (2 * precision_1k * recall_1k / (precision_1k + recall_1k)) if (precision_1k + recall_1k) else 0.0

    precision_100 = (n_correct_100 / n_mapped) if n_mapped else 0.0
    recall_100    = (n_correct_100 / n_total)  if n_total else 0.0
    f1_100 = (2 * precision_100 * recall_100 / (precision_100 + recall_100)) if (precision_100 + recall_100) else 0.0

    mapqs = [r["mapq"] for r in rows if r["mapped"]]
    return {
        "n_total": n_total,
        "n_mapped": n_mapped,
        "mapping_rate": (n_mapped / n_total) if n_total else 0.0,
        "n_correct_within_100bp": n_correct_100,
        "n_correct_within_1kb": n_correct_1k,
        "n_correct_strand_within_1kb": n_correct_strand,
        "precision_within_100bp": precision_100,
        "recall_within_100bp": recall_100,
        "f1_within_100bp": f1_100,
        "precision_within_1kb": precision_1k,
        "recall_within_1kb": recall_1k,
        "f1_within_1kb": f1_1k,
        "mean_mapq": (sum(mapqs) / len(mapqs)) if mapqs else 0.0,
    }


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description="Score mapper BAMs against truth.tsv for one tier directory.")
    ap.add_argument("--tier-dir", required=True, type=Path,
                    help="benchmarks/reports/<organism>/<tier>/")
    ap.add_argument("--mappers", required=True,
                    help="Comma-separated mapper names to consider")
    args = ap.parse_args(argv)

    tier_dir: Path = args.tier_dir
    truth_path = tier_dir / "truth.tsv"
    if not truth_path.is_file():
        print(f"ERROR: truth.tsv not found in {tier_dir}", file=sys.stderr)
        return 2
    truth = load_truth(truth_path)
    if not truth:
        print("ERROR: empty truth", file=sys.stderr)
        return 2

    mappers = [m.strip() for m in args.mappers.split(",") if m.strip()]
    per_mapper_calls: Dict[str, Dict[str, Dict]] = {}
    per_mapper_rows: Dict[str, List[Dict]] = {}

    summary: Dict[str, Dict] = {"per_mapper": {}, "mappers_run": [], "mappers_skipped": []}

    for mapper in mappers:
        bam = tier_dir / mapper / "rep0" / "alignments.bam"
        if not bam.is_file():
            print(f"[analyze_bench] {mapper}: no BAM at {bam} -> skipping", file=sys.stderr)
            summary["mappers_skipped"].append(mapper)
            continue
        calls = parse_bam(bam)
        per_mapper_calls[mapper] = calls

        rows: List[Dict] = []
        for rid, t in truth.items():
            score = score_one(t, calls.get(rid))
            score["read_id"] = rid
            score["region_class"] = t.get("region_class", "unique")
            rows.append(score)
        per_mapper_rows[mapper] = rows

        overall = aggregate(rows)
        by_class: Dict[str, Dict] = {}
        for cls in REGION_CLASSES:
            sub = [r for r in rows if r["region_class"] == cls]
            if sub:
                by_class[cls] = aggregate(sub)
        acc = {
            "mapper": mapper,
            "tier_dir": str(tier_dir),
            "overall": overall,
            "by_region_class": by_class,
        }
        out_path = tier_dir / mapper / "rep0" / "accuracy.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(acc, f, indent=2)
        print(f"[analyze_bench] wrote {out_path} "
              f"(f1@1kb={overall['f1_within_1kb']:.3f}, map_rate={overall['mapping_rate']:.3f})",
              file=sys.stderr)
        summary["per_mapper"][mapper] = acc
        summary["mappers_run"].append(mapper)

    # ---- comparison.tsv -------------------------------------------------------
    comp_path = tier_dir / "comparison.tsv"
    mapped_mappers = summary["mappers_run"]
    header = ["read_id", "truth_chrom", "truth_pos", "truth_strand", "region_class"]
    for m in mapped_mappers:
        header += [f"{m}_chrom", f"{m}_pos", f"{m}_strand", f"{m}_mapq",
                   f"{m}_within_100bp", f"{m}_within_1kb"]
    with open(comp_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for rid, t in truth.items():
            row = [rid, t["chrom"], str(t["pos"]), t["strand"], t.get("region_class", "unique")]
            for m in mapped_mappers:
                call = per_mapper_calls.get(m, {}).get(rid) or {}
                score = score_one(t, call)
                row += [
                    str(call.get("chrom") or ""),
                    str(call.get("pos") if call.get("pos") is not None else ""),
                    str(call.get("strand") or ""),
                    str(call.get("mapq", 0)),
                    str(score["within_100bp"]).lower(),
                    str(score["within_1kb"]).lower(),
                ]
            f.write("\t".join(row) + "\n")
    print(f"[analyze_bench] wrote {comp_path}", file=sys.stderr)

    summary_path = tier_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[analyze_bench] wrote {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
