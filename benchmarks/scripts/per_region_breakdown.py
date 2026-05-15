#!/usr/bin/env python3
"""Per-locus accuracy breakdown for cross-species benchmark BAMs.

Given an alignment BAM, a truth TSV, and a ``specific_loci.json`` for the
organism, this script joins each truth read to the named locus it originated
from (by coordinate intersection of its true position) and reports per-locus
hit / miss counts for every mapper present in the report tree.

Inputs
------
* ``--bam`` (one or more)  Path to ``alignments.bam`` files.  Each path is
  expected to live under ``reports/<organism>/<tier>/<mapper>/<rep>/`` so
  that the mapper name can be inferred — pass ``--mapper`` to override.
* ``--truth``              TSV with header ``# read_id<TAB>chrom<TAB>pos``.
  The same format already emitted by ``llmap generate-synth``.
* ``--loci``               JSON file with the named-locus catalogue for the
  organism.  Expected schema::

      [
        {"name": "IGHG4", "chrom": "chr14", "start": 105625000, "end": 105633000},
        ...
      ]

  Either a list (as above) or a dict keyed by name with chrom/start/end
  inside is accepted.
* ``--organism``           Free-text label, recorded in the output for
  cross-organism joins later.
* ``--out``                Output TSV path.  Defaults to
  ``benchmarks/reports/<organism>_per_locus.tsv``.
* ``--tolerance-bp``       Positional tolerance for "correctly placed",
  default 10 bp — matches the legacy ``compute.py`` setting.

Output
------
A long-format TSV, one row per (locus, mapper):

    locus    organism    chrom    start    end    mapper    n_truth_reads    n_placed_correct    accuracy
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

try:
    import pysam  # type: ignore
except ImportError:
    print("ERROR: pysam not installed (pip install pysam)", file=sys.stderr)
    raise


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def load_truth(path: Path) -> dict[str, tuple[str, int]]:
    """Read truth.tsv into {read_id: (chrom, pos)}."""
    out: dict[str, tuple[str, int]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rid, chrom, pos = parts[0], parts[1], parts[2]
            try:
                out[rid] = (chrom, int(pos))
            except ValueError:
                continue
    return out


def load_loci(path: Path) -> list[dict[str, Any]]:
    """Read a ``specific_loci.json`` and normalise to a flat list."""
    with path.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)

    loci: list[dict[str, Any]] = []
    if isinstance(raw, list):
        for item in raw:
            if not isinstance(item, dict):
                continue
            loci.append({
                "name":  item.get("name", item.get("id", "unnamed")),
                "chrom": str(item["chrom"]),
                "start": int(item["start"]),
                "end":   int(item["end"]),
            })
    elif isinstance(raw, dict):
        for name, item in raw.items():
            if not isinstance(item, dict):
                continue
            loci.append({
                "name":  name,
                "chrom": str(item["chrom"]),
                "start": int(item["start"]),
                "end":   int(item["end"]),
            })
    else:
        raise ValueError(f"unrecognised loci schema in {path}")
    return loci


# ---------------------------------------------------------------------------
# Indexing
# ---------------------------------------------------------------------------


def assign_locus(chrom: str, pos: int,
                 loci: list[dict[str, Any]]) -> str | None:
    """Linear scan — N is small (named loci, dozens not millions)."""
    for locus in loci:
        if locus["chrom"] == chrom and locus["start"] <= pos < locus["end"]:
            return locus["name"]
    return None


def truth_per_locus(truth: dict[str, tuple[str, int]],
                    loci: list[dict[str, Any]]
                    ) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Return (locus -> set(read_id), read_id -> locus_name)."""
    per: dict[str, set[str]] = {l["name"]: set() for l in loci}
    rid_to_locus: dict[str, str] = {}
    for rid, (chrom, pos) in truth.items():
        lname = assign_locus(chrom, pos, loci)
        if lname is None:
            continue
        per[lname].add(rid)
        rid_to_locus[rid] = lname
    return per, rid_to_locus


# ---------------------------------------------------------------------------
# BAM scoring
# ---------------------------------------------------------------------------


def infer_mapper(bam_path: Path, override: str | None) -> str:
    if override:
        return override
    # reports/<organism>/<tier>/<mapper>/<rep>/alignments.bam
    parts = bam_path.resolve().parts
    if len(parts) >= 3:
        return parts[-3]
    return bam_path.stem


def score_bam(bam_path: Path,
              truth: dict[str, tuple[str, int]],
              rid_to_locus: dict[str, str],
              tolerance_bp: int) -> dict[str, int]:
    """Return {locus_name: n_correctly_placed}."""
    counts: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            rid = aln.query_name
            if rid not in rid_to_locus:
                continue
            true_chrom, true_pos = truth[rid]
            mapped_chrom = aln.reference_name
            mapped_pos = aln.reference_start
            if mapped_chrom != true_chrom:
                continue
            if abs(mapped_pos - true_pos) > tolerance_bp:
                continue
            lname = rid_to_locus[rid]
            counts[lname] = counts.get(lname, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_tsv(rows: list[dict[str, Any]], out_path: Path) -> None:
    cols = ["locus", "organism", "chrom", "start", "end",
            "mapper", "n_truth_reads", "n_placed_correct", "accuracy"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (x["locus"], x["mapper"])):
            vals = []
            for c in cols:
                v = r.get(c)
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    vals.append("NaN")
                elif isinstance(v, float):
                    vals.append(f"{v:.6f}")
                else:
                    vals.append(str(v).replace("\t", " "))
            fh.write("\t".join(vals) + "\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", action="append", required=True,
                    help="alignment BAM (repeatable; one per mapper)")
    ap.add_argument("--truth", type=Path, required=True)
    ap.add_argument("--loci", type=Path, required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--mapper", default=None,
                    help="override mapper name (only sensible when --bam used once)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--tolerance-bp", type=int, default=10)
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    if not args.truth.exists():
        print(f"[ERROR] truth not found: {args.truth}", file=sys.stderr)
        return 2
    if not args.loci.exists():
        print(f"[ERROR] loci not found: {args.loci}", file=sys.stderr)
        return 2

    truth = load_truth(args.truth)
    loci = load_loci(args.loci)
    per_locus_reads, rid_to_locus = truth_per_locus(truth, loci)

    rows: list[dict[str, Any]] = []
    locus_meta = {l["name"]: l for l in loci}

    for bam_str in args.bam:
        bam_path = Path(bam_str)
        if not bam_path.exists():
            print(f"[WARN] BAM missing, skipping: {bam_path}", file=sys.stderr)
            continue
        mapper = infer_mapper(bam_path, args.mapper)
        counts = score_bam(bam_path, truth, rid_to_locus, args.tolerance_bp)
        for lname, meta in locus_meta.items():
            n_truth = len(per_locus_reads.get(lname, set()))
            n_correct = counts.get(lname, 0)
            accuracy = (n_correct / n_truth) if n_truth > 0 else math.nan
            rows.append({
                "locus":            lname,
                "organism":         args.organism,
                "chrom":            meta["chrom"],
                "start":            meta["start"],
                "end":              meta["end"],
                "mapper":           mapper,
                "n_truth_reads":    n_truth,
                "n_placed_correct": n_correct,
                "accuracy":         accuracy,
            })

    out_path = args.out
    if out_path is None:
        here = Path(__file__).resolve().parent.parent
        out_path = here / "reports" / f"{args.organism}_per_locus.tsv"
    write_tsv(rows, out_path)
    print(f"[OK] wrote {len(rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
