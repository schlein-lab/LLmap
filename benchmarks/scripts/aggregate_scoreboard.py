#!/usr/bin/env python3
"""LLmap cross-species benchmark scoreboard aggregator.

Scans ``benchmarks/reports/<organism>/<tier>/summary.json`` (produced by
``run_species_bench.sh`` -> ``analyze_bench.py``) and emits four artefacts
under ``benchmarks/``:

    1. scoreboard.tsv               long table, one row per
                                    (organism, tier, mapper, region_class)
    2. scoreboard_pivoted.tsv       wide table, F1 by mapper
    3. scoreboard_summary.md        human-readable head-to-head + top deltas
    4. scoreboard_failed_regions.json
                                    list of (organism, tier, region_class,
                                    mapper, failure_mode) tuples for any
                                    entry with F1 < 0.5, missing data,
                                    or empty-SAM/missing-BAM tier 5/6 cases.

Design goals
------------
* **Idempotent.**  Re-running on the same reports tree just replaces the
  artefacts in place.  No row duplication, no append mode.
* **Robust to missing data.**  If mapper X didn't run on tier Y, the cell is
  left as ``NaN`` in the pivot and a ``failed_regions`` entry with
  ``failure_mode = "missing"`` is recorded.  No crash.
* **Backwards compatible.**  The script *also* ingests legacy per-run files
  (``mapping_summary.json`` + ``ground_truth.json`` + ``resources.json``)
  emitted by the earlier ``compute.py`` flow, so the scoreboard works on the
  current ``T1/T2`` reports tree as well as future ``run_species_bench.sh``
  output.

Modern summary.json schema (cross-species, per analyze_bench.py)
----------------------------------------------------------------
::

    {
      "per_mapper": {
        "llmap": {
          "mapper": "llmap",
          "tier_dir": "benchmarks/reports/<organism>/<tier>",
          "overall": {
            "n_total": ..., "n_mapped": ..., "mapping_rate": ...,
            "precision_within_1kb": ..., "recall_within_1kb": ...,
            "f1_within_1kb": ..., "f1_within_100bp": ...,
            "mean_mapq": ...
          },
          "by_region_class": {
            "unique": {...overall fields...},
            "SD":     {...overall fields...},
            ...
          }
        },
        "minimap2": {...}
      },
      "mappers_run":     ["llmap","minimap2"],
      "mappers_skipped": []
    }

Per-mapper timing lives next to each BAM in
``<tier_dir>/<mapper>/rep0/timing.json``:
``{wallclock_seconds, user_cpu_seconds, peak_rss_bytes, bam_bytes}``.

Usage
-----
::

    python aggregate_scoreboard.py
    python aggregate_scoreboard.py --reports-dir /custom/reports --out-dir /tmp/scores
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Constants & helpers
# ---------------------------------------------------------------------------

# Long-table column order — keep stable so downstream notebooks can pin headers.
LONG_COLS: list[str] = [
    "organism",
    "tier",
    "mapper",
    "region_class",
    "precision",
    "recall",
    "F1",
    "mapq_correlation",
    "time_sec",
    "peak_mem_mb",
    "n_reads",
]

F1_FAIL_THRESHOLD: float = 0.5
LLMAP_KEY: str = "llmap"

# Region classes we explicitly call out in the markdown summary.
SPECIAL_REGION_CLASSES: list[str] = [
    "SD", "centromere", "telomere", "low_complexity",
    "SV_DEL", "SV_DUP", "SV_INV", "chimera",
]


def _nan_str(v: Any) -> str:
    """Render a value for TSV with stable NaN/None handling."""
    if v is None:
        return "NaN"
    if isinstance(v, float) and math.isnan(v):
        return "NaN"
    if isinstance(v, str):
        # No tabs / no newlines in TSV values, ever.
        return v.replace("\t", " ").replace("\n", " ")
    return str(v)


def _safe_float(v: Any) -> float:
    if v is None:
        return math.nan
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def _peak_rss_to_mb(record: dict[str, Any]) -> float:
    """Pull peak RSS from a record, normalising common spellings to MB."""
    if "peak_mem_mb" in record:
        return _safe_float(record["peak_mem_mb"])
    if "peak_rss_bytes" in record and record["peak_rss_bytes"]:
        return _safe_float(record["peak_rss_bytes"]) / (1024.0 * 1024.0)
    if "max_rss_kb" in record:
        return _safe_float(record["max_rss_kb"]) / 1024.0
    return math.nan


def _load_timing(tier_dir: Path, mapper: str) -> dict[str, Any]:
    """Load <tier_dir>/<mapper>/rep0/timing.json if present, else {}."""
    p = tier_dir / mapper / "rep0" / "timing.json"
    if not p.is_file():
        return {}
    try:
        with p.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[WARN] failed to read {p}: {exc}", file=sys.stderr)
        return {}


def _bam_empty(tier_dir: Path, mapper: str) -> bool:
    """True if the BAM is absent or suspiciously small (< 1024 bytes)."""
    bam = tier_dir / mapper / "rep0" / "alignments.bam"
    if not bam.exists():
        return True
    try:
        return bam.stat().st_size < 1024
    except OSError:
        return True


def _sam_leftover_empty(tier_dir: Path, mapper: str) -> bool:
    """True if a 0-byte aln.sam was left behind (run_llmap empty-SAM warning)."""
    sam = tier_dir / mapper / "rep0" / "aln.sam"
    if not sam.exists():
        return False
    try:
        return sam.stat().st_size == 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Discovery — modern schema
# ---------------------------------------------------------------------------


def _parse_tier_dir(summary_path: Path) -> tuple[str, str, Path]:
    """Extract (organism, tier, tier_dir) from a summary.json path."""
    tier_dir = summary_path.parent
    parts = tier_dir.parts
    # Path: .../reports/<organism>/<tier>/summary.json
    if len(parts) >= 2:
        tier = parts[-1]
        organism = parts[-2]
    else:
        tier = "unknown"
        organism = "unknown"
    if organism == "reports":
        organism = "unknown"
    return organism, tier, tier_dir


def _accuracy_to_row(acc: dict[str, Any], region_class: str) -> dict[str, float]:
    """Project the analyze_bench accuracy dict to (precision, recall, F1, n)."""
    return {
        "precision": _safe_float(acc.get("precision_within_1kb")),
        "recall":    _safe_float(acc.get("recall_within_1kb")),
        "F1":        _safe_float(acc.get("f1_within_1kb")),
        "n_reads":   _safe_float(acc.get("n_total")),
        "mean_mapq": _safe_float(acc.get("mean_mapq")),
        "region_class": region_class,
    }


def _records_from_modern_summary(summary_path: Path) -> list[dict[str, Any]]:
    """Expand a per_mapper summary into many flat records (one per region_class).

    Falls through gracefully when the JSON is the legacy single-mapper
    layout (no ``per_mapper`` key).
    """
    try:
        with summary_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[WARN] failed to read {summary_path}: {exc}", file=sys.stderr)
        return []

    organism, tier, tier_dir = _parse_tier_dir(summary_path)
    records: list[dict[str, Any]] = []

    if "per_mapper" in data:
        per_mapper = data["per_mapper"] or {}
        skipped: list[str] = list(data.get("mappers_skipped") or [])
        for mapper, acc in per_mapper.items():
            timing = _load_timing(tier_dir, mapper)
            time_sec = _safe_float(timing.get("wallclock_seconds"))
            peak_mem = _peak_rss_to_mb(timing)

            overall = acc.get("overall") or {}
            failure_mode: str | None = None
            if _sam_leftover_empty(tier_dir, mapper):
                failure_mode = "empty_sam"
            elif _bam_empty(tier_dir, mapper):
                failure_mode = "missing_bam"

            row = {
                "organism": organism,
                "tier": tier,
                "mapper": mapper,
                **_accuracy_to_row(overall, "all"),
                "mapq_correlation": math.nan,
                "time_sec": time_sec,
                "peak_mem_mb": peak_mem,
                "failure_mode": failure_mode,
                "_source": str(summary_path),
            }
            records.append(row)

            for rclass, sub_acc in (acc.get("by_region_class") or {}).items():
                sub_row = {
                    "organism": organism,
                    "tier": tier,
                    "mapper": mapper,
                    **_accuracy_to_row(sub_acc, rclass),
                    "mapq_correlation": math.nan,
                    "time_sec": time_sec,
                    "peak_mem_mb": peak_mem,
                    "failure_mode": None,
                    "_source": str(summary_path),
                }
                records.append(sub_row)

        # Emit explicit "missing"/"empty_sam" rows for mappers that didn't run.
        for mapper in skipped:
            failure_mode = "empty_sam" if _sam_leftover_empty(tier_dir, mapper) \
                else ("missing_bam" if _bam_empty(tier_dir, mapper) else "missing")
            records.append({
                "organism": organism,
                "tier": tier,
                "mapper": mapper,
                "region_class": "all",
                "precision": math.nan,
                "recall": math.nan,
                "F1": math.nan,
                "n_reads": math.nan,
                "mean_mapq": math.nan,
                "mapq_correlation": math.nan,
                "time_sec": math.nan,
                "peak_mem_mb": math.nan,
                "failure_mode": failure_mode,
                "_source": str(summary_path),
            })
        return records

    # Fallback: very old single-record summary.json (homo_sapiens T1/T2 legacy).
    legacy_record = {
        "organism": data.get("organism", organism),
        "tier": data.get("tier", tier),
        "mapper": data.get("mapper", "unknown"),
        "region_class": data.get("region_class", "all"),
        "precision": _safe_float(data.get("precision")),
        "recall": _safe_float(data.get("recall")),
        "F1": _safe_float(data.get("f1", data.get("F1"))),
        "n_reads": _safe_float(data.get("n_reads")),
        "mean_mapq": _safe_float(data.get("mean_mapq")),
        "mapq_correlation": _safe_float(data.get("mapq_correlation")),
        "time_sec": _safe_float(data.get("time_sec")),
        "peak_mem_mb": _peak_rss_to_mb(data),
        "failure_mode": data.get("failure_mode"),
        "_source": str(summary_path),
    }
    records.append(legacy_record)
    return records


# ---------------------------------------------------------------------------
# Discovery — legacy per-rep
# ---------------------------------------------------------------------------


def _load_legacy_run(run_dir: Path) -> dict[str, Any] | None:
    """Reconstruct a summary from the legacy per-run file triplet."""
    msum = run_dir / "mapping_summary.json"
    gtru = run_dir / "ground_truth.json"
    if not msum.exists():
        return None
    try:
        with msum.open("r", encoding="utf-8") as fh:
            mapping_summary = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None

    gt: dict[str, Any] = {}
    if gtru.exists():
        try:
            with gtru.open("r", encoding="utf-8") as fh:
                gt = json.load(fh)
        except (OSError, json.JSONDecodeError):
            gt = {}

    res: dict[str, Any] = {}
    res_path = run_dir / "resources.json"
    if res_path.exists():
        try:
            with res_path.open("r", encoding="utf-8") as fh:
                res = json.load(fh)
        except (OSError, json.JSONDecodeError):
            res = {}

    # Path: .../reports/<tier>/<mapper>/rep<N>/
    parts = run_dir.parts
    mapper = parts[-2] if len(parts) >= 2 else "unknown"
    tier = parts[-3] if len(parts) >= 3 else "unknown"

    return {
        "organism": "homo_sapiens",
        "tier": tier,
        "mapper": mapper,
        "region_class": "all",
        "precision": _safe_float(gt.get("precision")),
        "recall": _safe_float(gt.get("recall")),
        "F1": _safe_float(gt.get("f1")),
        "n_reads": _safe_float(mapping_summary.get("total_input_reads")),
        "mean_mapq": math.nan,
        "mapq_correlation": math.nan,
        "time_sec": _safe_float(res.get("wallclock_seconds")),
        "peak_mem_mb": _peak_rss_to_mb(res),
        "failure_mode": None,
        "_source": str(run_dir),
    }


def discover_records(reports_dir: Path) -> list[dict[str, Any]]:
    """Walk ``reports_dir`` and gather one record per leaf result.

    Idempotency contract: (organism, tier, mapper, region_class) is unique;
    last-write-wins within a single discovery pass.
    """
    records: list[dict[str, Any]] = []
    seen_keys: dict[tuple[str, str, str, str], int] = {}

    def _push(rec: dict[str, Any]) -> None:
        key = (rec["organism"], rec["tier"], rec["mapper"], rec["region_class"])
        if key in seen_keys:
            # last-write-wins: overwrite the prior index in-place
            records[seen_keys[key]] = rec
        else:
            seen_keys[key] = len(records)
            records.append(rec)

    # 1. Modern cross-species layout (preferred).
    for summary in sorted(reports_dir.rglob("summary.json")):
        for rec in _records_from_modern_summary(summary):
            _push(rec)

    # 2. Legacy per-rep layout — only ingest keys not already covered.
    for msum_path in sorted(reports_dir.rglob("mapping_summary.json")):
        run_dir = msum_path.parent
        rec = _load_legacy_run(run_dir)
        if rec is None:
            continue
        key = (rec["organism"], rec["tier"], rec["mapper"], rec["region_class"])
        if key in seen_keys:
            continue
        _push(rec)

    return records


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def write_long_tsv(records: list[dict[str, Any]], out_path: Path) -> None:
    rows: list[list[str]] = [LONG_COLS]
    # Stable sort: organism, tier, mapper, region_class.
    for r in sorted(records,
                    key=lambda x: (x["organism"], x["tier"],
                                   x["mapper"], x["region_class"])):
        rows.append([_nan_str(r.get(col)) for col in LONG_COLS])
    out_path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def write_pivot_tsv(records: list[dict[str, Any]], out_path: Path) -> None:
    """rows = organism x tier x region_class, cols = mapper, values = F1."""
    mappers = sorted({r["mapper"] for r in records})
    by_key: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    for r in records:
        key = (r["organism"], r["tier"], r["region_class"])
        by_key[key][r["mapper"]] = r["F1"]

    header = ["organism", "tier", "region_class", *mappers]
    rows: list[list[str]] = [header]
    for key in sorted(by_key):
        organism, tier, region_class = key
        row = [organism, tier, region_class]
        for m in mappers:
            row.append(_nan_str(by_key[key].get(m, math.nan)))
        rows.append(row)
    out_path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def write_failed_regions(records: list[dict[str, Any]],
                         out_path: Path,
                         expected_mappers: list[str]) -> None:
    """Emit failures: F1<0.5, explicit failure_mode, or missing combinations."""
    fails: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str, str]] = set()

    def _add(organism: str, tier: str, region_class: str,
             mapper: str, mode: str) -> None:
        key = (organism, tier, region_class, mapper, mode)
        if key in seen:
            return
        seen.add(key)
        fails.append({
            "organism": organism,
            "tier": tier,
            "region_class": region_class,
            "mapper": mapper,
            "failure_mode": mode,
        })

    # 1. Explicit + numeric failures.
    for r in records:
        mode = r.get("failure_mode")
        if mode:
            _add(r["organism"], r["tier"], r["region_class"], r["mapper"], mode)
            continue
        f1 = r["F1"]
        if isinstance(f1, float) and not math.isnan(f1) and f1 < F1_FAIL_THRESHOLD:
            _add(r["organism"], r["tier"], r["region_class"], r["mapper"],
                 f"low_F1({f1:.3f})")

    # 2. Missing-combination failures (only at region_class="all").
    present_all = {
        (r["organism"], r["tier"], r["mapper"])
        for r in records if r["region_class"] == "all"
    }
    organism_tiers = {(r["organism"], r["tier"]) for r in records}
    for organism, tier in sorted(organism_tiers):
        for mapper in expected_mappers:
            if (organism, tier, mapper) not in present_all:
                _add(organism, tier, "all", mapper, "missing")

    # Stable order for readability.
    fails.sort(key=lambda x: (x["organism"], x["tier"],
                              x["region_class"], x["mapper"]))
    out_path.write_text(
        json.dumps(fails, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------


def _f1_by(records: list[dict[str, Any]],
           predicate) -> dict[tuple[str, str, str], dict[str, float]]:
    out: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    for r in records:
        if not predicate(r):
            continue
        f1 = r["F1"]
        if isinstance(f1, float) and math.isnan(f1):
            continue
        key = (r["organism"], r["tier"], r["region_class"])
        out[key][r["mapper"]] = f1
    return out


def _head_to_head(records: list[dict[str, Any]],
                  opponent: str,
                  region_class: str | None) -> tuple[int, int, int]:
    """Return (wins, ties, total) for llmap vs opponent."""
    def pred(r: dict[str, Any]) -> bool:
        if region_class is not None and r["region_class"] != region_class:
            return False
        return True

    by_key = _f1_by(records, pred)
    wins = ties = total = 0
    for key, m in by_key.items():
        if LLMAP_KEY not in m or opponent not in m:
            continue
        total += 1
        if math.isclose(m[LLMAP_KEY], m[opponent], abs_tol=1e-6):
            ties += 1
        elif m[LLMAP_KEY] > m[opponent]:
            wins += 1
    return wins, ties, total


def _top_deltas(records: list[dict[str, Any]],
                n: int = 5,
                positive: bool = True
                ) -> list[tuple[float, str, str, str, str, float, float]]:
    """Return top-N (delta, organism, tier, region_class, best_other, llmap_f1, other_f1)."""
    by_key = _f1_by(records, lambda r: True)
    rows: list[tuple[float, str, str, str, str, float, float]] = []
    for (organism, tier, region_class), m in by_key.items():
        if LLMAP_KEY not in m:
            continue
        others = {k: v for k, v in m.items() if k != LLMAP_KEY}
        if not others:
            continue
        best_other_name = max(others, key=others.get)
        best_other = others[best_other_name]
        delta = m[LLMAP_KEY] - best_other
        rows.append(
            (delta, organism, tier, region_class, best_other_name,
             m[LLMAP_KEY], best_other)
        )
    rows.sort(reverse=positive)
    return rows[:n]


def write_markdown(records: list[dict[str, Any]], out_path: Path) -> None:
    lines: list[str] = []
    lines.append("# LLmap cross-species scoreboard\n")
    organisms_seen = sorted({r["organism"] for r in records})
    tiers_seen = sorted({(r["organism"], r["tier"]) for r in records})
    mappers_seen = sorted({r["mapper"] for r in records})
    lines.append(
        f"Aggregated over **{len(records)}** result rows "
        f"across **{len(organisms_seen)}** organisms, "
        f"**{len(tiers_seen)}** (organism, tier) cells, "
        f"**{len(mappers_seen)}** mappers ({', '.join(mappers_seen)}).\n"
    )

    # Head-to-head summary.
    lines.append("## Head-to-head (LLmap vs other mappers)\n")
    opponents = sorted({r["mapper"] for r in records
                        if r["mapper"] != LLMAP_KEY})
    for opp in opponents:
        for rc in [None, *SPECIAL_REGION_CLASSES]:
            wins, ties, total = _head_to_head(records, opp, rc)
            if total == 0:
                continue
            rc_lbl = "all regions (overall+by_class)" if rc is None \
                else f"region_class={rc}"
            losses = total - wins - ties
            lines.append(
                f"- **vs {opp}** [{rc_lbl}]: "
                f"LLmap wins **{wins}**, ties **{ties}**, "
                f"loses **{losses}** of **{total}** cells."
            )
    lines.append("")

    # Per-organism F1 table at region_class="all".
    lines.append("## Overall F1 by organism x tier (region_class=all)\n")
    lines.append("| organism | tier | " + " | ".join(mappers_seen) + " |")
    lines.append("|---|---|" + "".join(["---:|"] * len(mappers_seen)))
    by_ot: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    for r in records:
        if r["region_class"] != "all":
            continue
        f1 = r["F1"]
        if isinstance(f1, float) and math.isnan(f1):
            continue
        by_ot[(r["organism"], r["tier"])][r["mapper"]] = f1
    for (organism, tier) in sorted(by_ot):
        cells = []
        for m in mappers_seen:
            v = by_ot[(organism, tier)].get(m)
            cells.append("NA" if v is None else f"{v:.3f}")
        lines.append(f"| {organism} | {tier} | " + " | ".join(cells) + " |")
    lines.append("")

    # Top improvements.
    lines.append("## Top-5 LLmap wins (largest F1 gain over best other mapper)\n")
    lines.append("| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |")
    lines.append("|---|---|---|---|---:|---:|---:|")
    for delta, organism, tier, rc, best, lf1, of1 in _top_deltas(records, 5, True):
        lines.append(
            f"| {organism} | {tier} | {rc} | {best} | "
            f"{lf1:.3f} | {of1:.3f} | {delta:+.3f} |"
        )
    lines.append("")

    # Top regressions — feed Phase C tuning list.
    lines.append("## Top-5 LLmap regressions (Phase C tuning targets)\n")
    lines.append("| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |")
    lines.append("|---|---|---|---|---:|---:|---:|")
    for delta, organism, tier, rc, best, lf1, of1 in _top_deltas(records, 5, False):
        if delta >= 0:
            continue
        lines.append(
            f"| {organism} | {tier} | {rc} | {best} | "
            f"{lf1:.3f} | {of1:.3f} | {delta:+.3f} |"
        )
    lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reports-dir", type=Path, default=here / "reports")
    ap.add_argument("--out-dir", type=Path, default=here)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    reports_dir: Path = args.reports_dir
    out_dir: Path = args.out_dir
    if not reports_dir.exists():
        print(f"[ERROR] reports dir not found: {reports_dir}", file=sys.stderr)
        return 2
    out_dir.mkdir(parents=True, exist_ok=True)

    records = discover_records(reports_dir)
    if not records:
        print(f"[WARN] no records discovered under {reports_dir}", file=sys.stderr)

    expected_mappers = sorted({r["mapper"] for r in records})

    write_long_tsv(records, out_dir / "scoreboard.tsv")
    write_pivot_tsv(records, out_dir / "scoreboard_pivoted.tsv")
    write_markdown(records, out_dir / "scoreboard_summary.md")
    write_failed_regions(records,
                         out_dir / "scoreboard_failed_regions.json",
                         expected_mappers)

    print(f"[OK] wrote scoreboard artefacts to {out_dir} "
          f"({len(records)} records, {len(expected_mappers)} mappers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
