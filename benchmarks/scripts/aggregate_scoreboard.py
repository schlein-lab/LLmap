#!/usr/bin/env python3
"""LLmap cross-species benchmark scoreboard aggregator.

Scans ``benchmarks/reports/**/<tier>/summary.json`` (produced by
``run_species_bench.sh``) and emits four artefacts under ``benchmarks/``:

    1. scoreboard.tsv               long table, one row per
                                    (organism, tier, mapper, region_class)
    2. scoreboard_pivoted.tsv       wide table, F1 by mapper
    3. scoreboard_summary.md        human-readable head-to-head + top deltas
    4. scoreboard_failed_regions.json
                                    list of (organism, region, mapper, failure_mode)
                                    tuples for any entry with F1 < 0.5 or missing data.

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

Expected summary.json schema (cross-species)
--------------------------------------------
::

    {
      "organism":          "homo_sapiens",
      "tier":              "T2",
      "mapper":            "llmap",
      "region_class":      "SD",                # SD | centromere | unique | telomere | mito | ...
      "region":            "chr14:105500000-107300000",
      "precision":         0.91,
      "recall":            0.88,
      "f1":                0.895,
      "mapq_correlation":  0.74,
      "time_sec":          312.4,
      "peak_mem_mb":       1840.0,
      "n_reads":           500000,
      "failure_mode":      null                 # or "crash"/"timeout"/"low_recall"
    }

Any field except ``organism``, ``tier`` and ``mapper`` may be absent; the
aggregator inserts ``NaN`` / ``null`` placeholders as appropriate.

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
    "precision",
    "recall",
    "F1",
    "mapq_correlation",
    "region_class",
    "time_sec",
    "peak_mem_mb",
    "n_reads",
]

F1_FAIL_THRESHOLD: float = 0.5
LLMAP_KEY: str = "llmap"

# Region classes we explicitly call out in the markdown summary.  More can be
# added without code changes — they just won't get a dedicated bullet.
SPECIAL_REGION_CLASSES: list[str] = ["SD", "centromere", "telomere", "mito"]


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
    if "peak_rss_bytes" in record:
        return _safe_float(record["peak_rss_bytes"]) / (1024.0 * 1024.0)
    if "max_rss_kb" in record:
        return _safe_float(record["max_rss_kb"]) / 1024.0
    return math.nan


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def _load_modern_summary(summary_path: Path) -> dict[str, Any] | None:
    """Load a ``summary.json`` written by ``run_species_bench.sh``."""
    try:
        with summary_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[WARN] failed to read {summary_path}: {exc}", file=sys.stderr)
        return None

    # Infer organism / tier / mapper from path if missing.
    # Path convention: reports/<organism>/<tier>/<mapper>/summary.json
    parts = summary_path.parts
    inferred = {}
    if len(parts) >= 4:
        inferred["mapper"] = parts[-2]
        inferred["tier"] = parts[-3]
        inferred["organism"] = parts[-4] if parts[-4] != "reports" else "unknown"

    record = {
        "organism": data.get("organism", inferred.get("organism", "unknown")),
        "tier": data.get("tier", inferred.get("tier", "unknown")),
        "mapper": data.get("mapper", inferred.get("mapper", "unknown")),
        "region_class": data.get("region_class", "all"),
        "region": data.get("region", ""),
        "precision": _safe_float(data.get("precision")),
        "recall": _safe_float(data.get("recall")),
        "F1": _safe_float(data.get("f1", data.get("F1"))),
        "mapq_correlation": _safe_float(data.get("mapq_correlation")),
        "time_sec": _safe_float(data.get("time_sec")),
        "peak_mem_mb": _peak_rss_to_mb(data),
        "n_reads": _safe_float(data.get("n_reads")),
        "failure_mode": data.get("failure_mode"),
        "_source": str(summary_path),
    }
    return record


def _load_legacy_run(run_dir: Path) -> dict[str, Any] | None:
    """Reconstruct a summary from the legacy per-run file triplet.

    A legacy run is detected by the presence of ``mapping_summary.json`` and
    a sibling ``ground_truth.json`` (precision/recall/F1) in the same
    directory.  The organism is forced to ``"homo_sapiens"`` because the
    legacy flow never benchmarked anything else.
    """
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
        "region": "",
        "precision": _safe_float(gt.get("precision")),
        "recall": _safe_float(gt.get("recall")),
        "F1": _safe_float(gt.get("f1")),
        "mapq_correlation": math.nan,
        "time_sec": _safe_float(res.get("wallclock_seconds")),
        "peak_mem_mb": _peak_rss_to_mb(res),
        "n_reads": _safe_float(mapping_summary.get("total_input_reads")),
        "failure_mode": None,
        "_source": str(run_dir),
    }


def discover_records(reports_dir: Path) -> list[dict[str, Any]]:
    """Walk ``reports_dir`` and gather one record per leaf result."""
    records: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str, str, str]] = set()

    # 1. Modern cross-species layout (preferred).
    for summary in sorted(reports_dir.rglob("summary.json")):
        rec = _load_modern_summary(summary)
        if rec is None:
            continue
        key = (rec["organism"], rec["tier"], rec["mapper"], rec["region_class"])
        if key in seen_keys:
            # Idempotency: last write wins, but log it.
            print(f"[INFO] duplicate key {key} from {summary}, overwriting",
                  file=sys.stderr)
            records = [r for r in records
                       if (r["organism"], r["tier"], r["mapper"], r["region_class"])
                       != key]
        seen_keys.add(key)
        records.append(rec)

    # 2. Legacy per-rep layout — only ingest keys not already covered.
    for msum_path in sorted(reports_dir.rglob("mapping_summary.json")):
        run_dir = msum_path.parent
        rec = _load_legacy_run(run_dir)
        if rec is None:
            continue
        key = (rec["organism"], rec["tier"], rec["mapper"], rec["region_class"])
        if key in seen_keys:
            continue
        seen_keys.add(key)
        records.append(rec)

    return records


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def write_long_tsv(records: list[dict[str, Any]], out_path: Path) -> None:
    rows: list[list[str]] = [LONG_COLS]
    # Stable sort: organism, tier, region_class, mapper.
    for r in sorted(records,
                    key=lambda x: (x["organism"], x["tier"],
                                   x["region_class"], x["mapper"])):
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
    """Emit failures: F1<0.5, explicit failure_mode, or *missing* combinations."""
    fails: list[dict[str, Any]] = []

    # 1. Explicit / numeric failures.
    for r in records:
        f1 = r["F1"]
        mode = r.get("failure_mode")
        if mode:
            fails.append({
                "organism": r["organism"],
                "region": r.get("region") or r["region_class"],
                "mapper": r["mapper"],
                "failure_mode": mode,
            })
            continue
        if isinstance(f1, float) and not math.isnan(f1) and f1 < F1_FAIL_THRESHOLD:
            fails.append({
                "organism": r["organism"],
                "region": r.get("region") or r["region_class"],
                "mapper": r["mapper"],
                "failure_mode": f"low_F1({f1:.3f})",
            })

    # 2. Missing-combination failures: every (organism, tier, region_class)
    # key should have a row for each expected mapper.
    present: set[tuple[str, str, str, str]] = {
        (r["organism"], r["tier"], r["region_class"], r["mapper"])
        for r in records
    }
    keys: set[tuple[str, str, str]] = {
        (r["organism"], r["tier"], r["region_class"]) for r in records
    }
    for organism, tier, region_class in sorted(keys):
        for mapper in expected_mappers:
            if (organism, tier, region_class, mapper) not in present:
                fails.append({
                    "organism": organism,
                    "region": f"{tier}/{region_class}",
                    "mapper": mapper,
                    "failure_mode": "missing",
                })

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
                  region_class: str | None) -> tuple[int, int]:
    """Count (organism,tier) pairs where llmap beats `opponent`."""
    def pred(r: dict[str, Any]) -> bool:
        if region_class is not None and r["region_class"] != region_class:
            return False
        return True

    by_key = _f1_by(records, pred)
    wins = 0
    total = 0
    for key, m in by_key.items():
        if LLMAP_KEY not in m or opponent not in m:
            continue
        total += 1
        if m[LLMAP_KEY] > m[opponent]:
            wins += 1
    return wins, total


def _top_deltas(records: list[dict[str, Any]],
                n: int = 5,
                positive: bool = True) -> list[tuple[float, str, str, str, str, float, float]]:
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
    lines.append(f"Aggregated over **{len(records)}** result rows.\n")

    # Head-to-head summary.
    lines.append("## Head-to-head\n")
    opponents = sorted({r["mapper"] for r in records
                        if r["mapper"] != LLMAP_KEY})
    for opp in opponents:
        for rc in [None, *SPECIAL_REGION_CLASSES]:
            wins, total = _head_to_head(records, opp, rc)
            if total == 0:
                continue
            rc_lbl = "all regions" if rc is None else f"region_class={rc}"
            lines.append(
                f"- LLmap wins vs {opp} on {wins}/{total} (organism,tier) "
                f"pairs in {rc_lbl}."
            )
    lines.append("")

    # Top improvements.
    lines.append("## Top-5 LLmap improvements (largest F1 gain over best other)\n")
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
