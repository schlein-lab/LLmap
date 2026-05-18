#!/usr/bin/env python3
"""
build.py - Regenerate losslessmap.com static site from LLmap benchmark data.

Sources:
    benchmarks/scoreboard.tsv         (long TSV: organism, tier, mapper, region_class, F1, ...)
    benchmarks/scoreboard_pivoted.tsv (pivoted: organism, tier, region_class, mapper cols)
    benchmarks/scoreboard_summary.md  (highlights)
    benchmarks/scoreboard_failed_regions.json
    benchmarks/reports/<organism>/<tier>/summary.json
    benchmarks/reports/<organism>/<tier>/<mapper>/rep0/timing.json

Outputs:
    web/losslessmap.com/index.html
    web/losslessmap.com/scoreboard.html
    web/losslessmap.com/methodology.html
    web/losslessmap.com/organism/<name>.html
    web/losslessmap.com/styles.css
    web/losslessmap.com/data/scoreboard.json
    web/losslessmap.com/data/scoreboard.tsv

No external dependencies. Pure stdlib + f-strings.
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from html import escape
from pathlib import Path

# --------------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent              # web/losslessmap.com
REPO = ROOT.parent.parent                           # llmap-local/
BENCH = REPO / "benchmarks"
REPORTS = BENCH / "reports"
DATA_OUT = ROOT / "data"
ORG_OUT = ROOT / "organism"

SCOREBOARD_TSV = BENCH / "scoreboard.tsv"
PIVOTED_TSV = BENCH / "scoreboard_pivoted.tsv"
SUMMARY_MD = BENCH / "scoreboard_summary.md"
FAILED_JSON = BENCH / "scoreboard_failed_regions.json"

MAPPERS = ["llmap", "minimap2", "bwa-mem2", "winnowmap"]
TIER_ORDER = ["tier1", "tier2", "tier5", "tier6", "tier10", "T1", "T2"]
VALID_TIER_RE = re.compile(r"^(tier\d+|T\d+)$")
VALID_REGION_CLASSES = {
    "all", "unique", "SD", "centromere", "telomere", "low_complexity",
    "SV_DEL", "SV_DUP", "SV_INV",
}


# --------------------------------------------------------------------------- helpers
def parse_float(s: str) -> float | None:
    if s is None or s == "" or s == "NaN" or s == "NA":
        return None
    try:
        v = float(s)
        if math.isnan(v):
            return None
        return v
    except ValueError:
        return None


def load_scoreboard() -> list[dict]:
    rows = []
    with SCOREBOARD_TSV.open() as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            # skip header-as-row drift / invalid
            if r.get("organism") == "organism":
                continue
            # skip rows with mis-aligned columns (legacy schema where tier is a
            # bare integer and region_class is an int n_total).
            if not VALID_TIER_RE.match(r.get("tier", "")):
                continue
            if r.get("region_class", "") not in VALID_REGION_CLASSES:
                continue
            rec = {
                "organism": r["organism"],
                "tier": r["tier"],
                "mapper": r["mapper"],
                "region_class": r["region_class"],
                "precision": parse_float(r.get("precision", "")),
                "recall": parse_float(r.get("recall", "")),
                "F1": parse_float(r.get("F1", "")),
                "mapq_correlation": parse_float(r.get("mapq_correlation", "")),
                "time_sec": parse_float(r.get("time_sec", "")),
                "peak_mem_mb": parse_float(r.get("peak_mem_mb", "")),
                "n_reads": parse_float(r.get("n_reads", "")),
            }
            rows.append(rec)
    return rows


def f1_bucket(v: float | None) -> tuple[str, str]:
    """Return (css_class, glyph) for an F1 value."""
    if v is None:
        return ("na", "&mdash;")
    if v >= 0.99:
        return ("excellent", "✓")  # check
    if v >= 0.90:
        return ("good", "✓")
    if v >= 0.70:
        return ("warn", "⚠")       # warning
    return ("bad", "✗")            # cross


def fmt_f1(v: float | None) -> str:
    if v is None:
        return "&mdash;"
    return f"{v:.3f}"


def tier_sort_key(t: str) -> tuple[int, int, str]:
    """Sort: tier-prefix first (lowercase 'tier' group 0, 'T' group 1), then numeric."""
    m = re.match(r"^(tier|T)(\d+)$", t)
    if m:
        group = 0 if m.group(1) == "tier" else 1
        return (group, int(m.group(2)), t)
    return (2, 999, t)


# --------------------------------------------------------------------------- shared HTML chrome
CSS_HREF = "styles.css"


def page_shell(title: str, body: str, css_href: str = CSS_HREF, active: str = "") -> str:
    nav = []
    for href, label, key in [
        ("index.html", "Overview", "index"),
        ("scoreboard.html", "Scoreboard", "scoreboard"),
        ("methodology.html", "Methodology", "methodology"),
    ]:
        cls = ' class="active"' if active == key else ""
        nav.append(f'<a href="{href}"{cls}>{label}</a>')
    nav_html = " &middot; ".join(nav)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{escape(title)} &middot; losslessmap.com</title>
<link rel="stylesheet" href="{css_href}" />
</head>
<body>
<header>
  <h1><a href="index.html">losslessmap.com</a></h1>
  <p class="tagline">LLmap cross-species mapper scoreboard.</p>
  <nav>{nav_html}
    <span class="theme-toggle">
      <button id="theme-toggle" aria-label="Toggle dark mode">dark/light</button>
    </span>
  </nav>
</header>
<main>
{body}
</main>
<footer>
  <p>LLmap &mdash; open-source wave-particle adaptive mapper. GPU + AI optional. Per-region Bayesian priors. <a href="https://github.com/schlein-lab/LLmap">github.com/schlein-lab/LLmap</a></p>
  <p class="muted">Generated by <code>web/losslessmap.com/build.py</code> from <code>benchmarks/scoreboard.tsv</code>.</p>
</footer>
<script>
  // Theme toggle: persists in localStorage.
  (function() {{
    var root = document.documentElement;
    var stored = localStorage.getItem('llmap-theme');
    if (stored) root.setAttribute('data-theme', stored);
    var btn = document.getElementById('theme-toggle');
    if (btn) {{
      btn.addEventListener('click', function() {{
        var cur = root.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
        root.setAttribute('data-theme', cur);
        localStorage.setItem('llmap-theme', cur);
      }});
    }}
  }})();
</script>
</body>
</html>
"""


# --------------------------------------------------------------------------- styles
STYLES_CSS = """\
/* losslessmap.com - LLmap scoreboard - clean monospaced layout */
:root {
  --bg: #0f1115;
  --fg: #e8e8e8;
  --muted: #888;
  --accent: #6cf;
  --excellent: #2ecc71;
  --good: #f1c40f;
  --warn: #e67e22;
  --bad: #e74c3c;
  --na: #555;
  --row-alt: #161a20;
  --border: #2a2f38;
  --link: #6cf;
}
html[data-theme="light"] {
  --bg: #fafafa;
  --fg: #1a1a1a;
  --muted: #666;
  --accent: #0066cc;
  --excellent: #1f9d55;
  --good: #b88a00;
  --warn: #c25e00;
  --bad: #b3261e;
  --na: #aaa;
  --row-alt: #f0f0f0;
  --border: #d0d0d0;
  --link: #0066cc;
}
* { box-sizing: border-box; }
html, body {
  margin: 0;
  padding: 0;
  background: var(--bg);
  color: var(--fg);
  font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
  font-size: 14px;
  line-height: 1.5;
}
header {
  border-bottom: 1px solid var(--border);
  padding: 1.2em 2em 0.8em;
}
header h1 {
  margin: 0 0 0.2em;
  font-size: 1.4em;
}
header h1 a { color: var(--fg); text-decoration: none; }
header .tagline { margin: 0; color: var(--muted); font-size: 0.95em; }
header nav { margin-top: 0.6em; font-size: 0.95em; }
header nav a { color: var(--link); text-decoration: none; margin-right: 0.4em; }
header nav a.active { color: var(--accent); border-bottom: 1px solid var(--accent); }
header nav .theme-toggle { float: right; }
header nav .theme-toggle button {
  background: transparent; color: var(--muted);
  border: 1px solid var(--border); padding: 2px 8px; cursor: pointer;
  font-family: inherit; font-size: 0.9em; border-radius: 3px;
}
main { padding: 1.5em 2em; max-width: 1400px; }
footer { border-top: 1px solid var(--border); padding: 1em 2em; color: var(--muted); font-size: 0.9em; }
footer a { color: var(--link); }
h2 { border-bottom: 1px solid var(--border); padding-bottom: 0.3em; margin-top: 2em; }
h3 { color: var(--accent); margin-top: 1.5em; }
a { color: var(--link); }
code, pre { background: var(--row-alt); padding: 0 4px; border-radius: 3px; }
pre { padding: 0.8em 1em; overflow-x: auto; }
.muted { color: var(--muted); }
.summary-cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 1em;
  margin: 1em 0 2em;
}
.summary-card {
  border: 1px solid var(--border);
  padding: 0.8em 1em;
  border-radius: 4px;
  background: var(--row-alt);
}
.summary-card .num {
  font-size: 1.6em;
  font-weight: bold;
  color: var(--accent);
  display: block;
}
.summary-card .lbl {
  color: var(--muted);
  font-size: 0.85em;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 1em 0;
  font-size: 0.92em;
}
th, td {
  text-align: left;
  padding: 4px 8px;
  border-bottom: 1px solid var(--border);
}
th {
  background: var(--row-alt);
  font-weight: bold;
  color: var(--muted);
  cursor: pointer;
  user-select: none;
}
th.sortable:hover { color: var(--accent); }
tr:nth-child(even) td { background: var(--row-alt); }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
.cell-excellent { background: var(--excellent); color: #fff; }
.cell-good      { background: var(--good); color: #000; }
.cell-warn      { background: var(--warn); color: #000; }
.cell-bad       { background: var(--bad); color: #fff; }
.cell-na        { background: var(--na); color: #ddd; }
.cell-excellent, .cell-good, .cell-warn, .cell-bad, .cell-na {
  text-align: center;
  font-variant-numeric: tabular-nums;
  font-weight: bold;
}
.heatmap { font-size: 0.85em; }
.heatmap th, .heatmap td { padding: 4px 6px; min-width: 64px; }
.heatmap th.org-col { text-align: left; min-width: 140px; }
.heatmap td.org-col { text-align: left; font-weight: bold; }
.heatmap a { color: inherit; text-decoration: none; display: block; }
.legend { margin: 0.5em 0 1.5em; font-size: 0.85em; }
.legend span {
  display: inline-block; padding: 2px 8px; margin-right: 0.3em;
  border-radius: 3px; font-weight: bold;
}
.win  { color: var(--excellent); }
.loss { color: var(--bad); }
.tie  { color: var(--muted); }
.delta-pos { color: var(--excellent); font-weight: bold; }
.delta-neg { color: var(--bad); font-weight: bold; }
.cite {
  border-left: 3px solid var(--accent);
  padding: 0.6em 1em;
  background: var(--row-alt);
  margin: 1.5em 0;
  font-style: italic;
}
.cite strong { font-style: normal; color: var(--accent); }
ul.compact { line-height: 1.7; }
.controls { margin: 1em 0; }
.controls input, .controls select {
  background: var(--bg); color: var(--fg); border: 1px solid var(--border);
  padding: 4px 8px; font-family: inherit; font-size: inherit; border-radius: 3px;
}
"""


# --------------------------------------------------------------------------- index.html
def parse_summary_md(md: str) -> dict:
    """Pull out structured pieces from scoreboard_summary.md."""
    out: dict = {
        "header_paragraph": "",
        "head_to_head": [],   # list of (description, wins, ties, losses, total)
        "wins": [],           # list of dict rows
        "regressions": [],
    }
    lines = md.splitlines()

    # First paragraph after H1.
    for i, line in enumerate(lines):
        if line.startswith("Aggregated"):
            out["header_paragraph"] = line.strip()
            break

    # Head-to-head bullets.
    for line in lines:
        m = re.match(
            r"-\s+\*\*vs\s+(\w+)\*\*\s*\[(.*?)\]:\s+LLmap wins\s+\*\*(\d+)\*\*,\s+ties\s+\*\*(\d+)\*\*,\s+loses\s+\*\*(\d+)\*\*\s+of\s+\*\*(\d+)\*\*",
            line.strip(),
        )
        if m:
            out["head_to_head"].append({
                "other": m.group(1),
                "scope": m.group(2),
                "wins": int(m.group(3)),
                "ties": int(m.group(4)),
                "losses": int(m.group(5)),
                "total": int(m.group(6)),
            })

    # Wins / regressions tables (extract rows from markdown tables).
    def parse_md_table(start_idx: int) -> list[dict]:
        rows: list[dict] = []
        # Skip header + separator.
        i = start_idx
        # find header row
        while i < len(lines) and "|" not in lines[i]:
            i += 1
        if i >= len(lines):
            return rows
        header = [c.strip() for c in lines[i].strip().strip("|").split("|")]
        i += 2  # skip separator
        while i < len(lines) and lines[i].strip().startswith("|"):
            cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            row = dict(zip(header, cells))
            rows.append(row)
            i += 1
        return rows

    for i, line in enumerate(lines):
        if line.startswith("## Top-5 LLmap wins"):
            out["wins"] = parse_md_table(i + 1)
        elif line.startswith("## Top-5 LLmap regressions"):
            out["regressions"] = parse_md_table(i + 1)
    return out


def render_index(rows: list[dict], summary: dict, organisms: list[str], tiers: list[str]) -> str:
    # Headline win.
    headline = ""
    if summary["wins"]:
        w = summary["wins"][0]
        headline = (
            f'LLmap WINS by {w.get("delta", "?")} F1 on '
            f'<code>{escape(w.get("region_class", "?"))}</code> mapping '
            f'({escape(w.get("organism", "?"))} {escape(w.get("tier", "?"))})'
        )

    # Summary cards.
    # Compute row counts.
    n_rows = len(rows)
    n_orgs = len(organisms)
    n_cells = len(set((r["organism"], r["tier"]) for r in rows))
    n_mappers = len(set(r["mapper"] for r in rows))

    cards = f"""
<div class="summary-cards">
  <div class="summary-card"><span class="num">{n_orgs}</span><span class="lbl">organisms</span></div>
  <div class="summary-card"><span class="num">{n_cells}</span><span class="lbl">(org, tier) cells</span></div>
  <div class="summary-card"><span class="num">{n_mappers}</span><span class="lbl">mappers</span></div>
  <div class="summary-card"><span class="num">{n_rows}</span><span class="lbl">result rows</span></div>
</div>
"""

    # Head-to-head section.
    h2h_rows = []
    for h in summary["head_to_head"]:
        h2h_rows.append(
            f"<tr><td>vs <strong>{escape(h['other'])}</strong></td>"
            f"<td>{escape(h['scope'])}</td>"
            f"<td class='num win'>{h['wins']}</td>"
            f"<td class='num tie'>{h['ties']}</td>"
            f"<td class='num loss'>{h['losses']}</td>"
            f"<td class='num'>{h['total']}</td></tr>"
        )
    h2h_table = (
        '<table><thead><tr><th>matchup</th><th>scope</th><th class="num">wins</th>'
        '<th class="num">ties</th><th class="num">losses</th><th class="num">total</th></tr></thead>'
        f'<tbody>{"".join(h2h_rows)}</tbody></table>'
    )

    # Top wins/regressions tables.
    def render_delta_table(rows_md: list[dict], pos: bool) -> str:
        cls = "delta-pos" if pos else "delta-neg"
        out = [
            "<table><thead><tr>"
            "<th>organism</th><th>tier</th><th>region_class</th>"
            "<th>vs</th><th class='num'>LLmap F1</th><th class='num'>other F1</th><th class='num'>delta</th>"
            "</tr></thead><tbody>"
        ]
        for r in rows_md:
            org = r.get("organism", "")
            org_link = f'<a href="organism/{escape(org)}.html">{escape(org)}</a>' if org else ""
            out.append(
                f"<tr><td>{org_link}</td>"
                f"<td>{escape(r.get('tier',''))}</td>"
                f"<td>{escape(r.get('region_class',''))}</td>"
                f"<td>{escape(r.get('best_other',''))}</td>"
                f"<td class='num'>{escape(r.get('LLmap F1',''))}</td>"
                f"<td class='num'>{escape(r.get('other F1',''))}</td>"
                f"<td class='num {cls}'>{escape(r.get('delta',''))}</td></tr>"
            )
        out.append("</tbody></table>")
        return "".join(out)

    # Heatmap: organism x tier with F1 of LLmap (region_class=all).
    f1_lookup = {}
    for r in rows:
        if r["region_class"] != "all":
            continue
        f1_lookup[(r["organism"], r["tier"], r["mapper"])] = r["F1"]

    # Build two heatmaps: LLmap and minimap2.
    def render_heatmap(mapper: str, organisms: list[str], tiers: list[str]) -> str:
        head = [
            f"<table class='heatmap'><thead><tr><th class='org-col'>organism &times; tier</th>"
        ]
        for t in tiers:
            head.append(f"<th>{escape(t)}</th>")
        head.append("</tr></thead><tbody>")
        body = []
        for org in organisms:
            body.append(f"<tr><td class='org-col'><a href='organism/{escape(org)}.html'>{escape(org)}</a></td>")
            for t in tiers:
                v = f1_lookup.get((org, t, mapper))
                cls, glyph = f1_bucket(v)
                txt = fmt_f1(v)
                body.append(f"<td class='cell-{cls}' title='{escape(org)} {escape(t)} {mapper}'>{glyph} {txt}</td>")
            body.append("</tr>")
        body.append("</tbody></table>")
        return "".join(head) + "".join(body)

    legend = """
<div class="legend">
  <span class="cell-excellent">&#10003; F1 &ge; 0.99</span>
  <span class="cell-good">&#10003; F1 &ge; 0.90</span>
  <span class="cell-warn">&#9888; F1 &ge; 0.70</span>
  <span class="cell-bad">&#10007; F1 &lt; 0.70</span>
  <span class="cell-na">&mdash; NA</span>
</div>
"""

    body = f"""
<section>
  <p class="muted">{escape(summary["header_paragraph"])}</p>
  {("<p class='cite'><strong>Headline:</strong> " + headline + "</p>") if headline else ""}
  {cards}
</section>

<section>
  <h2>Head-to-head: LLmap vs other mappers</h2>
  {h2h_table}
</section>

<section>
  <h2>Top-5 LLmap wins</h2>
  <p class="muted">Largest F1 gain over the best non-LLmap mapper.</p>
  {render_delta_table(summary["wins"], pos=True)}
</section>

<section>
  <h2>Top-5 LLmap regressions (Phase C tuning queue)</h2>
  <p class="muted">Largest F1 deficits &mdash; current focus for Phase C tuning.</p>
  {render_delta_table(summary["regressions"], pos=False)}
</section>

<section>
  <h2>Cross-species F1 heatmap (region_class = all)</h2>
  {legend}
  <h3>LLmap</h3>
  {render_heatmap("llmap", organisms, tiers)}
  <h3>minimap2 (baseline)</h3>
  {render_heatmap("minimap2", organisms, tiers)}
</section>

<section>
  <h2>About LLmap</h2>
  <p class="cite">
    <strong>LLmap</strong> is the open-source mapper with native wave-particle adaptive &lambda;,
    GPU + AI optional, and per-region Bayesian priors. It targets long-read mapping in
    high-noise, high-SV, and otherwise difficult regions where minimap2's seed-chain-extend
    pipeline degrades (metagenomes, structural variants, repetitive plant/rodent genomes).
  </p>
  <p>Source: <a href="https://github.com/schlein-lab/LLmap">github.com/schlein-lab/LLmap</a>.
     Methodology: <a href="methodology.html">methodology.html</a>.</p>
</section>
"""
    return page_shell("Overview", body, active="index")


# --------------------------------------------------------------------------- scoreboard.html
def render_scoreboard(rows: list[dict]) -> str:
    # Pre-format rows sorted by organism, tier, mapper, region_class.
    s = sorted(
        rows,
        key=lambda r: (r["organism"], tier_sort_key(r["tier"]), r["region_class"], r["mapper"]),
    )

    def fmt_time(v):
        return "-" if v is None else f"{v:.1f}"

    def fmt_mem(v):
        return "-" if v is None else f"{v:.0f}"

    def fmt_int(v):
        return "-" if v is None else f"{int(v)}"

    body_rows = []
    for r in s:
        cls, glyph = f1_bucket(r["F1"])
        f1_txt = fmt_f1(r["F1"])
        org = escape(r["organism"])
        body_rows.append(
            "<tr>"
            f"<td><a href='organism/{org}.html'>{org}</a></td>"
            f"<td>{escape(r['tier'])}</td>"
            f"<td>{escape(r['mapper'])}</td>"
            f"<td>{escape(r['region_class'])}</td>"
            f"<td class='cell-{cls}'>{f1_txt}</td>"
            f"<td class='num'>{fmt_f1(r['precision'])}</td>"
            f"<td class='num'>{fmt_f1(r['recall'])}</td>"
            f"<td class='num'>{fmt_time(r['time_sec'])}</td>"
            f"<td class='num'>{fmt_mem(r['peak_mem_mb'])}</td>"
            f"<td class='num'>{fmt_int(r['n_reads'])}</td>"
            "</tr>"
        )
    table = f"""
<table id="scoreboard-table">
<thead><tr>
  <th class="sortable" data-col="0">organism</th>
  <th class="sortable" data-col="1">tier</th>
  <th class="sortable" data-col="2">mapper</th>
  <th class="sortable" data-col="3">region_class</th>
  <th class="sortable num" data-col="4">F1</th>
  <th class="sortable num" data-col="5">precision</th>
  <th class="sortable num" data-col="6">recall</th>
  <th class="sortable num" data-col="7">time (s)</th>
  <th class="sortable num" data-col="8">peak RSS (MB)</th>
  <th class="sortable num" data-col="9">n reads</th>
</tr></thead>
<tbody>
{"".join(body_rows)}
</tbody>
</table>
"""
    controls = """
<div class="controls">
  <label>Filter: <input id="filter" placeholder="organism / tier / mapper / region_class" size="48" /></label>
  <span class="muted">&nbsp;&middot;&nbsp; Click column headers to sort. Source: <a href="data/scoreboard.tsv">scoreboard.tsv</a>, <a href="data/scoreboard.json">scoreboard.json</a>.</span>
</div>
<script>
  // Lightweight sort + filter (no framework).
  (function() {
    var table = document.getElementById('scoreboard-table');
    if (!table) return;
    var tbody = table.tBodies[0];
    var ths = table.querySelectorAll('th.sortable');
    var sortDir = {};
    ths.forEach(function(th) {
      th.addEventListener('click', function() {
        var col = parseInt(th.getAttribute('data-col'), 10);
        var dir = sortDir[col] === 1 ? -1 : 1;
        sortDir[col] = dir;
        var numeric = th.classList.contains('num');
        var rows = Array.from(tbody.rows);
        rows.sort(function(a, b) {
          var av = a.cells[col].textContent.trim();
          var bv = b.cells[col].textContent.trim();
          if (numeric) {
            var ax = parseFloat(av); var bx = parseFloat(bv);
            if (isNaN(ax)) ax = -Infinity;
            if (isNaN(bx)) bx = -Infinity;
            return (ax - bx) * dir;
          }
          return av.localeCompare(bv) * dir;
        });
        rows.forEach(function(r) { tbody.appendChild(r); });
      });
    });
    var filt = document.getElementById('filter');
    if (filt) {
      filt.addEventListener('input', function() {
        var q = filt.value.toLowerCase();
        Array.from(tbody.rows).forEach(function(r) {
          r.style.display = r.textContent.toLowerCase().indexOf(q) === -1 ? 'none' : '';
        });
      });
    }
  })();
</script>
"""
    body = f"""
<section>
  <h2>Full scoreboard</h2>
  <p class="muted">All {len(rows)} (organism, tier, mapper, region_class) result rows. Color-coded by F1.</p>
  {controls}
  {table}
</section>
"""
    return page_shell("Scoreboard", body, active="scoreboard")


# --------------------------------------------------------------------------- per-organism
def load_timing(org: str, tier: str, mapper: str) -> dict | None:
    p = REPORTS / org / tier / mapper / "rep0" / "timing.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def render_organism(
    org: str,
    rows: list[dict],
    failed: list[dict],
) -> str:
    org_rows = [r for r in rows if r["organism"] == org]
    if not org_rows:
        body = f"<section><h2>{escape(org)}</h2><p>No data.</p></section>"
        return page_shell(org, body)

    tiers = sorted(set(r["tier"] for r in org_rows), key=tier_sort_key)
    region_classes = sorted(set(r["region_class"] for r in org_rows))
    mappers = sorted(set(r["mapper"] for r in org_rows))

    # Mapper-vs-mapper F1 per tier per region-class.
    # Build a table per tier.
    tier_sections = []
    for t in tiers:
        # rows for this tier
        # Map (region_class, mapper) -> F1.
        cell = {}
        for r in org_rows:
            if r["tier"] != t:
                continue
            cell[(r["region_class"], r["mapper"])] = r
        # header
        head = ["<table><thead><tr><th>region_class</th>"]
        for m in mappers:
            head.append(f"<th class='num'>{escape(m)} F1</th>")
        head.append("<th class='num'>best non-LLmap</th><th class='num'>LLmap delta</th></tr></thead><tbody>")
        rc_sorted = sorted(region_classes, key=lambda x: (0 if x == "all" else 1, x))
        for rc in rc_sorted:
            row_cells = [f"<td><strong>{escape(rc)}</strong></td>"]
            llmap_f1 = None
            best_other = None
            best_other_mapper = None
            for m in mappers:
                rec = cell.get((rc, m))
                v = rec["F1"] if rec else None
                cls, _ = f1_bucket(v)
                row_cells.append(f"<td class='cell-{cls}'>{fmt_f1(v)}</td>")
                if m == "llmap":
                    llmap_f1 = v
                else:
                    if v is not None and (best_other is None or v > best_other):
                        best_other = v
                        best_other_mapper = m
            if best_other is not None:
                row_cells.append(f"<td class='num'>{fmt_f1(best_other)} ({escape(best_other_mapper or '')})</td>")
            else:
                row_cells.append("<td class='num'>&mdash;</td>")
            if llmap_f1 is not None and best_other is not None:
                delta = llmap_f1 - best_other
                cls = "delta-pos" if delta > 0 else ("delta-neg" if delta < 0 else "")
                row_cells.append(f"<td class='num {cls}'>{delta:+.3f}</td>")
            else:
                row_cells.append("<td class='num'>&mdash;</td>")
            head.append("<tr>" + "".join(row_cells) + "</tr>")
        head.append("</tbody></table>")

        # Speed/memory comparison from timing.json
        timing_rows = []
        for m in mappers:
            tm = load_timing(org, t, m)
            if tm is None:
                timing_rows.append(
                    f"<tr><td>{escape(m)}</td>"
                    "<td class='num'>&mdash;</td><td class='num'>&mdash;</td>"
                    "<td class='num'>&mdash;</td><td class='num'>&mdash;</td></tr>"
                )
                continue
            wall = tm.get("wallclock_seconds")
            user = tm.get("user_cpu_seconds")
            mem = tm.get("peak_rss_bytes")
            mem_mb = (mem / (1024 * 1024)) if mem is not None else None
            bam = tm.get("bam_bytes")
            bam_kb = (bam / 1024) if bam is not None else None
            wall_s = "-" if wall is None else f"{wall:.1f}"
            user_s = "-" if user is None else f"{user:.1f}"
            mem_s = "-" if mem_mb is None else f"{mem_mb:.0f}"
            bam_s = "-" if bam_kb is None else f"{bam_kb:.0f}"
            timing_rows.append(
                f"<tr><td>{escape(m)}</td>"
                f"<td class='num'>{wall_s}</td>"
                f"<td class='num'>{user_s}</td>"
                f"<td class='num'>{mem_s}</td>"
                f"<td class='num'>{bam_s}</td></tr>"
            )
        timing_table = (
            "<h4>Speed &amp; memory</h4>"
            "<table><thead><tr>"
            "<th>mapper</th><th class='num'>wall (s)</th><th class='num'>user CPU (s)</th>"
            "<th class='num'>peak RSS (MB)</th><th class='num'>BAM (KB)</th>"
            "</tr></thead><tbody>"
            + "".join(timing_rows) +
            "</tbody></table>"
        )

        # Failures for this (org, tier).
        tf = [f for f in failed if f["organism"] == org and f["tier"] == t]
        fail_html = ""
        if tf:
            agg = defaultdict(list)
            for f in tf:
                agg[f["failure_mode"]].append(f"{f['mapper']}/{f['region_class']}")
            items = []
            for mode, lst in sorted(agg.items()):
                items.append(f"<li><code>{escape(mode)}</code>: {escape(', '.join(sorted(set(lst))))}</li>")
            fail_html = "<h4>Failure breakdown</h4><ul class='compact'>" + "".join(items) + "</ul>"

        tier_sections.append(
            f"<h3>{escape(t)}</h3>"
            + "".join(head)
            + timing_table
            + fail_html
        )

    body = f"""
<section>
  <p><a href="index.html">&larr; back to overview</a> &middot; <a href="../scoreboard.html">scoreboard</a></p>
  <h2>{escape(org)}</h2>
  <p class="muted">Per-tier, per-region-class breakdown. Mappers: {escape(', '.join(mappers))}.</p>
  {"".join(tier_sections)}
</section>
"""
    return page_shell(org, body, css_href="../" + CSS_HREF)


# --------------------------------------------------------------------------- methodology
def render_methodology() -> str:
    body = """
<section>
  <h2>Methodology</h2>

  <h3>Bench harness</h3>
  <p>We run a 4-mapper cross-species shootout on synthetic reads generated by
     <code>benchmarks/scripts/gen_synth_reads.py</code> with controlled error profiles
     and ground-truth alignment positions. Mappers under test:</p>
  <ul class="compact">
    <li><strong>LLmap</strong> &mdash; wave-particle adaptive &lambda;, GPU + AI optional.</li>
    <li><strong>minimap2</strong> &mdash; baseline long-read mapper (seed-chain-extend).</li>
    <li><strong>bwa-mem2</strong> &mdash; short-read baseline (mostly skipped on long-read tiers).</li>
    <li><strong>winnowmap</strong> &mdash; repeat-aware long-read mapper (run where applicable).</li>
  </ul>

  <h3>Tier ladder (gen_synth_reads.py tiers 1&ndash;10)</h3>
  <p>Each tier intensifies a stress dimension. Tiers 1, 2, 5, 6, 10 are the standard
     scoreboard checkpoints. <code>T1</code>/<code>T2</code> are the human-genome
     deep-dive cuts.</p>
  <ul class="compact">
    <li><strong>tier1</strong> &mdash; clean HiFi-like reads, low error, unique regions emphasized.</li>
    <li><strong>tier2</strong> &mdash; mixed error profile, balanced region coverage.</li>
    <li><strong>tier5</strong> &mdash; elevated indel rate, repetitive flanks.</li>
    <li><strong>tier6</strong> &mdash; structural variants (DEL/DUP/INV) layered on tier5 reads.</li>
    <li><strong>tier10</strong> &mdash; worst case: SVs + high error + telomere/centromere mass.</li>
    <li><strong>T1 / T2</strong> &mdash; human-genome-specific tiers (whole-genome, real chromosome lengths).</li>
  </ul>

  <h3>Region classes</h3>
  <p>Each truth read is annotated with the region class of its origin coordinate, so
     downstream metrics can be sliced by hardness band:</p>
  <ul class="compact">
    <li><code>unique</code> &mdash; single-copy reference regions.</li>
    <li><code>SD</code> &mdash; segmental duplications.</li>
    <li><code>centromere</code> &mdash; high-order repeats, alpha satellite.</li>
    <li><code>telomere</code> &mdash; terminal repeat blocks.</li>
    <li><code>low_complexity</code> &mdash; short-period tandem repeats.</li>
    <li><code>SV_DEL</code> / <code>SV_DUP</code> / <code>SV_INV</code> &mdash; reads spanning synthetic structural variants (tier6+).</li>
    <li><code>all</code> &mdash; aggregate across the whole truth set for that tier.</li>
  </ul>

  <h3>Metrics</h3>
  <p>For every (mapper, region_class) pair we record:</p>
  <ul class="compact">
    <li><strong>F1@100bp</strong> &mdash; harmonic mean of precision/recall where a read is
        "correct" if its predicted alignment is within 100&nbsp;bp of the truth coordinate
        (strand-aware). This is the strict F1 used for color-coding on the scoreboard.</li>
    <li><strong>F1@1kb</strong> &mdash; relaxed tolerance, useful for SV reads where the
        truth coordinate is only well-defined up to breakpoint resolution.</li>
    <li><strong>mapping rate</strong>, <strong>mean MAPQ</strong>, <strong>peak RSS</strong>, <strong>wallclock</strong>.</li>
  </ul>

  <h3>Color buckets</h3>
  <div class="legend">
    <span class="cell-excellent">&#10003; F1 &ge; 0.99</span>
    <span class="cell-good">&#10003; F1 &ge; 0.90</span>
    <span class="cell-warn">&#9888; F1 &ge; 0.70</span>
    <span class="cell-bad">&#10007; F1 &lt; 0.70</span>
    <span class="cell-na">&mdash; NA (mapper skipped)</span>
  </div>

  <h3>Reproducing the bench</h3>
  <p>From the repo root:</p>
  <pre>cd benchmarks
bash run_local_synthetic.sh         # local quick smoke
python3 orchestrate.py              # full local matrix
sbatch the HPC cluster_submit_t3_t6.sh       # SLURM heavy tiers (the HPC cluster)
python3 report.py                   # rebuild scoreboard.tsv + summary.md
python3 ../web/losslessmap.com/build.py   # rebuild this site</pre>

  <h3>Self-hosting / deployment</h3>
  <p>To self-host the scoreboard:</p>
  <pre>cd web/losslessmap.com &amp;&amp; python3 -m http.server 8080</pre>
  <p>For <code>losslessmap.com</code> deployment, either</p>
  <ul class="compact">
    <li><code>rsync -av web/losslessmap.com/ user@losslessmap.com:/var/www/losslessmap.com/</code>, or</li>
    <li>push <code>web/losslessmap.com/</code> to a <code>gh-pages</code> branch and point the
        domain's CNAME at the GitHub Pages site.</li>
  </ul>

  <h3>Data sources</h3>
  <ul class="compact">
    <li><a href="data/scoreboard.tsv">data/scoreboard.tsv</a> &mdash; long-format aggregate.</li>
    <li><a href="data/scoreboard.json">data/scoreboard.json</a> &mdash; same data as JSON for JS interactivity.</li>
    <li><code>benchmarks/scoreboard_pivoted.tsv</code> &mdash; pivoted view (organism &times; tier &times; region_class &rarr; mapper F1).</li>
    <li><code>benchmarks/reports/&lt;organism&gt;/&lt;tier&gt;/summary.json</code> &mdash; per-cell deep-dive.</li>
  </ul>
</section>
"""
    return page_shell("Methodology", body, active="methodology")


# --------------------------------------------------------------------------- main
def main() -> int:
    if not SCOREBOARD_TSV.exists():
        print(f"ERROR: {SCOREBOARD_TSV} not found", file=sys.stderr)
        return 1

    DATA_OUT.mkdir(parents=True, exist_ok=True)
    ORG_OUT.mkdir(parents=True, exist_ok=True)

    rows = load_scoreboard()
    summary_md_text = SUMMARY_MD.read_text() if SUMMARY_MD.exists() else ""
    summary = parse_summary_md(summary_md_text)

    failed: list[dict] = []
    if FAILED_JSON.exists():
        try:
            failed = json.loads(FAILED_JSON.read_text())
        except Exception as e:
            print(f"WARN: failed to parse {FAILED_JSON}: {e}", file=sys.stderr)
            failed = []

    organisms = sorted(set(r["organism"] for r in rows))
    tiers = sorted(set(r["tier"] for r in rows), key=tier_sort_key)

    # Copy raw data.
    shutil.copy(SCOREBOARD_TSV, DATA_OUT / "scoreboard.tsv")
    (DATA_OUT / "scoreboard.json").write_text(json.dumps(rows, indent=2))

    # CSS.
    (ROOT / "styles.css").write_text(STYLES_CSS)

    # Pages.
    (ROOT / "index.html").write_text(render_index(rows, summary, organisms, tiers))
    (ROOT / "scoreboard.html").write_text(render_scoreboard(rows))
    (ROOT / "methodology.html").write_text(render_methodology())

    for org in organisms:
        (ORG_OUT / f"{org}.html").write_text(render_organism(org, rows, failed))

    print(f"[build.py] wrote {len(organisms)} organism pages, {len(rows)} scoreboard rows.")
    print(f"[build.py] output: {ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
