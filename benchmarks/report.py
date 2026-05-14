#!/usr/bin/env python3
"""Report generator for LLmap benchmark campaign.

Reads benchmark results from benchmarks/reports/<task>/<tool>/rep<N>/
and generates:
    - Per-task README.md summaries
    - Cross-tool comparison tables (TSV/Markdown)
    - Plot generation stubs (PNG if matplotlib available)

Usage:
    python report.py                         # generate all reports
    python report.py --tasks T1 T2           # specific tasks only
    python report.py --output-dir ./results  # custom output location
    python report.py --format markdown       # markdown only (no TSV)
    python report.py --no-plots              # skip plot generation
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


@dataclass
class RunResult:
    """Results from a single benchmark run."""
    task: str
    tool: str
    replicate: int
    mapping_summary: Optional[dict] = None
    mapq_histogram: Optional[dict] = None
    ground_truth: Optional[dict] = None
    resources: Optional[dict] = None
    manifest: Optional[dict] = None


@dataclass
class ToolSummary:
    """Aggregated results for a tool across replicates."""
    tool: str
    task: str
    replicates: int = 0
    primary_mapped_mean: float = 0.0
    primary_mapped_std: float = 0.0
    mapping_rate_mean: float = 0.0
    mapping_rate_std: float = 0.0
    recall_mean: float = 0.0
    recall_std: float = 0.0
    precision_mean: float = 0.0
    precision_std: float = 0.0
    f1_mean: float = 0.0
    f1_std: float = 0.0
    wallclock_mean: float = 0.0
    wallclock_std: float = 0.0
    peak_rss_mb_mean: float = 0.0
    peak_rss_mb_std: float = 0.0


@dataclass
class TaskReport:
    """Complete report for a task."""
    task: str
    tools: list = field(default_factory=list)
    summaries: dict = field(default_factory=dict)  # tool -> ToolSummary
    concordance: dict = field(default_factory=dict)  # (tool_a, tool_b) -> dict
    runs: list = field(default_factory=list)  # list of RunResult


@dataclass
class ReportConfig:
    """Report generator configuration."""
    llmap_root: Path
    output_dir: Optional[Path] = None
    tasks_filter: list = field(default_factory=list)
    format: str = "all"  # "all", "markdown", "tsv"
    generate_plots: bool = True
    verbose: bool = False


def load_json_safe(path: Path) -> Optional[dict]:
    """Load JSON file, return None if missing or invalid."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def load_run_result(run_dir: Path, task: str, tool: str, rep: int) -> Optional[RunResult]:
    """Load all result files from a single run directory."""
    if not run_dir.exists():
        return None

    done_flag = run_dir / "done.flag"
    if not done_flag.exists():
        return None

    result = RunResult(task=task, tool=tool, replicate=rep)
    result.mapping_summary = load_json_safe(run_dir / "mapping_summary.json")
    result.mapq_histogram = load_json_safe(run_dir / "mapq_histogram.json")
    result.ground_truth = load_json_safe(run_dir / "ground_truth.json")
    result.resources = load_json_safe(run_dir / "resources.json")
    result.manifest = load_json_safe(run_dir / "manifest.json")

    if result.mapping_summary is None:
        return None

    return result


def aggregate_tool(runs: list) -> ToolSummary:
    """Aggregate multiple replicate runs into a summary."""
    if not runs:
        return ToolSummary(tool="", task="")

    first = runs[0]
    summary = ToolSummary(tool=first.tool, task=first.task, replicates=len(runs))

    primary_mapped = []
    mapping_rates = []
    recalls = []
    precisions = []
    f1s = []
    wallclocks = []
    peak_rss = []

    for r in runs:
        if r.mapping_summary:
            primary_mapped.append(r.mapping_summary.get("primary_mapped", 0))
            mapping_rates.append(r.mapping_summary.get("mapping_rate", 0.0))

        if r.ground_truth:
            recalls.append(r.ground_truth.get("recall", 0.0))
            precisions.append(r.ground_truth.get("precision", 0.0))
            f1s.append(r.ground_truth.get("f1", 0.0))

        if r.resources:
            wallclocks.append(r.resources.get("wallclock_seconds", 0.0))
            rss_bytes = r.resources.get("peak_rss_bytes", 0)
            peak_rss.append(rss_bytes / (1024 * 1024))

    def safe_mean(vals):
        return statistics.mean(vals) if vals else 0.0

    def safe_stdev(vals):
        return statistics.stdev(vals) if len(vals) > 1 else 0.0

    summary.primary_mapped_mean = safe_mean(primary_mapped)
    summary.primary_mapped_std = safe_stdev(primary_mapped)
    summary.mapping_rate_mean = safe_mean(mapping_rates)
    summary.mapping_rate_std = safe_stdev(mapping_rates)
    summary.recall_mean = safe_mean(recalls)
    summary.recall_std = safe_stdev(recalls)
    summary.precision_mean = safe_mean(precisions)
    summary.precision_std = safe_stdev(precisions)
    summary.f1_mean = safe_mean(f1s)
    summary.f1_std = safe_stdev(f1s)
    summary.wallclock_mean = safe_mean(wallclocks)
    summary.wallclock_std = safe_stdev(wallclocks)
    summary.peak_rss_mb_mean = safe_mean(peak_rss)
    summary.peak_rss_mb_std = safe_stdev(peak_rss)

    return summary


def load_task_report(reports_dir: Path, task: str) -> TaskReport:
    """Load all results for a task and aggregate by tool."""
    task_dir = reports_dir / task
    report = TaskReport(task=task)

    if not task_dir.exists():
        return report

    for tool_dir in sorted(task_dir.iterdir()):
        if not tool_dir.is_dir():
            continue
        tool = tool_dir.name
        report.tools.append(tool)

        runs = []
        for rep_dir in sorted(tool_dir.iterdir()):
            if not rep_dir.is_dir() or not rep_dir.name.startswith("rep"):
                continue
            try:
                rep = int(rep_dir.name[3:])
            except ValueError:
                continue

            run = load_run_result(rep_dir, task, tool, rep)
            if run:
                runs.append(run)
                report.runs.append(run)

        if runs:
            report.summaries[tool] = aggregate_tool(runs)

    concordance_file = task_dir / "concordance_matrix.json"
    if concordance_file.exists():
        report.concordance = load_json_safe(concordance_file) or {}

    return report


def format_number(val: float, decimals: int = 2) -> str:
    """Format a number for display."""
    if abs(val) >= 1000000:
        return f"{val/1000000:.1f}M"
    if abs(val) >= 1000:
        return f"{val/1000:.1f}K"
    return f"{val:.{decimals}f}"


def format_pct(val: float) -> str:
    """Format a percentage (0-1 scale) for display."""
    return f"{val*100:.1f}%"


def format_mean_std(mean: float, std: float, is_pct: bool = False) -> str:
    """Format mean +/- std for display."""
    if is_pct:
        return f"{mean*100:.1f} +/- {std*100:.1f}%"
    return f"{format_number(mean)} +/- {format_number(std)}"


def generate_task_readme(report: TaskReport, output_path: Path) -> None:
    """Generate a README.md for a task."""
    lines = [
        f"# Benchmark Results: Task {report.task}",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
    ]

    task_descriptions = {
        "T1": "Synthetic-truth WGS (1M simulated HiFi reads)",
        "T2": "Synthetic paralog stress (500k reads from duplicated regions)",
        "T3": "Real HiFi WGS (HG002 chr14+chr20 subset)",
        "T4": "Real Illumina WGS (HG002 chr14+chr20 subset)",
        "T5": "Iso-seq transcriptome (HPRC lymph FLNC)",
        "T6": "Targeted IGH-locus (paralog disambiguation challenge)",
    }

    desc = task_descriptions.get(report.task, "Unknown task")
    lines.extend([
        f"## Task Description",
        "",
        desc,
        "",
    ])

    if not report.summaries:
        lines.extend([
            "## Results",
            "",
            "No completed benchmark runs found.",
            "",
        ])
        output_path.write_text("\n".join(lines))
        return

    lines.extend([
        "## Summary Table",
        "",
        "| Tool | Replicates | Mapping Rate | Recall | Precision | F1 | Wallclock (s) | Peak RSS (MB) |",
        "|------|------------|--------------|--------|-----------|-----|---------------|---------------|",
    ])

    for tool in sorted(report.summaries.keys()):
        s = report.summaries[tool]
        row = (
            f"| {tool} "
            f"| {s.replicates} "
            f"| {format_pct(s.mapping_rate_mean)} "
            f"| {format_pct(s.recall_mean) if s.recall_mean else 'N/A'} "
            f"| {format_pct(s.precision_mean) if s.precision_mean else 'N/A'} "
            f"| {format_pct(s.f1_mean) if s.f1_mean else 'N/A'} "
            f"| {format_number(s.wallclock_mean, 1)} "
            f"| {format_number(s.peak_rss_mb_mean, 0)} |"
        )
        lines.append(row)

    lines.append("")

    lines.extend([
        "## Detailed Statistics",
        "",
    ])

    for tool in sorted(report.summaries.keys()):
        s = report.summaries[tool]
        lines.extend([
            f"### {tool}",
            "",
            f"- **Replicates**: {s.replicates}",
            f"- **Primary mapped**: {format_mean_std(s.primary_mapped_mean, s.primary_mapped_std)}",
            f"- **Mapping rate**: {format_mean_std(s.mapping_rate_mean, s.mapping_rate_std, is_pct=True)}",
        ])

        if s.recall_mean > 0:
            lines.extend([
                f"- **Recall**: {format_mean_std(s.recall_mean, s.recall_std, is_pct=True)}",
                f"- **Precision**: {format_mean_std(s.precision_mean, s.precision_std, is_pct=True)}",
                f"- **F1**: {format_mean_std(s.f1_mean, s.f1_std, is_pct=True)}",
            ])

        if s.wallclock_mean > 0:
            lines.extend([
                f"- **Wallclock**: {format_mean_std(s.wallclock_mean, s.wallclock_std)} seconds",
                f"- **Peak RSS**: {format_mean_std(s.peak_rss_mb_mean, s.peak_rss_mb_std)} MB",
            ])

        lines.append("")

    if report.concordance:
        lines.extend([
            "## Pairwise Concordance",
            "",
        ])

        for key, val in report.concordance.items():
            if isinstance(val, dict) and "name_a" in val and "name_b" in val:
                lines.extend([
                    f"### {val['name_a']} vs {val['name_b']}",
                    "",
                ])
                if "percentages" in val:
                    pcts = val["percentages"]
                    lines.extend([
                        f"- Same position: {pcts.get('same_position', 0):.1f}%",
                        f"- Same chrom, diff pos: {pcts.get('same_chrom_diff_pos', 0):.1f}%",
                        f"- Different chrom: {pcts.get('diff_chrom', 0):.1f}%",
                        f"- Only in A: {pcts.get('oneA', 0):.1f}%",
                        f"- Only in B: {pcts.get('oneB', 0):.1f}%",
                        f"- Both unmapped: {pcts.get('both_unmapped', 0):.1f}%",
                        "",
                    ])

    output_path.write_text("\n".join(lines))


def generate_comparison_tsv(reports: list, output_path: Path) -> None:
    """Generate a cross-tool comparison TSV."""
    headers = [
        "task", "tool", "replicates",
        "mapping_rate_mean", "mapping_rate_std",
        "recall_mean", "recall_std",
        "precision_mean", "precision_std",
        "f1_mean", "f1_std",
        "wallclock_mean", "wallclock_std",
        "peak_rss_mb_mean", "peak_rss_mb_std",
    ]

    lines = ["\t".join(headers)]

    for report in reports:
        for tool in sorted(report.summaries.keys()):
            s = report.summaries[tool]
            row = [
                report.task,
                tool,
                str(s.replicates),
                f"{s.mapping_rate_mean:.4f}",
                f"{s.mapping_rate_std:.4f}",
                f"{s.recall_mean:.4f}",
                f"{s.recall_std:.4f}",
                f"{s.precision_mean:.4f}",
                f"{s.precision_std:.4f}",
                f"{s.f1_mean:.4f}",
                f"{s.f1_std:.4f}",
                f"{s.wallclock_mean:.2f}",
                f"{s.wallclock_std:.2f}",
                f"{s.peak_rss_mb_mean:.1f}",
                f"{s.peak_rss_mb_std:.1f}",
            ]
            lines.append("\t".join(row))

    output_path.write_text("\n".join(lines) + "\n")


def generate_comparison_markdown(reports: list, output_path: Path) -> None:
    """Generate a cross-tool comparison Markdown table."""
    lines = [
        "# LLmap Benchmark Campaign: Cross-Tool Comparison",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Results Summary",
        "",
        "| Task | Tool | Mapping Rate | Recall | F1 | Wallclock (s) |",
        "|------|------|--------------|--------|-----|---------------|",
    ]

    for report in reports:
        for tool in sorted(report.summaries.keys()):
            s = report.summaries[tool]
            recall = format_pct(s.recall_mean) if s.recall_mean > 0 else "N/A"
            f1 = format_pct(s.f1_mean) if s.f1_mean > 0 else "N/A"
            wallclock = format_number(s.wallclock_mean, 1) if s.wallclock_mean > 0 else "N/A"
            row = f"| {report.task} | {tool} | {format_pct(s.mapping_rate_mean)} | {recall} | {f1} | {wallclock} |"
            lines.append(row)

    lines.append("")
    output_path.write_text("\n".join(lines))


def generate_plots(reports: list, output_dir: Path) -> bool:
    """Generate benchmark comparison plots."""
    if not HAS_MATPLOTLIB:
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    tasks_with_data = [r for r in reports if r.summaries]
    if not tasks_with_data:
        return False

    generate_mapping_rate_plot(tasks_with_data, output_dir / "mapping_rate.png")
    generate_f1_plot(tasks_with_data, output_dir / "f1_score.png")
    generate_performance_plot(tasks_with_data, output_dir / "wallclock.png")

    return True


def generate_mapping_rate_plot(reports: list, output_path: Path) -> None:
    """Generate mapping rate comparison bar chart."""
    fig, ax = plt.subplots(figsize=(10, 6))

    all_tools = sorted(set(t for r in reports for t in r.summaries.keys()))
    x_positions = range(len(reports))
    bar_width = 0.8 / len(all_tools) if all_tools else 0.8

    for i, tool in enumerate(all_tools):
        values = []
        errors = []
        for report in reports:
            if tool in report.summaries:
                s = report.summaries[tool]
                values.append(s.mapping_rate_mean * 100)
                errors.append(s.mapping_rate_std * 100)
            else:
                values.append(0)
                errors.append(0)

        offset = (i - len(all_tools) / 2 + 0.5) * bar_width
        positions = [x + offset for x in x_positions]
        ax.bar(positions, values, bar_width, yerr=errors, label=tool, capsize=3)

    ax.set_xlabel('Task')
    ax.set_ylabel('Mapping Rate (%)')
    ax.set_title('Mapping Rate by Task and Tool')
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels([r.task for r in reports])
    ax.legend(loc='lower right')
    ax.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def generate_f1_plot(reports: list, output_path: Path) -> None:
    """Generate F1 score comparison bar chart (only for tasks with ground truth)."""
    gt_reports = [r for r in reports if any(s.f1_mean > 0 for s in r.summaries.values())]
    if not gt_reports:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    all_tools = sorted(set(t for r in gt_reports for t in r.summaries.keys()))
    x_positions = range(len(gt_reports))
    bar_width = 0.8 / len(all_tools) if all_tools else 0.8

    for i, tool in enumerate(all_tools):
        values = []
        errors = []
        for report in gt_reports:
            if tool in report.summaries:
                s = report.summaries[tool]
                values.append(s.f1_mean * 100)
                errors.append(s.f1_std * 100)
            else:
                values.append(0)
                errors.append(0)

        offset = (i - len(all_tools) / 2 + 0.5) * bar_width
        positions = [x + offset for x in x_positions]
        ax.bar(positions, values, bar_width, yerr=errors, label=tool, capsize=3)

    ax.set_xlabel('Task')
    ax.set_ylabel('F1 Score (%)')
    ax.set_title('F1 Score by Task and Tool (Synthetic Tasks with Ground Truth)')
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels([r.task for r in gt_reports])
    ax.legend(loc='lower right')
    ax.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def generate_performance_plot(reports: list, output_path: Path) -> None:
    """Generate wallclock time comparison bar chart."""
    perf_reports = [r for r in reports if any(s.wallclock_mean > 0 for s in r.summaries.values())]
    if not perf_reports:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    all_tools = sorted(set(t for r in perf_reports for t in r.summaries.keys()))
    x_positions = range(len(perf_reports))
    bar_width = 0.8 / len(all_tools) if all_tools else 0.8

    for i, tool in enumerate(all_tools):
        values = []
        errors = []
        for report in perf_reports:
            if tool in report.summaries:
                s = report.summaries[tool]
                values.append(s.wallclock_mean)
                errors.append(s.wallclock_std)
            else:
                values.append(0)
                errors.append(0)

        offset = (i - len(all_tools) / 2 + 0.5) * bar_width
        positions = [x + offset for x in x_positions]
        ax.bar(positions, values, bar_width, yerr=errors, label=tool, capsize=3)

    ax.set_xlabel('Task')
    ax.set_ylabel('Wallclock Time (seconds)')
    ax.set_title('Execution Time by Task and Tool')
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels([r.task for r in perf_reports])
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def generate_reports(config: ReportConfig) -> dict:
    """Main report generation entry point."""
    reports_dir = config.llmap_root / "benchmarks" / "reports"
    output_dir = config.output_dir or reports_dir

    results = {
        "timestamp": datetime.now().isoformat(),
        "tasks_processed": [],
        "readmes_generated": [],
        "tables_generated": [],
        "plots_generated": [],
        "errors": [],
    }

    all_tasks = ["T1", "T2", "T3", "T4", "T5", "T6"]
    if config.tasks_filter:
        all_tasks = [t for t in all_tasks if t in config.tasks_filter]

    all_reports = []

    for task in all_tasks:
        report = load_task_report(reports_dir, task)
        all_reports.append(report)
        results["tasks_processed"].append(task)

        if report.summaries:
            task_output_dir = output_dir / task
            task_output_dir.mkdir(parents=True, exist_ok=True)

            readme_path = task_output_dir / "README.md"
            generate_task_readme(report, readme_path)
            results["readmes_generated"].append(str(readme_path))

            if config.verbose:
                print(f"Generated: {readme_path}")

    reports_with_data = [r for r in all_reports if r.summaries]

    if reports_with_data:
        output_dir.mkdir(parents=True, exist_ok=True)

        if config.format in ("all", "tsv"):
            tsv_path = output_dir / "comparison.tsv"
            generate_comparison_tsv(reports_with_data, tsv_path)
            results["tables_generated"].append(str(tsv_path))
            if config.verbose:
                print(f"Generated: {tsv_path}")

        if config.format in ("all", "markdown"):
            md_path = output_dir / "comparison.md"
            generate_comparison_markdown(reports_with_data, md_path)
            results["tables_generated"].append(str(md_path))
            if config.verbose:
                print(f"Generated: {md_path}")

        if config.generate_plots and HAS_MATPLOTLIB:
            figures_dir = output_dir / "figures"
            if generate_plots(reports_with_data, figures_dir):
                results["plots_generated"].extend([
                    str(figures_dir / "mapping_rate.png"),
                    str(figures_dir / "f1_score.png"),
                    str(figures_dir / "wallclock.png"),
                ])
                if config.verbose:
                    print(f"Generated plots in: {figures_dir}")
        elif config.generate_plots and not HAS_MATPLOTLIB:
            results["errors"].append("matplotlib not installed, skipping plots")
            if config.verbose:
                print("Warning: matplotlib not installed, skipping plots")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Report generator for LLmap benchmark campaign"
    )
    parser.add_argument("--tasks", nargs="+", default=[],
                        help="Filter to specific tasks (e.g., T1 T2)")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory (default: benchmarks/reports)")
    parser.add_argument("--format", choices=["all", "markdown", "tsv"],
                        default="all", help="Output format (default: all)")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip plot generation")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print progress")
    parser.add_argument("--json", action="store_true",
                        help="Output results as JSON")
    args = parser.parse_args()

    llmap_root = Path(__file__).parent.parent.resolve()

    config = ReportConfig(
        llmap_root=llmap_root,
        output_dir=args.output_dir,
        tasks_filter=args.tasks,
        format=args.format,
        generate_plots=not args.no_plots,
        verbose=args.verbose,
    )

    results = generate_reports(config)

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print(f"\nReport generation summary:")
        print(f"  Tasks processed: {len(results['tasks_processed'])}")
        print(f"  READMEs generated: {len(results['readmes_generated'])}")
        print(f"  Tables generated: {len(results['tables_generated'])}")
        print(f"  Plots generated: {len(results['plots_generated'])}")
        if results['errors']:
            print(f"  Warnings: {len(results['errors'])}")
            for err in results['errors']:
                print(f"    - {err}")


if __name__ == "__main__":
    main()
