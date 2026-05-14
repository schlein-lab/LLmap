#!/usr/bin/env python3
"""Unit tests for benchmarks/report.py."""

import json
import tempfile
import unittest
from pathlib import Path

from report import (
    ReportConfig,
    RunResult,
    TaskReport,
    ToolSummary,
    aggregate_tool,
    format_mean_std,
    format_number,
    format_pct,
    generate_comparison_markdown,
    generate_comparison_tsv,
    generate_reports,
    generate_task_readme,
    load_json_safe,
    load_run_result,
    load_task_report,
)


class TestFormatHelpers(unittest.TestCase):
    """Tests for formatting helper functions."""

    def test_format_number_small(self):
        self.assertEqual(format_number(0.0), "0.00")
        self.assertEqual(format_number(0.12345), "0.12")
        self.assertEqual(format_number(42.5), "42.50")

    def test_format_number_thousands(self):
        self.assertEqual(format_number(1500), "1.5K")
        self.assertEqual(format_number(42000), "42.0K")

    def test_format_number_millions(self):
        self.assertEqual(format_number(1500000), "1.5M")
        self.assertEqual(format_number(42000000), "42.0M")

    def test_format_number_decimals(self):
        self.assertEqual(format_number(3.14159, decimals=3), "3.142")
        self.assertEqual(format_number(3.14159, decimals=0), "3")

    def test_format_pct(self):
        self.assertEqual(format_pct(0.0), "0.0%")
        self.assertEqual(format_pct(0.5), "50.0%")
        self.assertEqual(format_pct(1.0), "100.0%")
        self.assertEqual(format_pct(0.956), "95.6%")

    def test_format_mean_std(self):
        result = format_mean_std(100.0, 5.0)
        self.assertIn("100", result)
        self.assertIn("5", result)

    def test_format_mean_std_pct(self):
        result = format_mean_std(0.95, 0.02, is_pct=True)
        self.assertIn("95", result)
        self.assertIn("2", result)
        self.assertIn("%", result)


class TestLoadJsonSafe(unittest.TestCase):
    """Tests for load_json_safe function."""

    def test_load_valid_json(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump({"key": "value", "num": 42}, f)
            f.flush()
            result = load_json_safe(Path(f.name))

        self.assertEqual(result, {"key": "value", "num": 42})
        Path(f.name).unlink()

    def test_load_missing_file(self):
        result = load_json_safe(Path("/nonexistent/file.json"))
        self.assertIsNone(result)

    def test_load_invalid_json(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            f.write("not valid json {{{")
            f.flush()
            result = load_json_safe(Path(f.name))

        self.assertIsNone(result)
        Path(f.name).unlink()


class TestRunResult(unittest.TestCase):
    """Tests for RunResult dataclass."""

    def test_create_run_result(self):
        r = RunResult(task="T1", tool="llmap", replicate=0)
        self.assertEqual(r.task, "T1")
        self.assertEqual(r.tool, "llmap")
        self.assertEqual(r.replicate, 0)
        self.assertIsNone(r.mapping_summary)

    def test_run_result_with_data(self):
        r = RunResult(
            task="T1",
            tool="minimap2",
            replicate=1,
            mapping_summary={"primary_mapped": 1000, "mapping_rate": 0.95},
            ground_truth={"recall": 0.92, "precision": 0.94, "f1": 0.93},
        )
        self.assertEqual(r.mapping_summary["primary_mapped"], 1000)
        self.assertEqual(r.ground_truth["recall"], 0.92)


class TestLoadRunResult(unittest.TestCase):
    """Tests for load_run_result function."""

    def test_load_missing_dir(self):
        result = load_run_result(Path("/nonexistent"), "T1", "llmap", 0)
        self.assertIsNone(result)

    def test_load_no_done_flag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "run"
            run_dir.mkdir()
            result = load_run_result(run_dir, "T1", "llmap", 0)
            self.assertIsNone(result)

    def test_load_complete_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            (run_dir / "done.flag").touch()
            (run_dir / "mapping_summary.json").write_text(json.dumps({
                "primary_mapped": 950,
                "mapping_rate": 0.95,
            }))
            (run_dir / "ground_truth.json").write_text(json.dumps({
                "recall": 0.92,
                "precision": 0.94,
                "f1": 0.93,
            }))

            result = load_run_result(run_dir, "T1", "llmap", 0)

            self.assertIsNotNone(result)
            self.assertEqual(result.task, "T1")
            self.assertEqual(result.tool, "llmap")
            self.assertEqual(result.mapping_summary["primary_mapped"], 950)
            self.assertEqual(result.ground_truth["recall"], 0.92)


class TestAggregateRuns(unittest.TestCase):
    """Tests for aggregate_tool function."""

    def test_aggregate_empty(self):
        result = aggregate_tool([])
        self.assertEqual(result.tool, "")
        self.assertEqual(result.replicates, 0)

    def test_aggregate_single_run(self):
        runs = [
            RunResult(
                task="T1", tool="llmap", replicate=0,
                mapping_summary={"primary_mapped": 1000, "mapping_rate": 0.95},
                ground_truth={"recall": 0.92, "precision": 0.94, "f1": 0.93},
                resources={"wallclock_seconds": 120.5, "peak_rss_bytes": 1024*1024*500},
            ),
        ]

        result = aggregate_tool(runs)

        self.assertEqual(result.tool, "llmap")
        self.assertEqual(result.replicates, 1)
        self.assertAlmostEqual(result.mapping_rate_mean, 0.95)
        self.assertAlmostEqual(result.recall_mean, 0.92)
        self.assertAlmostEqual(result.wallclock_mean, 120.5)
        self.assertAlmostEqual(result.peak_rss_mb_mean, 500.0)

    def test_aggregate_multiple_runs(self):
        runs = [
            RunResult(
                task="T1", tool="minimap2", replicate=0,
                mapping_summary={"primary_mapped": 1000, "mapping_rate": 0.94},
                resources={"wallclock_seconds": 100.0, "peak_rss_bytes": 1024*1024*400},
            ),
            RunResult(
                task="T1", tool="minimap2", replicate=1,
                mapping_summary={"primary_mapped": 1020, "mapping_rate": 0.96},
                resources={"wallclock_seconds": 110.0, "peak_rss_bytes": 1024*1024*420},
            ),
            RunResult(
                task="T1", tool="minimap2", replicate=2,
                mapping_summary={"primary_mapped": 1010, "mapping_rate": 0.95},
                resources={"wallclock_seconds": 105.0, "peak_rss_bytes": 1024*1024*410},
            ),
        ]

        result = aggregate_tool(runs)

        self.assertEqual(result.replicates, 3)
        self.assertAlmostEqual(result.mapping_rate_mean, 0.95, places=3)
        self.assertAlmostEqual(result.wallclock_mean, 105.0, places=1)
        self.assertGreater(result.mapping_rate_std, 0)


class TestToolSummary(unittest.TestCase):
    """Tests for ToolSummary dataclass."""

    def test_default_values(self):
        s = ToolSummary(tool="llmap", task="T1")
        self.assertEqual(s.replicates, 0)
        self.assertEqual(s.mapping_rate_mean, 0.0)
        self.assertEqual(s.f1_mean, 0.0)


class TestTaskReport(unittest.TestCase):
    """Tests for TaskReport dataclass."""

    def test_empty_report(self):
        report = TaskReport(task="T1")
        self.assertEqual(report.task, "T1")
        self.assertEqual(report.tools, [])
        self.assertEqual(report.summaries, {})

    def test_report_with_data(self):
        report = TaskReport(
            task="T1",
            tools=["llmap", "minimap2"],
            summaries={
                "llmap": ToolSummary(tool="llmap", task="T1", replicates=3),
                "minimap2": ToolSummary(tool="minimap2", task="T1", replicates=3),
            },
        )
        self.assertEqual(len(report.tools), 2)
        self.assertIn("llmap", report.summaries)


class TestLoadTaskReport(unittest.TestCase):
    """Tests for load_task_report function."""

    def test_load_missing_task(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reports_dir = Path(tmpdir)
            result = load_task_report(reports_dir, "T1")
            self.assertEqual(result.task, "T1")
            self.assertEqual(result.tools, [])

    def test_load_task_with_runs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reports_dir = Path(tmpdir)
            task_dir = reports_dir / "T1" / "llmap" / "rep0"
            task_dir.mkdir(parents=True)

            (task_dir / "done.flag").touch()
            (task_dir / "mapping_summary.json").write_text(json.dumps({
                "primary_mapped": 950,
                "mapping_rate": 0.95,
            }))

            result = load_task_report(reports_dir, "T1")

            self.assertEqual(result.task, "T1")
            self.assertIn("llmap", result.tools)
            self.assertIn("llmap", result.summaries)
            self.assertEqual(result.summaries["llmap"].replicates, 1)


class TestGenerateTaskReadme(unittest.TestCase):
    """Tests for generate_task_readme function."""

    def test_generate_empty_readme(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report = TaskReport(task="T1")
            output = Path(tmpdir) / "README.md"

            generate_task_readme(report, output)

            content = output.read_text()
            self.assertIn("Task T1", content)
            self.assertIn("No completed benchmark runs", content)

    def test_generate_readme_with_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report = TaskReport(
                task="T1",
                tools=["llmap"],
                summaries={
                    "llmap": ToolSummary(
                        tool="llmap", task="T1", replicates=3,
                        mapping_rate_mean=0.95, mapping_rate_std=0.01,
                        recall_mean=0.92, recall_std=0.02,
                        f1_mean=0.93, f1_std=0.01,
                        wallclock_mean=120.0, wallclock_std=5.0,
                    ),
                },
            )
            output = Path(tmpdir) / "README.md"

            generate_task_readme(report, output)

            content = output.read_text()
            self.assertIn("Task T1", content)
            self.assertIn("llmap", content)
            self.assertIn("95", content)  # mapping rate
            self.assertIn("Summary Table", content)


class TestGenerateComparisonTsv(unittest.TestCase):
    """Tests for generate_comparison_tsv function."""

    def test_generate_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = [
                TaskReport(
                    task="T1",
                    summaries={
                        "llmap": ToolSummary(
                            tool="llmap", task="T1", replicates=3,
                            mapping_rate_mean=0.95, recall_mean=0.92,
                        ),
                    },
                ),
            ]
            output = Path(tmpdir) / "comparison.tsv"

            generate_comparison_tsv(reports, output)

            content = output.read_text()
            lines = content.strip().split("\n")
            self.assertEqual(len(lines), 2)  # header + 1 data row
            self.assertIn("task", lines[0])
            self.assertIn("T1", lines[1])
            self.assertIn("llmap", lines[1])


class TestGenerateComparisonMarkdown(unittest.TestCase):
    """Tests for generate_comparison_markdown function."""

    def test_generate_markdown(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = [
                TaskReport(
                    task="T1",
                    summaries={
                        "llmap": ToolSummary(
                            tool="llmap", task="T1", replicates=3,
                            mapping_rate_mean=0.95, f1_mean=0.93,
                        ),
                        "minimap2": ToolSummary(
                            tool="minimap2", task="T1", replicates=3,
                            mapping_rate_mean=0.94, f1_mean=0.91,
                        ),
                    },
                ),
            ]
            output = Path(tmpdir) / "comparison.md"

            generate_comparison_markdown(reports, output)

            content = output.read_text()
            self.assertIn("LLmap Benchmark Campaign", content)
            self.assertIn("| T1 | llmap", content)
            self.assertIn("| T1 | minimap2", content)


class TestReportConfig(unittest.TestCase):
    """Tests for ReportConfig dataclass."""

    def test_default_config(self):
        config = ReportConfig(llmap_root=Path("/tmp"))
        self.assertEqual(config.format, "all")
        self.assertTrue(config.generate_plots)
        self.assertFalse(config.verbose)

    def test_custom_config(self):
        config = ReportConfig(
            llmap_root=Path("/tmp"),
            tasks_filter=["T1", "T2"],
            format="markdown",
            generate_plots=False,
        )
        self.assertEqual(config.tasks_filter, ["T1", "T2"])
        self.assertEqual(config.format, "markdown")
        self.assertFalse(config.generate_plots)


class TestGenerateReports(unittest.TestCase):
    """Tests for generate_reports function."""

    def test_generate_empty_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            llmap_root = Path(tmpdir)
            (llmap_root / "benchmarks" / "reports").mkdir(parents=True)

            config = ReportConfig(
                llmap_root=llmap_root,
                generate_plots=False,
            )

            results = generate_reports(config)

            self.assertEqual(len(results["tasks_processed"]), 6)  # T1-T6
            self.assertEqual(len(results["readmes_generated"]), 0)  # no data

    def test_generate_with_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            llmap_root = Path(tmpdir)
            reports_dir = llmap_root / "benchmarks" / "reports"

            run_dir = reports_dir / "T1" / "llmap" / "rep0"
            run_dir.mkdir(parents=True)
            (run_dir / "done.flag").touch()
            (run_dir / "mapping_summary.json").write_text(json.dumps({
                "primary_mapped": 950,
                "mapping_rate": 0.95,
            }))

            config = ReportConfig(
                llmap_root=llmap_root,
                tasks_filter=["T1"],
                generate_plots=False,
            )

            results = generate_reports(config)

            self.assertEqual(len(results["tasks_processed"]), 1)
            self.assertEqual(len(results["readmes_generated"]), 1)
            self.assertEqual(len(results["tables_generated"]), 2)  # TSV + MD


class TestIntegration(unittest.TestCase):
    """Integration tests for full report generation workflow."""

    def test_full_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            llmap_root = Path(tmpdir)
            reports_dir = llmap_root / "benchmarks" / "reports"

            for tool in ["llmap", "minimap2"]:
                for rep in range(2):
                    run_dir = reports_dir / "T1" / tool / f"rep{rep}"
                    run_dir.mkdir(parents=True)
                    (run_dir / "done.flag").touch()
                    (run_dir / "mapping_summary.json").write_text(json.dumps({
                        "primary_mapped": 900 + rep * 50,
                        "mapping_rate": 0.90 + rep * 0.05,
                    }))
                    (run_dir / "ground_truth.json").write_text(json.dumps({
                        "recall": 0.88 + rep * 0.04,
                        "precision": 0.90 + rep * 0.03,
                        "f1": 0.89 + rep * 0.03,
                    }))

            config = ReportConfig(
                llmap_root=llmap_root,
                tasks_filter=["T1"],
                generate_plots=False,
                verbose=True,
            )

            results = generate_reports(config)

            self.assertEqual(len(results["readmes_generated"]), 1)
            self.assertEqual(len(results["tables_generated"]), 2)

            readme = reports_dir / "T1" / "README.md"
            self.assertTrue(readme.exists())
            content = readme.read_text()
            self.assertIn("llmap", content)
            self.assertIn("minimap2", content)

            tsv = reports_dir / "comparison.tsv"
            self.assertTrue(tsv.exists())


if __name__ == "__main__":
    unittest.main()
