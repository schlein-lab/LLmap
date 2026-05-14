#!/usr/bin/env python3
"""Tests for the SLURM orchestrator.

Run: python -m pytest test_orchestrate.py -v
"""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from orchestrate import (
    BENCHMARK_MATRIX,
    BenchmarkCell,
    OrchestratorConfig,
    SlurmJob,
    create_job,
    detect_mode,
    expand_matrix,
    generate_metrics_script,
    generate_sbatch_script,
    is_cell_done,
    orchestrate,
)


class TestBenchmarkMatrix:
    """Tests for the benchmark matrix definition."""

    def test_matrix_has_expected_tasks(self):
        """Matrix includes all 6 benchmark tasks."""
        tasks = {row[0] for row in BENCHMARK_MATRIX}
        assert tasks == {"T1", "T2", "T3", "T4", "T5", "T6"}

    def test_matrix_has_expected_tools(self):
        """Matrix includes all expected tools."""
        tools = {row[1] for row in BENCHMARK_MATRIX}
        assert "llmap" in tools
        assert "minimap2" in tools
        assert "winnowmap2" in tools
        assert "bwa_mem2" in tools
        assert "star" in tools
        assert "bowtie2" in tools

    def test_t1_synthetic_has_llmap_minimap2_winnowmap(self):
        """T1 includes the expected tool set."""
        t1_tools = {row[1] for row in BENCHMARK_MATRIX if row[0] == "T1"}
        assert t1_tools == {"llmap", "minimap2", "winnowmap2"}

    def test_t4_illumina_has_short_read_tools(self):
        """T4 includes short-read alignment tools."""
        t4_tools = {row[1] for row in BENCHMARK_MATRIX if row[0] == "T4"}
        assert "bwa_mem2" in t4_tools
        assert "bowtie2" in t4_tools

    def test_t5_isoseq_has_star(self):
        """T5 transcriptome task includes STAR."""
        t5_tools = {row[1] for row in BENCHMARK_MATRIX if row[0] == "T5"}
        assert "star" in t5_tools


class TestExpandMatrix:
    """Tests for matrix expansion."""

    def test_expand_all_cells(self):
        """Expand produces cells for all matrix rows and replicates."""
        config = OrchestratorConfig(
            llmap_root=Path("/fake"),
            mode="local",
            replicates=3,
        )
        cells = expand_matrix(config)
        assert len(cells) == len(BENCHMARK_MATRIX) * 3

    def test_expand_single_replicate(self):
        """Expand with 1 replicate produces one cell per row."""
        config = OrchestratorConfig(
            llmap_root=Path("/fake"),
            mode="local",
            replicates=1,
        )
        cells = expand_matrix(config)
        assert len(cells) == len(BENCHMARK_MATRIX)

    def test_expand_filter_by_task(self):
        """Expand respects task filter."""
        config = OrchestratorConfig(
            llmap_root=Path("/fake"),
            mode="local",
            tasks_filter=["T1"],
            replicates=1,
        )
        cells = expand_matrix(config)
        assert all(c.task == "T1" for c in cells)
        assert len(cells) == 3  # llmap, minimap2, winnowmap2

    def test_expand_filter_by_tool(self):
        """Expand respects tool filter."""
        config = OrchestratorConfig(
            llmap_root=Path("/fake"),
            mode="local",
            tools_filter=["llmap"],
            replicates=1,
        )
        cells = expand_matrix(config)
        assert all(c.tool == "llmap" for c in cells)
        assert len(cells) == 6  # one per task

    def test_expand_combined_filters(self):
        """Expand respects both task and tool filters."""
        config = OrchestratorConfig(
            llmap_root=Path("/fake"),
            mode="local",
            tasks_filter=["T1", "T2"],
            tools_filter=["llmap", "minimap2"],
            replicates=2,
        )
        cells = expand_matrix(config)
        # T1: llmap, minimap2. T2: llmap, minimap2. Each with 2 reps.
        assert len(cells) == 8


class TestBenchmarkCell:
    """Tests for BenchmarkCell dataclass."""

    def test_cell_creation(self):
        """Cell stores all required fields."""
        cell = BenchmarkCell(
            task="T1",
            tool="llmap",
            preset="map-hifi",
            ref_id="synth_t1",
            reads_id="fastq",
            replicate=0,
        )
        assert cell.task == "T1"
        assert cell.tool == "llmap"
        assert cell.replicate == 0


class TestCreateJob:
    """Tests for SLURM job creation."""

    def test_create_job_basic(self):
        """create_job produces a valid SlurmJob."""
        cell = BenchmarkCell(
            task="T1",
            tool="llmap",
            preset="map-hifi",
            ref_id="synth_t1",
            reads_id="fastq",
            replicate=0,
        )
        config = OrchestratorConfig(
            llmap_root=Path("/fake/llmap"),
            mode="hummel",
        )
        job = create_job(cell, config)
        assert job.name == "bench_T1_llmap_0"
        assert job.output_dir == Path("/fake/llmap/benchmarks/reports/T1/llmap/rep0")
        assert job.script_path == Path("/fake/llmap/benchmarks/reports/T1/llmap/rep0/submit.sbatch")

    def test_create_job_gpu_for_real_data_llmap(self):
        """LLmap on real data tasks gets GPU."""
        cell = BenchmarkCell(
            task="T3",
            tool="llmap",
            preset="map-hifi",
            ref_id="grch38",
            reads_id="hg002_hifi",
            replicate=0,
        )
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        job = create_job(cell, config)
        assert job.gpu is True

    def test_create_job_no_gpu_for_synthetic(self):
        """Synthetic tasks don't need GPU."""
        cell = BenchmarkCell(
            task="T1",
            tool="llmap",
            preset="map-hifi",
            ref_id="synth_t1",
            reads_id="fastq",
            replicate=0,
        )
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        job = create_job(cell, config)
        assert job.gpu is False

    def test_create_job_no_gpu_for_other_tools(self):
        """Non-LLmap tools don't use GPU."""
        cell = BenchmarkCell(
            task="T3",
            tool="minimap2",
            preset="map-hifi",
            ref_id="grch38",
            reads_id="hg002_hifi",
            replicate=0,
        )
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        job = create_job(cell, config)
        assert job.gpu is False


class TestGenerateSbatchScript:
    """Tests for sbatch script generation."""

    def test_generate_script_contains_directives(self):
        """Generated script has required SLURM directives."""
        cell = BenchmarkCell(
            task="T1",
            tool="llmap",
            preset="map-hifi",
            ref_id="synth_t1",
            reads_id="fastq",
            replicate=0,
        )
        config = OrchestratorConfig(
            llmap_root=Path("/fake/llmap"),
            mode="hummel",
            threads=16,
            seed_base=42,
        )
        job = SlurmJob(
            name="bench_T1_llmap_0",
            script_path=Path("/fake/submit.sbatch"),
            output_dir=Path("/fake/output"),
        )
        script = generate_sbatch_script(job, cell, config)

        assert "#SBATCH --job-name=bench_T1_llmap_0" in script
        assert "#SBATCH --cpus-per-task=16" in script
        assert "#SBATCH --mem=64G" in script
        assert 'export TOOL="llmap"' in script
        assert 'export TASK="T1"' in script
        assert 'export SEED="42"' in script

    def test_generate_script_gpu_directive(self):
        """Script includes GPU directive when requested."""
        cell = BenchmarkCell(
            task="T3", tool="llmap", preset="map-hifi",
            ref_id="grch38", reads_id="hg002", replicate=0,
        )
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        job = SlurmJob(
            name="test", script_path=Path("/fake"), output_dir=Path("/fake"),
            gpu=True,
        )
        script = generate_sbatch_script(job, cell, config)
        assert "#SBATCH --gres=gpu:1" in script

    def test_generate_script_dependency_directive(self):
        """Script includes dependency when specified."""
        cell = BenchmarkCell(
            task="T1", tool="llmap", preset="map-hifi",
            ref_id="synth_t1", reads_id="fastq", replicate=0,
        )
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        job = SlurmJob(
            name="test", script_path=Path("/fake"), output_dir=Path("/fake"),
            dependency="afterok:12345",
        )
        script = generate_sbatch_script(job, cell, config)
        assert "#SBATCH --dependency=afterok:12345" in script


class TestGenerateMetricsScript:
    """Tests for metrics job script generation."""

    def test_generate_metrics_script_basic(self):
        """Metrics script has expected structure."""
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        script = generate_metrics_script("T1", config, [])

        assert "#SBATCH --job-name=metrics_T1" in script
        assert "python compute.py --task T1" in script
        assert "python concordance.py --task T1" in script

    def test_generate_metrics_script_with_deps(self):
        """Metrics script includes dependencies."""
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="hummel")
        script = generate_metrics_script("T1", config, ["123", "456", "789"])

        assert "#SBATCH --dependency=afterok:123,456,789" in script


class TestIsCellDone:
    """Tests for cell completion detection."""

    def test_is_done_when_flag_exists(self):
        """Cell is done when done.flag exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "done.flag").touch()
            assert is_cell_done(output_dir) is True

    def test_not_done_when_flag_missing(self):
        """Cell is not done when done.flag is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            assert is_cell_done(output_dir) is False


class TestDetectMode:
    """Tests for mode auto-detection."""

    @patch("shutil.which")
    def test_detect_hummel_when_sbatch_available(self, mock_which):
        """Detects hummel mode when sbatch is available."""
        mock_which.return_value = "/usr/bin/sbatch"
        assert detect_mode() == "hummel"

    @patch("shutil.which")
    def test_detect_local_when_sbatch_missing(self, mock_which):
        """Detects local mode when sbatch is unavailable."""
        mock_which.return_value = None
        assert detect_mode() == "local"


class TestOrchestrateDryRun:
    """Tests for orchestration in dry-run mode."""

    def test_dry_run_no_jobs_submitted(self):
        """Dry run does not submit any jobs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = OrchestratorConfig(
                llmap_root=Path(tmpdir),
                mode="hummel",
                dry_run=True,
                tasks_filter=["T1"],
                tools_filter=["llmap"],
                replicates=1,
            )
            results = orchestrate(config)

            assert results["dry_run"] is True
            assert results["cells_submitted"] == 1
            assert results["cells_failed"] == 0
            assert results["job_ids"] == []

    def test_dry_run_local_mode(self):
        """Dry run in local mode reports what would run."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = OrchestratorConfig(
                llmap_root=Path(tmpdir),
                mode="local",
                dry_run=True,
                tasks_filter=["T1"],
                tools_filter=["llmap"],
                replicates=1,
            )
            results = orchestrate(config)

            assert results["mode"] == "local"
            assert results["dry_run"] is True
            assert results["cells_submitted"] == 1


class TestOrchestrateSkipDone:
    """Tests for skipping completed cells."""

    def test_skip_completed_cells(self):
        """Orchestrator skips cells with done.flag."""
        with tempfile.TemporaryDirectory() as tmpdir:
            llmap_root = Path(tmpdir)
            reports_dir = llmap_root / "benchmarks" / "reports" / "T1" / "llmap" / "rep0"
            reports_dir.mkdir(parents=True)
            (reports_dir / "done.flag").touch()

            config = OrchestratorConfig(
                llmap_root=llmap_root,
                mode="hummel",
                dry_run=True,
                tasks_filter=["T1"],
                tools_filter=["llmap"],
                replicates=1,
            )
            results = orchestrate(config)

            assert results["cells_skipped"] == 1
            assert results["cells_submitted"] == 0


class TestDependencyChaining:
    """Tests for job dependency chain formation."""

    def test_task_job_ids_grouped(self):
        """Jobs are grouped by task for dependency chaining."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = OrchestratorConfig(
                llmap_root=Path(tmpdir),
                mode="hummel",
                dry_run=True,
                tasks_filter=["T1"],
                replicates=2,
                add_dependencies=True,
            )
            results = orchestrate(config)

            # In dry-run, job_ids is empty but cells are counted
            assert results["cells_submitted"] == 6  # 3 tools * 2 replicates


class TestOrchestratorConfig:
    """Tests for OrchestratorConfig."""

    def test_config_defaults(self):
        """Config has sensible defaults."""
        config = OrchestratorConfig(llmap_root=Path("/fake"), mode="local")
        assert config.dry_run is False
        assert config.replicates == 3
        assert config.threads == 16
        assert config.seed_base == 42
        assert config.add_dependencies is False

    def test_config_custom_values(self):
        """Config accepts custom values."""
        config = OrchestratorConfig(
            llmap_root=Path("/custom"),
            mode="hummel",
            dry_run=True,
            replicates=5,
            threads=32,
            seed_base=100,
        )
        assert config.llmap_root == Path("/custom")
        assert config.mode == "hummel"
        assert config.dry_run is True
        assert config.replicates == 5
        assert config.threads == 32
        assert config.seed_base == 100


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
