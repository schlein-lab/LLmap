#!/usr/bin/env python3
"""SLURM submission orchestrator for LLmap benchmark campaign.

Generates and optionally submits SLURM sbatch jobs for the benchmark matrix.
Supports dry-run mode, job dependency chains, and idempotent re-runs.

Usage:
    python orchestrate.py                          # auto-detect mode
    python orchestrate.py --dry-run                # print jobs, don't submit
    python orchestrate.py --local                  # run locally (no SLURM)
    python orchestrate.py --tasks T1 T2            # subset of tasks
    python orchestrate.py --tools llmap minimap2   # subset of tools
    python orchestrate.py --deps                   # add metric job dependencies
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional


@dataclass
class BenchmarkCell:
    """Single cell in the benchmark matrix."""
    task: str
    tool: str
    preset: str
    ref_id: str
    reads_id: str
    replicate: int = 0


@dataclass
class SlurmJob:
    """SLURM job specification."""
    name: str
    script_path: Path
    output_dir: Path
    cpus: int = 16
    mem_gb: int = 64
    time_hours: int = 4
    partition: str = "standard"
    gpu: bool = False
    dependency: Optional[str] = None  # "afterok:jobid" format


@dataclass
class OrchestratorConfig:
    """Orchestrator configuration."""
    llmap_root: Path
    mode: str  # "auto", "local", "hummel"
    dry_run: bool = False
    tasks_filter: list = field(default_factory=list)
    tools_filter: list = field(default_factory=list)
    replicates: int = 3
    threads: int = 16
    seed_base: int = 42
    add_dependencies: bool = False


# Benchmark matrix definition
BENCHMARK_MATRIX = [
    # T1: Synthetic WGS
    ("T1", "llmap", "map-hifi", "synth_t1", "fastq"),
    ("T1", "minimap2", "map-hifi", "synth_t1", "fastq"),
    ("T1", "winnowmap2", "map-pb", "synth_t1", "fastq"),
    # T2: Synthetic paralog stress
    ("T2", "llmap", "map-hifi", "synth_t2", "fastq"),
    ("T2", "minimap2", "map-hifi", "synth_t2", "fastq"),
    ("T2", "winnowmap2", "map-pb", "synth_t2", "fastq"),
    # T3: Real HiFi WGS
    ("T3", "llmap", "map-hifi", "grch38", "hg002_hifi_chr14_chr20"),
    ("T3", "minimap2", "map-hifi", "grch38", "hg002_hifi_chr14_chr20"),
    ("T3", "winnowmap2", "map-pb", "grch38", "hg002_hifi_chr14_chr20"),
    # T4: Real Illumina WGS
    ("T4", "llmap", "sr", "grch38", "hg002_illumina_chr14_chr20"),
    ("T4", "minimap2", "sr", "grch38", "hg002_illumina_chr14_chr20"),
    ("T4", "bwa_mem2", "-", "grch38", "hg002_illumina_chr14_chr20"),
    ("T4", "bowtie2", "-", "grch38", "hg002_illumina_chr14_chr20"),
    # T5: Iso-seq transcriptome
    ("T5", "llmap", "splice", "gencode_v46", "hprc_isoseq_hg00290"),
    ("T5", "minimap2", "splice", "gencode_v46", "hprc_isoseq_hg00290"),
    ("T5", "star", "-", "gencode_v46", "hprc_isoseq_hg00290"),
    # T6: Targeted IGH-locus
    ("T6", "llmap", "map-hifi", "igh_locus", "hg002_hifi_chr14_chr20"),
    ("T6", "minimap2", "map-hifi", "igh_locus", "hg002_hifi_chr14_chr20"),
    ("T6", "winnowmap2", "map-pb", "igh_locus", "hg002_hifi_chr14_chr20"),
]


def detect_mode() -> str:
    """Auto-detect execution mode based on environment."""
    if shutil.which("sbatch"):
        return "hummel"
    return "local"


def is_cell_done(output_dir: Path) -> bool:
    """Check if a benchmark cell is already complete."""
    return (output_dir / "done.flag").exists()


def generate_sbatch_script(job: SlurmJob, cell: BenchmarkCell, config: OrchestratorConfig) -> str:
    """Generate SLURM sbatch script content."""
    seed = config.seed_base + cell.replicate
    runners_dir = config.llmap_root / "benchmarks" / "runners"
    runner_script = runners_dir / f"run_{cell.tool}.sh"

    gpu_directive = ""
    if job.gpu:
        gpu_directive = "#SBATCH --gres=gpu:1"

    dependency_directive = ""
    if job.dependency:
        dependency_directive = f"#SBATCH --dependency={job.dependency}"

    script = f"""#!/bin/bash
#SBATCH --job-name={job.name}
#SBATCH --output={job.output_dir}/slurm.out
#SBATCH --error={job.output_dir}/slurm.err
#SBATCH --cpus-per-task={job.cpus}
#SBATCH --mem={job.mem_gb}G
#SBATCH --time={job.time_hours}:00:00
#SBATCH --partition={job.partition}
{gpu_directive}
{dependency_directive}

set -euo pipefail

export TOOL="{cell.tool}"
export TASK="{cell.task}"
export REPLICATE="{cell.replicate}"
export PRESET="{cell.preset}"
export REF_ID="{cell.ref_id}"
export READS_ID="{cell.reads_id}"
export THREADS="{config.threads}"
export SEED="{seed}"
export OUTPUT_DIR="{job.output_dir}"

echo "Starting benchmark: $TASK/$TOOL/rep$REPLICATE"
echo "Time: $(date -Iseconds)"

bash "{runner_script}"

echo "Completed: $TASK/$TOOL/rep$REPLICATE"
echo "Time: $(date -Iseconds)"
"""
    return script.strip() + "\n"


def generate_metrics_script(task: str, config: OrchestratorConfig, depend_jobs: list) -> str:
    """Generate SLURM script for metrics aggregation."""
    reports_dir = config.llmap_root / "benchmarks" / "reports" / task
    metrics_dir = config.llmap_root / "benchmarks" / "metrics"

    deps = ",".join(depend_jobs) if depend_jobs else ""
    dep_directive = f"#SBATCH --dependency=afterok:{deps}" if deps else ""

    script = f"""#!/bin/bash
#SBATCH --job-name=metrics_{task}
#SBATCH --output={reports_dir}/metrics_slurm.out
#SBATCH --error={reports_dir}/metrics_slurm.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
{dep_directive}

set -euo pipefail

cd "{metrics_dir}"

echo "Computing metrics for task {task}"
echo "Time: $(date -Iseconds)"

python compute.py --task {task} --reports-dir "{reports_dir}"
python concordance.py --task {task} --reports-dir "{reports_dir}"

echo "Metrics complete for task {task}"
echo "Time: $(date -Iseconds)"
"""
    return script.strip() + "\n"


def expand_matrix(config: OrchestratorConfig) -> list:
    """Expand the benchmark matrix into individual cells."""
    cells = []
    for task, tool, preset, ref_id, reads_id in BENCHMARK_MATRIX:
        if config.tasks_filter and task not in config.tasks_filter:
            continue
        if config.tools_filter and tool not in config.tools_filter:
            continue
        for rep in range(config.replicates):
            cells.append(BenchmarkCell(
                task=task,
                tool=tool,
                preset=preset,
                ref_id=ref_id,
                reads_id=reads_id,
                replicate=rep,
            ))
    return cells


def create_job(cell: BenchmarkCell, config: OrchestratorConfig) -> SlurmJob:
    """Create a SLURM job specification for a benchmark cell."""
    reports_dir = config.llmap_root / "benchmarks" / "reports"
    output_dir = reports_dir / cell.task / cell.tool / f"rep{cell.replicate}"
    script_path = output_dir / "submit.sbatch"

    gpu = (cell.tool == "llmap" and cell.task in ("T3", "T4", "T5", "T6"))

    return SlurmJob(
        name=f"bench_{cell.task}_{cell.tool}_{cell.replicate}",
        script_path=script_path,
        output_dir=output_dir,
        gpu=gpu,
    )


def submit_job(job: SlurmJob, script_content: str, dry_run: bool) -> Optional[str]:
    """Submit a SLURM job and return the job ID."""
    job.output_dir.mkdir(parents=True, exist_ok=True)
    job.script_path.write_text(script_content)

    if dry_run:
        print(f"[DRY-RUN] Would submit: {job.script_path}")
        return None

    result = subprocess.run(
        ["sbatch", str(job.script_path)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"ERROR submitting {job.name}: {result.stderr}", file=sys.stderr)
        return None

    # Parse job ID from "Submitted batch job 12345"
    output = result.stdout.strip()
    if "Submitted batch job" in output:
        job_id = output.split()[-1]
        print(f"Submitted {job.name}: job {job_id}")
        return job_id

    print(f"WARNING: Could not parse job ID from: {output}")
    return None


def run_local(cell: BenchmarkCell, config: OrchestratorConfig, dry_run: bool) -> bool:
    """Run a benchmark cell locally (no SLURM)."""
    reports_dir = config.llmap_root / "benchmarks" / "reports"
    output_dir = reports_dir / cell.task / cell.tool / f"rep{cell.replicate}"

    if is_cell_done(output_dir):
        print(f"[SKIP] {cell.task}/{cell.tool}/rep{cell.replicate} already done")
        return True

    output_dir.mkdir(parents=True, exist_ok=True)

    seed = config.seed_base + cell.replicate
    runners_dir = config.llmap_root / "benchmarks" / "runners"
    runner_script = runners_dir / f"run_{cell.tool}.sh"

    env = os.environ.copy()
    env.update({
        "TOOL": cell.tool,
        "TASK": cell.task,
        "REPLICATE": str(cell.replicate),
        "PRESET": cell.preset,
        "REF_ID": cell.ref_id,
        "READS_ID": cell.reads_id,
        "THREADS": str(config.threads),
        "SEED": str(seed),
        "OUTPUT_DIR": str(output_dir),
    })

    cmd = ["bash", str(runner_script)]

    if dry_run:
        print(f"[DRY-RUN] Would run: {' '.join(cmd)}")
        print(f"          Output: {output_dir}")
        return True

    print(f"[LOCAL] {cell.task}/{cell.tool}/rep{cell.replicate}")
    result = subprocess.run(cmd, env=env)
    return result.returncode == 0


def orchestrate(config: OrchestratorConfig) -> dict:
    """Main orchestration entry point."""
    cells = expand_matrix(config)
    results = {
        "mode": config.mode,
        "dry_run": config.dry_run,
        "timestamp": datetime.now().isoformat(),
        "cells_total": len(cells),
        "cells_skipped": 0,
        "cells_submitted": 0,
        "cells_failed": 0,
        "job_ids": [],
        "task_job_ids": {},  # task -> list of job IDs for dependency chaining
    }

    if config.mode == "local":
        for cell in cells:
            output_dir = (config.llmap_root / "benchmarks" / "reports" /
                         cell.task / cell.tool / f"rep{cell.replicate}")
            if is_cell_done(output_dir):
                results["cells_skipped"] += 1
                continue

            success = run_local(cell, config, config.dry_run)
            if success or config.dry_run:
                results["cells_submitted"] += 1
            else:
                results["cells_failed"] += 1

    elif config.mode == "hummel":
        for cell in cells:
            job = create_job(cell, config)

            if is_cell_done(job.output_dir):
                print(f"[SKIP] {cell.task}/{cell.tool}/rep{cell.replicate} already done")
                results["cells_skipped"] += 1
                continue

            script = generate_sbatch_script(job, cell, config)
            job_id = submit_job(job, script, config.dry_run)

            if job_id or config.dry_run:
                results["cells_submitted"] += 1
                if job_id:
                    results["job_ids"].append(job_id)
                    if cell.task not in results["task_job_ids"]:
                        results["task_job_ids"][cell.task] = []
                    results["task_job_ids"][cell.task].append(job_id)
            else:
                results["cells_failed"] += 1

        if config.add_dependencies and results["task_job_ids"]:
            results["metrics_jobs"] = {}
            for task, job_ids in results["task_job_ids"].items():
                reports_dir = config.llmap_root / "benchmarks" / "reports" / task
                reports_dir.mkdir(parents=True, exist_ok=True)

                script = generate_metrics_script(task, config, job_ids)
                script_path = reports_dir / "metrics_submit.sbatch"
                script_path.write_text(script)

                if config.dry_run:
                    print(f"[DRY-RUN] Would submit metrics job for {task} "
                          f"(depends on {len(job_ids)} alignment jobs)")
                else:
                    res = subprocess.run(
                        ["sbatch", str(script_path)],
                        capture_output=True,
                        text=True,
                    )
                    if res.returncode == 0 and "Submitted batch job" in res.stdout:
                        metrics_id = res.stdout.strip().split()[-1]
                        results["metrics_jobs"][task] = metrics_id
                        print(f"Submitted metrics job for {task}: {metrics_id}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="SLURM orchestrator for LLmap benchmark campaign"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Print jobs without submitting")
    parser.add_argument("--local", action="store_true",
                        help="Run locally (no SLURM)")
    parser.add_argument("--hummel", action="store_true",
                        help="Submit to SLURM")
    parser.add_argument("--tasks", nargs="+", default=[],
                        help="Filter to specific tasks (e.g., T1 T2)")
    parser.add_argument("--tools", nargs="+", default=[],
                        help="Filter to specific tools (e.g., llmap minimap2)")
    parser.add_argument("--replicates", type=int, default=3,
                        help="Number of replicates (default: 3)")
    parser.add_argument("--threads", type=int, default=16,
                        help="Threads per job (default: 16)")
    parser.add_argument("--deps", action="store_true",
                        help="Add metrics job dependencies")
    parser.add_argument("--json", action="store_true",
                        help="Output results as JSON")
    args = parser.parse_args()

    llmap_root = Path(__file__).parent.parent.resolve()

    if args.local:
        mode = "local"
    elif args.hummel:
        mode = "hummel"
    else:
        mode = detect_mode()

    config = OrchestratorConfig(
        llmap_root=llmap_root,
        mode=mode,
        dry_run=args.dry_run,
        tasks_filter=args.tasks,
        tools_filter=args.tools,
        replicates=args.replicates,
        threads=args.threads,
        add_dependencies=args.deps,
    )

    results = orchestrate(config)

    if args.json:
        print(json.dumps(results, indent=2, default=str))
    else:
        print(f"\nOrchestration summary:")
        print(f"  Mode: {results['mode']}")
        print(f"  Dry run: {results['dry_run']}")
        print(f"  Total cells: {results['cells_total']}")
        print(f"  Skipped (done): {results['cells_skipped']}")
        print(f"  Submitted: {results['cells_submitted']}")
        print(f"  Failed: {results['cells_failed']}")
        if results.get("job_ids"):
            print(f"  Job IDs: {', '.join(results['job_ids'][:10])}"
                  f"{'...' if len(results['job_ids']) > 10 else ''}")


if __name__ == "__main__":
    main()
