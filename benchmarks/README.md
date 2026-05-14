# LLmap Benchmarks

Phase 11 comparative benchmark campaign. See [SPEC.md](SPEC.md) for the full specification (tools, tasks, metrics, datasets, output layout).

## Quick Start

### Local synthetic benchmarks (T1, T2)

```bash
# Run synthetic benchmarks locally (CPU only)
./run_local_synthetic.sh

# Or step by step:
./runners/check_versions.sh                          # Verify tools
llmap generate-synth --task t1 --output datasets/cache/synth_t1/
llmap generate-synth --task t2 --paralog igh,mhc --output datasets/cache/synth_t2/
./runners/submit_all.sh --local --task T1 --task T2
./aggregate_results.sh
```

### Real-data benchmarks on Hummel (T3-T6)

**See [HUMMEL_SUBMISSION.md](HUMMEL_SUBMISSION.md) for detailed instructions.**

```bash
# On Hummel
cd ~/llmap/benchmarks

# Generate and submit all T3-T6 jobs
./hummel_submit_t3_t6.sh --submit

# Monitor status
./hummel_submit_t3_t6.sh --status
squeue -u $USER
```

### After all jobs complete

```bash
# Aggregate and generate reports
./aggregate_results.sh
python report.py --tasks T1 T2 T3 T4 T5 T6
```

## Status

- **T1, T2** (synthetic): Complete on local host
- **T3-T6** (real data): SLURM scripts ready, awaiting manual submission

Implementation tracked in [STATE.md](../STATE.md) Phase 11.

## Layout

```
benchmarks/
  SPEC.md            full specification
  README.md          this file
  datasets/
    tools.yaml       pinned tool versions
    datasets.yaml    pinned dataset locations
    build.sh         (TODO 11.2) dataset preparation
    cache/           generated/extracted artifacts (gitignored)
  runners/
    _template.sh     shared runner contract
    run_<tool>.sh    per-tool wrappers
    submit_all.sh    matrix orchestrator (local + sbatch)
    check_versions.sh
  metrics/
    compute.py       per-run mapping summary + ground-truth eval
    concordance.py   pairwise tool agreement
    plot_*.py        (TODO 11.7) plot generators
  reports/           run outputs (gitignored except aggregated summaries)
```
