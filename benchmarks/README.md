# LLmap Benchmarks

Phase 11 comparative benchmark campaign. See [SPEC.md](SPEC.md) for the full specification (tools, tasks, metrics, datasets, output layout).

## Quick start

```bash
# 1. Verify tool versions
./runners/check_versions.sh

# 2. Generate synthetic datasets with ground truth (T1, T2)
llmap generate-synth --task t1 --output datasets/cache/synth_t1/
llmap generate-synth --task t2 --paralog igh,mhc --output datasets/cache/synth_t2/

# 3. Run a single cell (dry-run print first)
./runners/submit_all.sh --task T1 --dry-run

# 4. Run synthetic tasks locally (CPU only)
./runners/submit_all.sh --local --task T1 --task T2

# 5. Submit real-data tasks to Hummel
./runners/submit_all.sh --hummel --task T3 --task T4 --task T5 --task T6

# 6. Compute metrics (after BAMs exist)
for run in reports/T*/*/rep*; do
    python3 metrics/compute.py --bam "$run/alignments.bam" \
        --total $(grep total_input "$run/manifest.json" | ...) \
        --truth datasets/cache/synth_t1/truth.tsv \
        --out-dir "$run/"
done

# 7. Pairwise concordance
python3 metrics/concordance.py \
    --bam-a reports/T1/llmap/rep0/alignments.bam --name-a llmap \
    --bam-b reports/T1/minimap2/rep0/alignments.bam --name-b minimap2 \
    --out-dir reports/T1/concordance_llmap_vs_minimap2/
```

## Status

Scaffolding committed. Implementation tracked in [STATE.md](../STATE.md) Phase 11.

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
