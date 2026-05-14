# LLmap Comparative Benchmark Campaign — Phase 11

This is the specification for the Phase 11 benchmark campaign. Goal: a systematic, reproducible, publication-ready comparison of LLmap against the five dominant production mappers across a representative set of mapping tasks. All commands, dataset versions, tool versions, and parameter choices are pinned here so the campaign can be re-run on any cluster.

---

## Tools under test

| Tool | Version | Role | URL |
|---|---|---|---|
| **LLmap** | 1.0.0 | this work | https://github.com/schlein-lab/LLmap |
| **minimap2** | 2.28 | long + short, gold standard | https://github.com/lh3/minimap2 |
| **BWA-MEM2** | 2.2.1 | short-read gold standard | https://github.com/bwa-mem2/bwa-mem2 |
| **Winnowmap2** | 2.03 | long-read repeat-aware (LLmap's main competitor for SD regions) | https://github.com/marbl/Winnowmap |
| **STAR** | 2.7.11b | RNA-seq spliced alignment | https://github.com/alexdobin/STAR |
| **Bowtie2** | 2.5.4 | short-read, BWT-based alternative algorithm | https://github.com/BenLangmead/bowtie2 |

All tool versions are pinned via [benchmarks/datasets/tools.yaml](datasets/tools.yaml) (conda + module fallback). Each runner reports `<tool> --version` into its result manifest.

---

## Tasks

Six tasks span the operating space of a production mapper. Not every tool runs every task — a tool that does not natively support a data type is simply absent from that row in the results matrix.

| ID | Task | Input | Reference | Tools |
|---|---|---|---|---|
| **T1** | Synthetic-truth WGS | 1M simulated HiFi reads, 30× equivalent, exact origin recorded | GRCh38 chr20 + chr14 (IGH locus) | LLmap, minimap2, Winnowmap2 |
| **T2** | Synthetic paralog stress | 500k simulated reads from IGHG1/2/3/4 + IGHGP, NPHP1/NPHP1B, MHC class I duplicates — ground truth per-paralog assignment | IGH locus + MHC | LLmap, minimap2, Winnowmap2 |
| **T3** | Real HiFi WGS | HG002 HiFi reads, chr14 + chr20 subset (~5M reads) | GRCh38 | LLmap, minimap2, Winnowmap2 |
| **T4** | Real Illumina WGS | HG002 Illumina 2×150 reads, chr14 + chr20 subset (~10M pairs) | GRCh38 | LLmap, minimap2, BWA-MEM2, Bowtie2 |
| **T5** | Iso-seq transcriptome | HPRC HG00290 lymph FLNC reads (paralog-rich for IGH constant exons) | GRCh38 transcripts (GENCODE v46) | LLmap, minimap2, STAR |
| **T6** | Targeted IGH-locus | Real reads pre-selected to IGH region (chr14:105.5–107.3 Mb) — paralog disambiguation challenge | IGH locus extract + IGHGP / IGHGY / IGHEP1 / IGHEP2 pseudogenes | LLmap, minimap2, Winnowmap2 |

For each task, all tools are invoked with their published "recommended for this data type" presets:

- minimap2: `-x map-hifi`, `-x map-ont`, `-x sr`, `-x splice` as appropriate
- LLmap: `-x map-hifi` etc. (same surface)
- Winnowmap2: `-W repetitive_k15.txt -ax map-pb` per its README
- BWA-MEM2: `mem -t N`
- STAR: `--outSAMtype BAM SortedByCoordinate --outSJtype Standard` etc.
- Bowtie2: `--very-sensitive`

The exact command line per task per tool is in [benchmarks/runners/](runners/) and is logged into every output manifest.

---

## Metrics

Five metric families. Each task records all metrics that apply. The metric-computation code lives in [benchmarks/metrics/](metrics/) and produces a single `result.json` per (tool × task × replicate) run.

### 1. Performance

| Metric | Unit | How measured |
|---|---|---|
| `wallclock_seconds` | s | `/usr/bin/time -v` |
| `cpu_seconds` | s | `time -v %S + %U` |
| `peak_rss_bytes` | bytes | `time -v %M × 1024` |
| `max_threads_used` | int | from `pidstat` sample during run |
| `gpu_seconds` | s | `nvidia-smi pmon` sample (LLmap only when CUDA on) |
| `gpu_peak_mem_bytes` | bytes | `nvidia-smi --query-gpu=memory.used` |
| `output_file_bytes` | bytes | `stat -c %s` |

All runs pinned to N=16 threads on identical Hummel partitions (default `--cpus-per-task=16 --mem=64G`) to keep comparisons fair. Three replicates per (tool, task).

### 2. Mapping coverage

| Metric | Definition |
|---|---|
| `primary_mapped` | count of records with `FLAG & 0x4 == 0` and not secondary/supplementary |
| `secondary_count` | count of `FLAG & 0x100` records |
| `supplementary_count` | count of `FLAG & 0x800` records |
| `unmapped` | count of `FLAG & 0x4` records |
| `mapping_rate` | `primary_mapped / total_input_reads` |
| `drop_rate` | `1 - (primary_mapped + secondary + supplementary) / total_input_reads` — captures information loss, not just "unmapped" |

LLmap is special-cased: a read that LLmap reports as `uninformative` is counted as `primary_mapped=1` with `mapq=0` in BAM space *and* as a probability-distribution row in the Parquet output. Both views are reported separately so the comparison is honest.

### 3. MAPQ distribution

Per tool: histogram of MAPQ values 0…60. Reported as ECDF and as the distribution mean/median/p90. For LLmap the BAM MAPQ comes from the WaveCollapse posterior:

```
MAPQ = round(-10 * log10(1 - max_b P(b|r)))    capped at 60
```

The Parquet output retains the full posterior; comparing only the collapsed MAPQ is the lossy view.

### 4. Pairwise concordance

For every (toolA, toolB) pair on the same task, compute the agreement matrix on primary-mapped reads:

| Class | Definition |
|---|---|
| `same_position` | both tools map to same chromosome, position within ±10 bp |
| `same_chrom_diff_pos` | same chrom, distance > 10 bp |
| `diff_chrom` | different chrom |
| `oneA_oneB` | mapped by exactly one of the two |
| `both_unmapped` | both report unmapped |

The asymmetric matrix is the primary "where do they disagree?" view. For paralog-rich loci (T2, T6) the disagreement count is expected to be very high — that is the point.

### 5. Ground-truth accuracy

For synthetic tasks (T1, T2) the simulator stamps the true origin into each read name. The accuracy metrics:

| Metric | Definition |
|---|---|
| `recall` | fraction of input reads with primary alignment within ±10 bp of truth |
| `precision` | fraction of mapped reads whose primary alignment is within ±10 bp of truth |
| `f1` | 2·P·R/(P+R) |
| `paralog_correct` | fraction of paralog reads (T2) assigned to the correct paralog member |
| `paralog_calibration` | (LLmap only) Brier score between LLmap posterior and one-hot truth |

For T5 (transcriptome), splice junction F1:

| Metric | Definition |
|---|---|
| `junction_recall` | fraction of GENCODE introns supported by ≥1 read with correctly placed split |
| `junction_precision` | fraction of called junctions that match GENCODE |

For T3 / T4 / T6 (real data, no per-read truth) we substitute proxy truths:

- T3/T4: consensus of all tools on uniquely-mappable regions (90% of GRCh38 outside repeats) is the proxy truth. Per-tool deviation from this consensus is reported. NOT used for accuracy claims, only for divergence analysis.
- T6: per-read paralog truth derived from PSV genotyping (see `pseudocaller` integration). 

---

## Datasets

Pinned data sources in [benchmarks/datasets/](datasets/):

| ID | Source | Size | Location on Hummel |
|---|---|---|---|
| `ref_grch38` | NCBI GRCh38.p14 | 3.1 GB | `/beegfs/u/bbg6775/humangenetik/shared/references/GRCh38/` |
| `ref_chm13` | T2T CHM13 v2.0 | 3.0 GB | `/beegfs/u/bbg6775/humangenetik/shared/references/CHM13/` |
| `ref_igh_locus` | GRCh38 chr14:105.5–107.3 Mb extract | 1.8 MB | generated by `make_igh_locus.sh` |
| `ref_gencode_v46` | GENCODE v46 transcripts | 250 MB | `/beegfs/u/bbg6775/humangenetik/shared/references/gencode/v46/` |
| `hg002_hifi` | HG002 HiFi (Q20+), chr14+chr20 subset | 12 GB | `/beegfs/u/bbg6775/humangenetik/shared/HG002/hifi/` |
| `hg002_illumina` | HG002 NovaSeq 2×150 30×, chr14+chr20 subset | 8 GB | `/beegfs/u/bbg6775/humangenetik/shared/HG002/illumina/` |
| `hprc_isoseq_hg00290` | HPRC lymph FLNC | 1.7 GB | `/beegfs/u/bbg6775/humangenetik/shared/transcriptome_longread/hprc_isoseq/HG00290/` |
| `synth_truth_wgs` | generated by `llmap generate-synth --task t1` | ~10 GB | generated |
| `synth_paralog_stress` | generated by `llmap generate-synth --task t2 --paralog` | ~5 GB | generated |

Generation commands are in `benchmarks/datasets/build.sh`.

---

## Output layout

Each run produces a directory:

```
benchmarks/reports/<run_id>/
  manifest.json                          # tool, version, command, dataset SHA256, host, seed
  resources.json                         # wallclock, RAM, GPU, etc.
  alignments.bam                         # standardized to sorted+indexed BAM
  mapping_summary.json                   # primary/secondary/unmapped counts
  mapq_histogram.json                    # 0..60 bin counts
  ground_truth.json                      # recall/precision/F1/paralog (if applicable)
  per_read_truth.parquet                 # per-read: read_id, truth_pos, called_pos, distance, correct
```

Cross-tool aggregation produces:

```
benchmarks/reports/<task>/
  concordance_matrix.json                # all pairwise (toolA, toolB) agreement
  performance_table.tsv                  # one row per tool, all metrics
  divergence_examples/                   # 100 reads where tools disagree, for inspection
  README.md                              # human-readable summary
```

Final publication-ready outputs land in:

```
docs/BENCHMARKS.md                       # text + tables, ready to commit
docs/figures/                            # SVG/PNG plots
docs/benchmarks_supplement.pdf           # full appendix
```

---

## Execution model

Local-host work: dataset generation scripts, runner skeletons, metric computation code, plot generation. All testable on CPU.

Hummel-2 work: every tool run that requires real data. Each (tool, task, replicate) is one SLURM job submitted from `benchmarks/runners/submit_all.sh`. The job array uses a single GPU partition for LLmap CUDA runs and the CPU partition for all others.

Submission idempotency: each job writes a `done.flag` on success. `submit_all.sh` skips runs whose flag already exists. Failed runs are deleted (the flag + output dir) and re-submitted.

Wallclock budget for full matrix: ~24 h on 16 nodes (rough estimate, dominated by full-WGS tasks T3/T4). Phase 11 substeps fit incremental work into the autonomous-driver loop.

---

## Phase 11 substep plan

| Substep | Deliverable | Effort estimate |
|---|---|---|
| 11.1 | This SPEC.md + dataset/tool registries + runner template | done |
| 11.2 | Synthetic-truth data generator (extend Phase 0.4 with truth stamping for T1, T2) | 2 iter |
| 11.3 | Tool installation manifest (conda env + module load) + version-verification script | 1 iter |
| 11.4 | Per-tool runner scripts (minimap2.sh, bwa_mem2.sh, winnowmap2.sh, star.sh, bowtie2.sh, llmap.sh) | 3 iter |
| 11.5 | Metrics collector: SAM/BAM → mapping_summary, MAPQ histogram, ground_truth | 2 iter |
| 11.6 | Concordance analyzer: pairwise agreement matrices, divergence example extraction | 2 iter |
| 11.7 | Report generator: per-task README.md + cross-tool tables + plots | 2 iter |
| 11.8 | SLURM submission orchestrator (submit_all.sh + per-task arrays) | 1 iter |
| 11.9 | Run synthetic tasks T1, T2 (local CPU, no Hummel needed) | 2 iter |
| 11.10 | Submit real-data tasks T3–T6 to Hummel (manual step, user submits) | 1 iter + manual |
| 11.11 | Aggregate results → docs/BENCHMARKS.md + supplement | 2 iter |
| 11.12 | Identify performance / quality regressions → LLmap improvement issues | 1 iter |

Total: ~18 iterations (~4.5 hours of autonomous-loop time + manual Hummel submission).

---

## What "publication-ready" means here

The output of Phase 11 must let an external reader reproduce every number:

1. All tool versions pinned, with SHA256 of the binary recorded.
2. All datasets pinned, with SHA256 of the input file recorded.
3. All command lines recorded verbatim in each `manifest.json`.
4. All random seeds (simulators, downsampling) recorded.
5. Three independent replicates per (tool, task) cell. Variance reported as min/median/max wallclock and bootstrap CIs on accuracy.
6. Every plot has a script in `benchmarks/metrics/plot_*.py` that regenerates it from `result.json` files.
7. Tied / overlapping results are not hidden — the publication-ready table is honest. Where LLmap loses, that loss is the input for Phase 12 (improvement targets).

If a reader downloads the repo, installs the tools per the manifest, downloads the datasets per the build script, and runs `benchmarks/run_all.sh`, they get the same numbers (up to variance from non-deterministic threading and GPU scheduling).
