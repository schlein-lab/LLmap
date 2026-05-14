# Hummel Submission Guide for T3-T6 Benchmarks

This guide explains how to submit the real-data benchmark tasks (T3-T6) to the Hummel-2 HPC cluster.

## Prerequisites

1. **SSH access to Hummel-2**: Ensure you can log in to `hummel-login`
2. **Clone or sync the repository** on Hummel:
   ```bash
   ssh hummel-login
   cd ~/llmap
   git pull origin main
   ```
3. **Tools installed**: Verify all tools are available:
   ```bash
   cd ~/llmap/benchmarks
   ./runners/check_versions.sh
   ```
4. **Datasets available**: Confirm paths in `datasets/datasets.yaml` exist on beegfs

## Quick Start

### Option 1: Using the submission script (recommended)

```bash
# On Hummel
cd ~/llmap/benchmarks

# Generate sbatch scripts (dry-run, shows what will be created)
./hummel_submit_t3_t6.sh --generate-only

# Submit all T3-T6 jobs to SLURM
./hummel_submit_t3_t6.sh --submit

# Monitor status
./hummel_submit_t3_t6.sh --status
```

### Option 2: Using the Python orchestrator

```bash
cd ~/llmap/benchmarks

# Dry-run (shows what would be submitted)
python orchestrate.py --hummel --tasks T3 T4 T5 T6 --dry-run

# Submit
python orchestrate.py --hummel --tasks T3 T4 T5 T6

# With dependency chaining for metrics
python orchestrate.py --hummel --tasks T3 T4 T5 T6 --deps
```

### Option 3: Manual submission

```bash
cd ~/llmap/benchmarks

# Generate scripts first
./hummel_submit_t3_t6.sh --generate-only

# Submit individually
for script in slurm_jobs/*.sbatch; do
  sbatch "$script"
done
```

## Task Overview

| Task | Description | Tools | Data | Time Est. |
|------|-------------|-------|------|-----------|
| T3 | Real HiFi WGS | llmap, minimap2, winnowmap2 | HG002 HiFi chr14+chr20 | ~2h/tool |
| T4 | Real Illumina WGS | llmap, minimap2, bwa_mem2, bowtie2 | HG002 Illumina chr14+chr20 | ~1h/tool |
| T5 | Iso-seq transcriptome | llmap, minimap2, star | HPRC HG00290 FLNC | ~1h/tool |
| T6 | Targeted IGH-locus | llmap, minimap2, winnowmap2 | HG002 HiFi (IGH region) | ~30m/tool |

Total: ~30 jobs (10 tool×task cells × 3 replicates)

## Monitoring

```bash
# Check queue status
squeue -u $USER

# Watch specific job
squeue -j <job_id>

# View job output (live)
tail -f benchmarks/reports/T3/llmap/rep0/slurm.out

# Check completion status
./hummel_submit_t3_t6.sh --status

# List completed runs
find benchmarks/reports -name "done.flag" | wc -l
```

## Resource Configuration

Default per-job resources (configurable in `hummel_submit_t3_t6.sh`):

- **CPUs**: 16 cores
- **Memory**: 64 GB
- **Time**: 8 hours
- **GPU**: 1 GPU for LLmap on T3-T6 (gpu partition)
- **Partition**: `standard` for CPU-only tools, `gpu` for LLmap

## Output Layout

Each run produces:

```
benchmarks/reports/<task>/<tool>/rep<N>/
├── alignments.bam      # Sorted+indexed BAM
├── alignments.bam.bai  # BAM index
├── manifest.json       # Tool version, command, dataset SHAs
├── resources.json      # Wallclock, CPU, peak RSS, output size
├── slurm.out           # SLURM stdout
├── slurm.err           # SLURM stderr
└── done.flag           # Completion marker (idempotency)
```

## After Completion

Once all jobs finish:

```bash
# Aggregate metrics
./aggregate_results.sh

# Generate reports
python report.py --task T3 T4 T5 T6

# Copy results back locally
scp -r hummel-login:~/llmap/benchmarks/reports/ benchmarks/
```

## Troubleshooting

### Job failed
```bash
# Check error log
cat benchmarks/reports/<task>/<tool>/rep<N>/slurm.err

# Remove done.flag to retry
rm benchmarks/reports/<task>/<tool>/rep<N>/done.flag

# Resubmit
sbatch slurm_jobs/<task>_<tool>_rep<N>.sbatch
```

### Missing datasets
```bash
# Verify paths exist
ls -l /beegfs/u/bbg6775/humangenetik/shared/references/GRCh38/
ls -l /beegfs/u/bbg6775/humangenetik/shared/HG002/hifi/
```

### Tool not found
```bash
# Load modules
module load minimap2/2.28
module load samtools/1.20

# Or use conda
conda activate llmap-bench
```

## Data Paths (Hummel)

From `datasets/datasets.yaml`:

| Dataset | Path |
|---------|------|
| GRCh38 reference | `/beegfs/u/bbg6775/humangenetik/shared/references/GRCh38/GRCh38.p14.fa` |
| GENCODE v46 | `/beegfs/u/bbg6775/humangenetik/shared/references/gencode/v46/gencode.v46.transcripts.fa` |
| HG002 HiFi | `/beegfs/u/bbg6775/humangenetik/shared/HG002/hifi/HG002.HiFi.chr14_chr20.fastq.gz` |
| HG002 Illumina R1 | `/beegfs/u/bbg6775/humangenetik/shared/HG002/illumina/HG002.NovaSeq.R1.chr14_chr20.fastq.gz` |
| HG002 Illumina R2 | `/beegfs/u/bbg6775/humangenetik/shared/HG002/illumina/HG002.NovaSeq.R2.chr14_chr20.fastq.gz` |
| HPRC Iso-seq | `/beegfs/u/bbg6775/humangenetik/shared/transcriptome_longread/hprc_isoseq/HG00290/flnc.bam` |
