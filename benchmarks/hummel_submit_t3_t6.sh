#!/usr/bin/env bash
# Generate and optionally submit SLURM jobs for T3-T6 benchmark tasks on Hummel.
#
# Prerequisites:
#   - Run this script on Hummel (hummel-login) or copy the generated scripts there
#   - Ensure tools are installed (run: ./runners/check_versions.sh)
#   - Ensure datasets exist (see datasets/datasets.yaml for paths)
#
# Usage:
#   ./hummel_submit_t3_t6.sh --generate-only   # Generate scripts without submitting
#   ./hummel_submit_t3_t6.sh --submit          # Generate and submit to SLURM
#   ./hummel_submit_t3_t6.sh --status          # Check status of submitted jobs
#
# After submission, monitor with:
#   squeue -u $USER
#   tail -f benchmarks/reports/T3/llmap/rep0/slurm.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPORTS_DIR="$SCRIPT_DIR/reports"
SLURM_DIR="$SCRIPT_DIR/slurm_jobs"
RUNNERS_DIR="$SCRIPT_DIR/runners"

# Configuration
THREADS=16
MEM_GB=64
TIME_HOURS=8
REPLICATES=3
SEED_BASE=42

# Hummel dataset paths (from datasets/datasets.yaml)
REF_GRCH38="/beegfs/u/bbg6775/humangenetik/shared/references/GRCh38/GRCh38.p14.fa"
REF_GENCODE="/beegfs/u/bbg6775/humangenetik/shared/references/gencode/v46/gencode.v46.transcripts.fa"
REF_IGH_LOCUS="$SCRIPT_DIR/datasets/cache/igh_locus.fa"

READS_HG002_HIFI="/beegfs/u/bbg6775/humangenetik/shared/HG002/hifi/HG002.HiFi.chr14_chr20.fastq.gz"
READS_HG002_ILLUMINA_R1="/beegfs/u/bbg6775/humangenetik/shared/HG002/illumina/HG002.NovaSeq.R1.chr14_chr20.fastq.gz"
READS_HG002_ILLUMINA_R2="/beegfs/u/bbg6775/humangenetik/shared/HG002/illumina/HG002.NovaSeq.R2.chr14_chr20.fastq.gz"
READS_HPRC_ISOSEQ="/beegfs/u/bbg6775/humangenetik/shared/transcriptome_longread/hprc_isoseq/HG00290/flnc.bam"

# Partitions
CPU_PARTITION="standard"
GPU_PARTITION="gpu"

MODE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --generate-only) MODE="generate"; shift ;;
    --submit)        MODE="submit"; shift ;;
    --status)        MODE="status"; shift ;;
    -h|--help)
      echo "Usage: $0 [--generate-only|--submit|--status]"
      echo ""
      echo "Options:"
      echo "  --generate-only  Generate sbatch scripts without submitting"
      echo "  --submit         Generate and submit to SLURM"
      echo "  --status         Check status of running jobs"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$MODE" ]]; then
  echo "Please specify --generate-only, --submit, or --status"
  exit 1
fi

if [[ "$MODE" == "status" ]]; then
  echo "=== Benchmark Job Status ==="
  squeue -u "${USER:-$(whoami)}" -o "%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" 2>/dev/null || echo "(squeue not available)"
  echo ""
  echo "=== Completed runs (done.flag exists) ==="
  find "$REPORTS_DIR" -name "done.flag" -type f | while read flag; do
    dir=$(dirname "$flag")
    echo "  ✓ $dir"
  done
  echo ""
  echo "=== Pending runs (no done.flag) ==="
  for task in T3 T4 T5 T6; do
    task_dir="$REPORTS_DIR/$task"
    [[ -d "$task_dir" ]] || continue
    for tool_dir in "$task_dir"/*/; do
      [[ -d "$tool_dir" ]] || continue
      tool=$(basename "$tool_dir")
      for rep in 0 1 2; do
        rep_dir="$tool_dir/rep$rep"
        if [[ -d "$rep_dir" && ! -f "$rep_dir/done.flag" ]]; then
          echo "  ○ $task/$tool/rep$rep"
        fi
      done
    done
  done
  exit 0
fi

mkdir -p "$SLURM_DIR"

# T3-T6 matrix: task tool preset ref reads_r1 [reads_r2] partition gpu
MATRIX=(
  # T3: Real HiFi WGS (HG002)
  "T3 llmap map-hifi $REF_GRCH38 $READS_HG002_HIFI - $GPU_PARTITION 1"
  "T3 minimap2 map-hifi $REF_GRCH38 $READS_HG002_HIFI - $CPU_PARTITION 0"
  "T3 winnowmap2 map-pb $REF_GRCH38 $READS_HG002_HIFI - $CPU_PARTITION 0"

  # T4: Real Illumina WGS (HG002)
  "T4 llmap sr $REF_GRCH38 $READS_HG002_ILLUMINA_R1 $READS_HG002_ILLUMINA_R2 $GPU_PARTITION 1"
  "T4 minimap2 sr $REF_GRCH38 $READS_HG002_ILLUMINA_R1 $READS_HG002_ILLUMINA_R2 $CPU_PARTITION 0"
  "T4 bwa_mem2 - $REF_GRCH38 $READS_HG002_ILLUMINA_R1 $READS_HG002_ILLUMINA_R2 $CPU_PARTITION 0"
  "T4 bowtie2 - $REF_GRCH38 $READS_HG002_ILLUMINA_R1 $READS_HG002_ILLUMINA_R2 $CPU_PARTITION 0"

  # T5: Iso-seq transcriptome (HPRC HG00290)
  "T5 llmap splice $REF_GENCODE $READS_HPRC_ISOSEQ - $GPU_PARTITION 1"
  "T5 minimap2 splice $REF_GENCODE $READS_HPRC_ISOSEQ - $CPU_PARTITION 0"
  "T5 star - $REF_GENCODE $READS_HPRC_ISOSEQ - $CPU_PARTITION 0"

  # T6: Targeted IGH-locus (same reads as T3, but IGH-locus reference)
  "T6 llmap map-hifi $REF_IGH_LOCUS $READS_HG002_HIFI - $GPU_PARTITION 1"
  "T6 minimap2 map-hifi $REF_IGH_LOCUS $READS_HG002_HIFI - $CPU_PARTITION 0"
  "T6 winnowmap2 map-pb $REF_IGH_LOCUS $READS_HG002_HIFI - $CPU_PARTITION 0"
)

GENERATED=0
SKIPPED=0
SUBMITTED=0

generate_sbatch() {
  local task=$1 tool=$2 preset=$3 ref=$4 reads_r1=$5 reads_r2=$6 partition=$7 use_gpu=$8 rep=$9

  local out_dir="$REPORTS_DIR/$task/$tool/rep$rep"
  local script_path="$SLURM_DIR/${task}_${tool}_rep${rep}.sbatch"
  local seed=$((SEED_BASE + rep))

  # Skip if already done
  if [[ -f "$out_dir/done.flag" ]]; then
    echo "[SKIP] $task/$tool/rep$rep - already complete"
    ((SKIPPED++)) || true
    return
  fi

  mkdir -p "$out_dir"

  local gpu_directive=""
  if [[ "$use_gpu" == "1" ]]; then
    gpu_directive="#SBATCH --gres=gpu:1"
  fi

  local reads_r2_export=""
  if [[ "$reads_r2" != "-" ]]; then
    reads_r2_export="export READS_R2=\"$reads_r2\""
  fi

  cat > "$script_path" <<SBATCH
#!/bin/bash
#SBATCH --job-name=bench_${task}_${tool}_${rep}
#SBATCH --output=${out_dir}/slurm.out
#SBATCH --error=${out_dir}/slurm.err
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM_GB}G
#SBATCH --time=${TIME_HOURS}:00:00
#SBATCH --partition=${partition}
${gpu_directive}

set -euo pipefail

echo "=== LLmap Benchmark: ${task}/${tool}/rep${rep} ==="
echo "Started: \$(date -Iseconds)"
echo "Host: \$(hostname)"
echo "Partition: \$SLURM_JOB_PARTITION"

# Environment setup
module load samtools/1.20 2>/dev/null || true

export TOOL="${tool}"
export TASK="${task}"
export REPLICATE="${rep}"
export PRESET="${preset}"
export REF="${ref}"
export READS_R1="${reads_r1}"
${reads_r2_export}
export THREADS="${THREADS}"
export SEED="${seed}"
export OUTPUT_DIR="${out_dir}"

# Verify inputs exist
echo "Checking inputs..."
[[ -f "\$REF" ]] || { echo "ERROR: Reference not found: \$REF"; exit 1; }
[[ -f "\$READS_R1" ]] || { echo "ERROR: Reads not found: \$READS_R1"; exit 1; }
if [[ -n "\${READS_R2:-}" ]]; then
  [[ -f "\$READS_R2" ]] || { echo "ERROR: R2 reads not found: \$READS_R2"; exit 1; }
fi

# Run the tool-specific runner
echo "Running ${tool}..."
cd "${LLMAP_ROOT}"
bash "${RUNNERS_DIR}/run_${tool}.sh"

echo "Completed: \$(date -Iseconds)"
SBATCH

  chmod +x "$script_path"
  echo "[GENERATED] $script_path"
  ((GENERATED++)) || true
}

echo "=== Generating SLURM scripts for T3-T6 ==="
echo "Output directory: $SLURM_DIR"
echo ""

for entry in "${MATRIX[@]}"; do
  read -r task tool preset ref reads_r1 reads_r2 partition use_gpu <<< "$entry"
  for rep in $(seq 0 $((REPLICATES - 1))); do
    generate_sbatch "$task" "$tool" "$preset" "$ref" "$reads_r1" "$reads_r2" "$partition" "$use_gpu" "$rep"
  done
done

echo ""
echo "=== Summary ==="
echo "Generated: $GENERATED scripts"
echo "Skipped (already done): $SKIPPED"
echo "Scripts location: $SLURM_DIR/"

if [[ "$MODE" == "submit" ]]; then
  echo ""
  echo "=== Submitting to SLURM ==="

  if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on Hummel?"
    echo "Copy the scripts to Hummel and run sbatch manually:"
    echo "  scp -r $SLURM_DIR hummel-login:~/llmap/benchmarks/"
    exit 1
  fi

  for script in "$SLURM_DIR"/*.sbatch; do
    [[ -f "$script" ]] || continue
    name=$(basename "$script" .sbatch)

    # Check if this job is already done
    task=$(echo "$name" | cut -d_ -f1)
    tool=$(echo "$name" | cut -d_ -f2)
    rep=$(echo "$name" | sed 's/.*rep//')

    if [[ -f "$REPORTS_DIR/$task/$tool/rep$rep/done.flag" ]]; then
      continue
    fi

    echo "Submitting: $script"
    sbatch "$script"
    ((SUBMITTED++)) || true
  done

  echo ""
  echo "Submitted: $SUBMITTED jobs"
  echo ""
  echo "Monitor with:"
  echo "  squeue -u \$USER"
  echo "  $0 --status"
else
  echo ""
  echo "=== Next steps ==="
  echo "1. Copy scripts to Hummel (if running locally):"
  echo "     scp -r $SLURM_DIR hummel-login:~/llmap/benchmarks/"
  echo ""
  echo "2. On Hummel, submit all jobs:"
  echo "     for script in ~/llmap/benchmarks/slurm_jobs/*.sbatch; do sbatch \"\$script\"; done"
  echo ""
  echo "   Or use the orchestrator:"
  echo "     python orchestrate.py --hummel --tasks T3 T4 T5 T6"
  echo ""
  echo "3. Monitor progress:"
  echo "     squeue -u \$USER"
  echo "     $0 --status"
fi
