#!/usr/bin/env bash
# Orchestrates the full benchmark matrix.
# Runs idempotently — skips (tool, task, replicate) cells where done.flag exists.
#
# Local-host mode: runs synthetic tasks T1, T2 inline.
# Hummel mode: submits each cell as a sbatch job under benchmarks/slurm/.
#
# Usage:
#   ./submit_all.sh                       # auto-detect (sbatch present → Hummel)
#   ./submit_all.sh --local               # force local
#   ./submit_all.sh --task T1 --task T2   # subset
#   ./submit_all.sh --dry-run             # print, do not execute

set -euo pipefail

LLMAP_ROOT="${LLMAP_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
BENCH_ROOT="$LLMAP_ROOT/benchmarks"
RUNNERS="$BENCH_ROOT/runners"
REPORTS="$BENCH_ROOT/reports"
DATASETS="$BENCH_ROOT/datasets"

# Defaults pulled from datasets/tools.yaml
THREADS=16
REPLICATES=3
SEED_BASE=42

TASKS_FILTER=()
MODE="auto"
DRY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local)   MODE="local"; shift ;;
    --hummel)  MODE="hummel"; shift ;;
    --dry-run) DRY=1; shift ;;
    --task)    TASKS_FILTER+=("$2"); shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Auto-detect mode
if [[ "$MODE" == "auto" ]]; then
  if command -v sbatch >/dev/null 2>&1; then
    MODE="hummel"
  else
    MODE="local"
  fi
fi

# Matrix definition. Each row:
#   TASK TOOL PRESET REF READS_R1 [READS_R2]
read -r -d '' MATRIX <<'EOF' || true
T1 llmap       map-hifi  synth_t1 fastq
T1 minimap2    map-hifi  synth_t1 fastq
T1 winnowmap2  map-pb    synth_t1 fastq
T2 llmap       map-hifi  synth_t2 fastq
T2 minimap2    map-hifi  synth_t2 fastq
T2 winnowmap2  map-pb    synth_t2 fastq
T3 llmap       map-hifi  grch38   hg002_hifi_chr14_chr20
T3 minimap2    map-hifi  grch38   hg002_hifi_chr14_chr20
T3 winnowmap2  map-pb    grch38   hg002_hifi_chr14_chr20
T4 llmap       sr        grch38   hg002_illumina_chr14_chr20
T4 minimap2    sr        grch38   hg002_illumina_chr14_chr20
T4 bwa_mem2    -         grch38   hg002_illumina_chr14_chr20
T4 bowtie2     -         grch38   hg002_illumina_chr14_chr20
T5 llmap       splice    gencode_v46 hprc_isoseq_hg00290
T5 minimap2    splice    gencode_v46 hprc_isoseq_hg00290
T5 star        -         gencode_v46 hprc_isoseq_hg00290
T6 llmap       map-hifi  igh_locus hg002_hifi_chr14_chr20
T6 minimap2    map-hifi  igh_locus hg002_hifi_chr14_chr20
T6 winnowmap2  map-pb    igh_locus hg002_hifi_chr14_chr20
EOF

run_cell() {
  local task=$1 tool=$2 preset=$3 ref_id=$4 reads_id=$5

  if [[ ${#TASKS_FILTER[@]} -gt 0 ]]; then
    local found=0
    for t in "${TASKS_FILTER[@]}"; do [[ "$t" == "$task" ]] && found=1; done
    [[ $found -eq 1 ]] || return 0
  fi

  for rep in $(seq 0 $((REPLICATES - 1))); do
    local out="$REPORTS/$task/$tool/rep$rep"
    if [[ -f "$out/done.flag" ]]; then
      [[ $DRY -eq 0 ]] && echo "[$task/$tool/rep$rep] already done, skip"
      continue
    fi
    mkdir -p "$out"

    local seed=$((SEED_BASE + rep))
    local cmd="TOOL=$tool TASK=$task REPLICATE=$rep PRESET=$preset \
              REF_ID=$ref_id READS_ID=$reads_id THREADS=$THREADS SEED=$seed \
              OUTPUT_DIR=$out bash $RUNNERS/run_${tool}.sh"

    if [[ $DRY -eq 1 ]]; then
      echo "[DRY] $cmd"
      continue
    fi

    case "$MODE" in
      local)
        echo "[local] $task/$tool/rep$rep"
        eval "$cmd"
        ;;
      hummel)
        local sbatch_script="$out/submit.sbatch"
        cat > "$sbatch_script" <<SBATCH
#!/bin/bash
#SBATCH --job-name=bench_${task}_${tool}_${rep}
#SBATCH --output=$out/slurm.out
#SBATCH --error=$out/slurm.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=64G
#SBATCH --time=4:00:00
$cmd
SBATCH
        sbatch "$sbatch_script"
        ;;
    esac
  done
}

while read -r task tool preset ref_id reads_id; do
  [[ -z "${task:-}" || "${task:0:1}" == "#" ]] && continue
  run_cell "$task" "$tool" "$preset" "$ref_id" "$reads_id"
done <<< "$MATRIX"
