#!/usr/bin/env bash
# Run T1 and T2 synthetic benchmarks locally.
# Usage: ./run_local_synthetic.sh [--dry-run] [--tools TOOLS] [--replicates N]
#
# This script sets up environment variables and calls the runner scripts
# for locally-generated synthetic datasets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LLMAP_ROOT="$(dirname "$SCRIPT_DIR")"
DATASETS_CACHE="$SCRIPT_DIR/datasets/cache"
REPORTS_DIR="$SCRIPT_DIR/reports"
RUNNERS_DIR="$SCRIPT_DIR/runners"

DRY_RUN=false
TOOLS="llmap,minimap2"  # Tools to run (comma-separated)
REPLICATES=1            # Number of replicates (1 for quick local test)
THREADS=4               # Local CPU threads

while [[ $# -gt 0 ]]; do
  case $1 in
    --dry-run) DRY_RUN=true; shift ;;
    --tools) TOOLS="$2"; shift 2 ;;
    --replicates) REPLICATES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Check prerequisites
if ! command -v samtools &>/dev/null; then
  echo "ERROR: samtools not found" >&2
  exit 1
fi

# Dataset paths
declare -A TASK_REF
declare -A TASK_READS
declare -A TASK_PRESET

TASK_REF[T1]="$DATASETS_CACHE/synth_t1/reference.fa"
TASK_READS[T1]="$DATASETS_CACHE/synth_t1/reads.fastq"
TASK_PRESET[T1]="map-hifi"

TASK_REF[T2]="$DATASETS_CACHE/synth_t2/reference.fa"
TASK_READS[T2]="$DATASETS_CACHE/synth_t2/reads.fastq"
TASK_PRESET[T2]="map-hifi"

# Check datasets exist
for task in T1 T2; do
  if [[ ! -f "${TASK_REF[$task]}" ]]; then
    echo "ERROR: Reference not found: ${TASK_REF[$task]}" >&2
    echo "Run: llmap generate-synth --task $(echo $task | tr '[:upper:]' '[:lower:]') --output $DATASETS_CACHE/synth_$(echo $task | tr '[:upper:]' '[:lower:]')" >&2
    exit 1
  fi
  if [[ ! -f "${TASK_READS[$task]}" ]]; then
    echo "ERROR: Reads not found: ${TASK_READS[$task]}" >&2
    exit 1
  fi
done

run_benchmark() {
  local task="$1"
  local tool="$2"
  local rep="$3"

  local output_dir="$REPORTS_DIR/$task/$tool/rep$rep"

  # Skip if already done
  if [[ -f "$output_dir/done.flag" ]]; then
    echo "[SKIP] $task/$tool/rep$rep already complete"
    return 0
  fi

  # Check tool availability
  local tool_bin=""
  case "$tool" in
    llmap)
      tool_bin="$LLMAP_ROOT/build/src/llmap"
      if [[ ! -x "$tool_bin" ]]; then
        echo "[SKIP] $task/$tool: binary not found at $tool_bin"
        return 0
      fi
      ;;
    minimap2)
      if ! command -v minimap2 &>/dev/null; then
        echo "[SKIP] $task/$tool: minimap2 not in PATH"
        return 0
      fi
      ;;
    winnowmap2)
      if ! command -v winnowmap &>/dev/null && ! command -v winnowmap2 &>/dev/null; then
        echo "[SKIP] $task/$tool: winnowmap2 not in PATH"
        return 0
      fi
      ;;
    *)
      echo "[SKIP] $task/$tool: unsupported for local run"
      return 0
      ;;
  esac

  local ref_path="${TASK_REF[$task]}"
  local reads_path="${TASK_READS[$task]}"
  local preset="${TASK_PRESET[$task]}"

  echo "[RUN] $task/$tool/rep$rep"

  if $DRY_RUN; then
    echo "  REF=$ref_path"
    echo "  READS_R1=$reads_path"
    echo "  PRESET=$preset"
    echo "  OUTPUT_DIR=$output_dir"
    return 0
  fi

  mkdir -p "$output_dir"

  # Export environment for runner
  export TOOL="$tool"
  export TASK="$task"
  export REPLICATE="$rep"
  export REF="$ref_path"
  export READS_R1="$reads_path"
  export READS_R2=""
  export PRESET="$preset"
  export THREADS="$THREADS"
  export OUTPUT_DIR="$output_dir"
  export SEED="$((42 + rep))"
  export LLMAP_BIN="$LLMAP_ROOT/build/src/llmap"

  # Run the tool-specific runner
  local runner="$RUNNERS_DIR/run_${tool}.sh"
  if [[ -f "$runner" ]]; then
    if bash "$runner"; then
      echo "[DONE] $task/$tool/rep$rep"
    else
      echo "[FAIL] $task/$tool/rep$rep"
      return 1
    fi
  else
    echo "[SKIP] $task/$tool: no runner at $runner"
  fi
}

# Main execution
echo "=== LLmap Local Synthetic Benchmark ==="
echo "Tasks: T1, T2"
echo "Tools: $TOOLS"
echo "Replicates: $REPLICATES"
echo "Threads: $THREADS"
echo "Dry run: $DRY_RUN"
echo ""

IFS=',' read -ra TOOL_LIST <<< "$TOOLS"

total=0
completed=0
skipped=0
failed=0

for task in T1 T2; do
  for tool in "${TOOL_LIST[@]}"; do
    for rep in $(seq 0 $((REPLICATES - 1))); do
      ((total++)) || true
      if run_benchmark "$task" "$tool" "$rep"; then
        if [[ -f "$REPORTS_DIR/$task/$tool/rep$rep/done.flag" ]] || $DRY_RUN; then
          ((completed++)) || true
        else
          ((skipped++)) || true
        fi
      else
        ((failed++)) || true
      fi
    done
  done
done

echo ""
echo "=== Summary ==="
echo "Total cells: $total"
echo "Completed: $completed"
echo "Skipped: $skipped"
echo "Failed: $failed"
