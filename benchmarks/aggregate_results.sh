#!/usr/bin/env bash
# Aggregate benchmark results and generate reports.
# Computes per-run metrics, generates per-task reports.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
METRICS_DIR="$SCRIPT_DIR/metrics"
REPORTS_DIR="$SCRIPT_DIR/reports"
DATASETS_CACHE="$SCRIPT_DIR/datasets/cache"

echo "=== Aggregating Benchmark Results ==="

# Count reads per task for normalization
count_reads() {
  local fastq="$1"
  # Count lines that start with @ (FASTQ headers)
  grep -c '^@' "$fastq" 2>/dev/null || echo 0
}

# Task T1 - synthetic WGS
T1_FASTQ="$DATASETS_CACHE/synth_t1/reads.fastq"
T1_TRUTH="$DATASETS_CACHE/synth_t1/truth.tsv"
T1_TOTAL=$(count_reads "$T1_FASTQ")

# Task T2 - synthetic paralog stress
T2_FASTQ="$DATASETS_CACHE/synth_t2/reads.fastq"
T2_TRUTH="$DATASETS_CACHE/synth_t2/truth.tsv"
T2_TOTAL=$(count_reads "$T2_FASTQ")

echo "T1 total reads: $T1_TOTAL"
echo "T2 total reads: $T2_TOTAL"

compute_metrics() {
  local task="$1"
  local tool="$2"
  local rep="$3"
  local total="$4"
  local truth="$5"

  local run_dir="$REPORTS_DIR/$task/$tool/rep$rep"
  local bam="$run_dir/alignments.bam"

  if [[ ! -f "$bam" ]]; then
    echo "[SKIP] $task/$tool/rep$rep: no BAM file"
    return 0
  fi

  if [[ -f "$run_dir/mapping_summary.json" ]]; then
    echo "[SKIP] $task/$tool/rep$rep: metrics already computed"
    return 0
  fi

  echo "[METRICS] $task/$tool/rep$rep"

  local truth_arg=""
  if [[ -f "$truth" ]]; then
    truth_arg="--truth $truth"
  fi

  python3 "$METRICS_DIR/compute.py" \
    --bam "$bam" \
    --out-dir "$run_dir" \
    --total "$total" \
    $truth_arg
}

# Compute metrics for all runs
for task_dir in "$REPORTS_DIR"/T*/; do
  task=$(basename "$task_dir")
  for tool_dir in "$task_dir"/*/; do
    [[ -d "$tool_dir" ]] || continue
    tool=$(basename "$tool_dir")
    for rep_dir in "$tool_dir"/rep*/; do
      [[ -d "$rep_dir" ]] || continue
      rep=$(basename "$rep_dir" | sed 's/rep//')

      case "$task" in
        T1) compute_metrics "$task" "$tool" "$rep" "$T1_TOTAL" "$T1_TRUTH" ;;
        T2) compute_metrics "$task" "$tool" "$rep" "$T2_TOTAL" "$T2_TRUTH" ;;
        *)  echo "[SKIP] Unknown task: $task" ;;
      esac
    done
  done
done

echo ""
echo "=== Generating Reports ==="

# Generate reports
echo "[REPORT] Generating for T1 T2"
python3 "$SCRIPT_DIR/report.py" --tasks T1 T2 --output-dir "$REPORTS_DIR" --verbose --no-plots

echo ""
echo "=== Summary ==="
for task in T1 T2; do
  if [[ -f "$REPORTS_DIR/$task/README.md" ]]; then
    echo "$task: README.md generated"
  else
    echo "$task: no report"
  fi
done
