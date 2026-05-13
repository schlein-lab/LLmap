#!/bin/bash
#SBATCH --job-name=llmap_validate_real
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=llmap_validate_%j.out
#SBATCH --error=llmap_validate_%j.err

# LLmap Real Reference Validation SLURM Script
# Phase 5.5 — GPU validation against hg38 chr14 IGH locus
#
# Usage:
#   sbatch --export=ALL,REF=/path/to/chr14_igh.fa,READS=/path/to/reads.fq \
#          scripts/slurm_validate_real.sh

set -euo pipefail

module load cuda/12.3
module load gcc/13.2

echo "=== LLmap Real Reference Validation ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start: $(date)"
echo ""

# Validate environment
if [[ -z "${REF:-}" ]]; then
    echo "ERROR: REF environment variable not set"
    exit 1
fi

if [[ -z "${READS:-}" ]]; then
    echo "ERROR: READS environment variable not set"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference file not found: $REF"
    exit 1
fi

if [[ ! -f "$READS" ]]; then
    echo "ERROR: Reads file not found: $READS"
    exit 1
fi

# Optional baseline BAM
BASELINE_ARG=""
if [[ -n "${BASELINE:-}" ]] && [[ -f "$BASELINE" ]]; then
    BASELINE_ARG="--baseline $BASELINE"
fi

# Optional ground truth BED
TRUTH_ARG=""
if [[ -n "${TRUTH:-}" ]] && [[ -f "$TRUTH" ]]; then
    TRUTH_ARG="--truth $TRUTH"
fi

# Output directory
OUTPUT_DIR="${OUTPUT_DIR:-./validation_output}"
mkdir -p "$OUTPUT_DIR"

# Work directory
WORKDIR="${WORKDIR:-$(pwd)}"
cd "$WORKDIR"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Configuration:"
echo "  Reference:  $REF"
echo "  Reads:      $READS"
echo "  Baseline:   ${BASELINE:-none}"
echo "  Truth:      ${TRUTH:-none}"
echo "  Output:     $OUTPUT_DIR"
echo "  GPU:        enabled"
echo ""

# Run validation
./build/llmap validate-real \
    --reference "$REF" \
    --reads "$READS" \
    $BASELINE_ARG \
    $TRUTH_ARG \
    --output-dir "$OUTPUT_DIR" \
    --k 15 \
    --w 10 \
    --min-identity 0.70 \
    --gpu \
    --gpu-device 0

EXIT_CODE=$?

echo ""
echo "=== Validation Complete ==="
echo "Exit code: $EXIT_CODE"
echo "End: $(date)"

exit $EXIT_CODE
