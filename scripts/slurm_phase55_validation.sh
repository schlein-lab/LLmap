#!/bin/bash
#SBATCH --job-name=llmap_p55_gpu_val
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=llmap_p55_%j.out
#SBATCH --error=llmap_p55_%j.err

# LLmap Phase 5.5 — Full GPU Validation Pipeline
#
# This script:
# 1. Builds LLmap with CUDA enabled
# 2. Generates synthetic IGH locus data
# 3. Runs minimap2 baseline alignment
# 4. Runs LLmap validation
# 5. Reports kill-switch verdict
#
# Usage:
#   sbatch scripts/slurm_phase55_validation.sh

set -euo pipefail

echo "=== LLmap Phase 5.5 GPU Validation ==="
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURM_NODELIST"
echo "Start:     $(date)"
echo ""

# Load modules
module load cuda/12.3 2>/dev/null || module load cuda/12.0 2>/dev/null || echo "CUDA module load warning"
module load gcc/13.2 2>/dev/null || module load gcc/12.2 2>/dev/null || echo "GCC module load warning"
module load cmake/3.28 2>/dev/null || module load cmake/3.20 2>/dev/null || echo "CMake module load warning"

# Report environment
echo "=== Environment ==="
echo "CUDA:  $(which nvcc 2>/dev/null || echo 'not found')"
echo "GCC:   $(which g++ 2>/dev/null || echo 'not found')"
echo "CMake: $(which cmake 2>/dev/null || echo 'not found')"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "GPU info unavailable"
echo ""

WORKDIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$WORKDIR"

DATADIR="/beegfs/u/bbg6775/llmap/validation_p55_${SLURM_JOB_ID}"
mkdir -p "$DATADIR"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Step 1: Build with CUDA
echo "=== Step 1: Building LLmap with CUDA ==="
rm -rf build_cuda
mkdir -p build_cuda
cd build_cuda

cmake -DLLMAP_ENABLE_CUDA=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DLLMAP_USE_NATIVE_ARCH=ON \
      ..

cmake --build . -j${SLURM_CPUS_PER_TASK:-16}

if [[ ! -x src/llmap ]]; then
    echo "ERROR: Build failed — llmap executable not found"
    exit 1
fi
echo "Build successful: $(./src/llmap --version)"
echo ""

# Step 2: Generate synthetic data
echo "=== Step 2: Generating Synthetic Data ==="
./src/llmap generate-synth \
    --out-prefix "$DATADIR/synth" \
    --preset balanced \
    --seed 42

ls -la "$DATADIR"/synth_*
echo ""

# Step 3: Run minimap2 baseline (if available)
echo "=== Step 3: Running Minimap2 Baseline ==="
if command -v minimap2 &>/dev/null; then
    minimap2 -ax map-hifi \
        -t ${SLURM_CPUS_PER_TASK:-16} \
        "$DATADIR/synth_reference.fa" \
        "$DATADIR/synth_reads.fq" \
        > "$DATADIR/minimap2.sam" 2>/dev/null

    # Convert to simple mapped/unmapped
    awk 'BEGIN{OFS="\t"} !/^@/{print $1, ($2==4 ? "unmapped" : "mapped")}' \
        "$DATADIR/minimap2.sam" > "$DATADIR/minimap2_status.tsv"

    MM2_MAPPED=$(grep -c "mapped" "$DATADIR/minimap2_status.tsv" || echo "0")
    MM2_TOTAL=$(wc -l < "$DATADIR/minimap2_status.tsv")
    echo "Minimap2: $MM2_MAPPED / $MM2_TOTAL reads mapped"
    BASELINE_ARG="--baseline $DATADIR/minimap2.sam"
else
    echo "minimap2 not found — skipping baseline comparison"
    BASELINE_ARG=""
fi
echo ""

# Step 4: Run LLmap validation
echo "=== Step 4: Running LLmap Validation ==="
./src/llmap validate-real \
    --reference "$DATADIR/synth_reference.fa" \
    --reads "$DATADIR/synth_reads.fq" \
    $BASELINE_ARG \
    --k 15 \
    --w 10 \
    --min-identity 0.70 \
    --gpu \
    --gpu-device 0 \
    | tee "$DATADIR/validation_result.txt"

EXIT_CODE=$?
echo ""

# Report
echo "=== Summary ==="
echo "Output dir:  $DATADIR"
echo "Exit code:   $EXIT_CODE"
echo "End:         $(date)"
echo ""

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "VERDICT: PASS — Kill-switch validation succeeded"
else
    echo "VERDICT: FAIL — Kill-switch validation failed"
fi

exit $EXIT_CODE
