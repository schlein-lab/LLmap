#!/bin/bash
# LLmap — the HPC cluster CPU-only build + Transcript-Mode test driver.
#
# Mirrors the local CPU-only baseline (Agent 2): no CUDA/FAISS/ONNX/Claude/bench.
# Run from the repo root on a the HPC cluster compute/front node (NOT login — build is CPU
# heavy). Builds, runs the new Transcript-Mode unit tests, and (if FLNC data +
# reference are given) does a real `--mode transcript` smoke run.
#
# Usage:
#   scripts/build_test_the HPC cluster.sh [BUILD_DIR] [FLNC_FASTQ] [REFERENCE_FA]
#
# Examples:
#   scripts/build_test_the HPC cluster.sh build
#   scripts/build_test_the HPC cluster.sh build /path/flnc.fastq /path/ref.fa
set -euo pipefail

BUILD_DIR="${1:-build}"
FLNC="${2:-}"
REF="${3:-}"
JOBS="$(nproc)"

echo "== [1/4] configure (CPU-only) =="
cmake -S . -B "$BUILD_DIR" \
  -DLLMAP_ENABLE_CUDA=OFF \
  -DLLMAP_ENABLE_FOUNDATION=OFF \
  -DLLMAP_ENABLE_FAISS=OFF \
  -DLLMAP_ENABLE_CLAUDE=OFF \
  -DLLMAP_ENABLE_BENCH=OFF

echo "== [2/4] build (-j$JOBS) =="
cmake --build "$BUILD_DIR" -j"$JOBS"

echo "== [3/4] Transcript-Mode unit tests =="
ctest --test-dir "$BUILD_DIR" -R 'input_sniffer|transcript_stage' --output-on-failure

echo "== [3b] full ctest (optional, comment out if too heavy) =="
ctest --test-dir "$BUILD_DIR" --output-on-failure || {
  echo "WARN: full ctest had failures — inspect above"; }

LLMAP="$BUILD_DIR/src/llmap"
echo "== llmap binary: $LLMAP =="
"$LLMAP" --version || true
"$LLMAP" align --help | grep -A2 -- '--mode' || true

if [[ -n "$FLNC" && -n "$REF" ]]; then
  echo "== [4/4] real FLNC --mode transcript smoke run =="
  OUT="$BUILD_DIR/transcript_smoke.bam"
  "$LLMAP" align --mode transcript -r "$FLNC" --reference "$REF" \
    -o "$OUT" --bam -v 2>&1 | tee "$BUILD_DIR/transcript_smoke.log" | grep -iE 'mode-detect|transcript' || true
  echo "-- spliced (N-op) CIGARs in output (first 5) --"
  if command -v samtools >/dev/null 2>&1; then
    samtools view "$OUT" 2>/dev/null | awk '$6 ~ /N/ {print $1, $4, $6}' | head -5
    echo "spliced reads: $(samtools view "$OUT" 2>/dev/null | awk '$6 ~ /N/' | wc -l)"
  else
    echo "(samtools not found — inspect $OUT manually)"
  fi
else
  echo "== [4/4] skipped real FLNC run (no FLNC_FASTQ + REFERENCE_FA given) =="
fi

echo "== DONE =="
