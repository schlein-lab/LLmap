#!/usr/bin/env bash
# LLmap runner.
# Writes BOTH the BAM (for like-for-like comparison) and the Parquet
# probabilistic output (LLmap-only metric: posterior calibration).
set -euo pipefail

TOOL="llmap"
TOOL_VERSION_CMD="${LLMAP_BIN:-llmap} --version"
PARQUET_OUT="${OUTPUT_DIR}/probabilities.parquet"

# Build / use index
INDEX="${OUTPUT_DIR}/ref.llmap.idx"
if [[ ! -f "$INDEX" ]]; then
  "${LLMAP_BIN:-llmap}" index -r "$REF" -o "$INDEX" -k 19 -w 19
fi

if [[ -n "${READS_R2:-}" ]]; then
  ALIGN_CMD="${LLMAP_BIN:-llmap} align -i $INDEX -r $READS_R1 -2 $READS_R2 -x $PRESET --llm off --threads $THREADS --bam - --parquet $PARQUET_OUT"
else
  ALIGN_CMD="${LLMAP_BIN:-llmap} align -i $INDEX -r $READS_R1 -x $PRESET --llm off --threads $THREADS --bam - --parquet $PARQUET_OUT"
fi

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
