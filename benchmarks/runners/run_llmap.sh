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

# LLmap align: -x is the READ-TYPE PRESET (map-hifi/map-ont/sr), --reference is the FASTA.
# Earlier versions of this runner passed `-x $REF` which silently broke under the new
# arg parser (preset name "..reference.fa" → rejected) — fixed 2026-05-15.
# Output SAM to file (llmap status goes to stdout, so we must output to file first)
SAM_OUT="${OUTPUT_DIR}/alignments_raw.sam"
LLMAP_PRESET="${LLMAP_PRESET:-map-hifi}"
if [[ -n "${READS_R2:-}" ]]; then
  ALIGN_CMD="${LLMAP_BIN:-llmap} align -x $LLMAP_PRESET --reference $REF -r $READS_R1 -o $SAM_OUT --threads $THREADS --sam --parquet $PARQUET_OUT 1>&2 && cat $SAM_OUT && rm -f $SAM_OUT"
else
  ALIGN_CMD="${LLMAP_BIN:-llmap} align -x $LLMAP_PRESET --reference $REF -r $READS_R1 -o $SAM_OUT --threads $THREADS --sam --parquet $PARQUET_OUT 1>&2 && cat $SAM_OUT && rm -f $SAM_OUT"
fi

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
