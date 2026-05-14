#!/usr/bin/env bash
# Bowtie2 runner — short-read only.
set -euo pipefail

TOOL="bowtie2"
TOOL_VERSION_CMD="bowtie2 --version"

# Bowtie2 needs a pre-built index, not the FASTA
INDEX_PREFIX="${BOWTIE2_INDEX_PREFIX:-${REF%.fa}}"
if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]]; then
  echo "[$TOOL] Bowtie2 index missing for prefix $INDEX_PREFIX — building" >&2
  bowtie2-build --threads "$THREADS" "$REF" "$INDEX_PREFIX"
fi

if [[ -n "${READS_R2:-}" ]]; then
  ALIGN_CMD="bowtie2 --very-sensitive -p $THREADS -x $INDEX_PREFIX -1 $READS_R1 -2 $READS_R2"
else
  ALIGN_CMD="bowtie2 --very-sensitive -p $THREADS -x $INDEX_PREFIX -U $READS_R1"
fi

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
