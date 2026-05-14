#!/usr/bin/env bash
# BWA-MEM2 runner — short-read only.
set -euo pipefail

TOOL="bwa_mem2"
TOOL_VERSION_CMD="bwa-mem2 version"

# Index check (creates if absent)
if [[ ! -f "${REF}.bwt.2bit.64" ]]; then
  echo "[$TOOL] bwa-mem2 index missing for $REF — building" >&2
  bwa-mem2 index "$REF"
fi

if [[ -n "${READS_R2:-}" ]]; then
  ALIGN_CMD="bwa-mem2 mem -t $THREADS $REF $READS_R1 $READS_R2"
else
  ALIGN_CMD="bwa-mem2 mem -t $THREADS $REF $READS_R1"
fi

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
