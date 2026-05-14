#!/usr/bin/env bash
# minimap2 runner — long + short + splice.
set -euo pipefail

TOOL="minimap2"
TOOL_VERSION_CMD="minimap2 --version"

if [[ -n "${READS_R2:-}" ]]; then
  # paired short-read
  ALIGN_CMD="minimap2 -t $THREADS -ax $PRESET $REF $READS_R1 $READS_R2"
else
  ALIGN_CMD="minimap2 -t $THREADS -ax $PRESET $REF $READS_R1"
fi

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
