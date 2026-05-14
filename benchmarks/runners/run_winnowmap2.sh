#!/usr/bin/env bash
# Winnowmap2 runner — long-read repeat-aware.
# Requires precomputed repetitive k-mer list (see make_winnowmap_kmers.sh).
set -euo pipefail

TOOL="winnowmap2"
TOOL_VERSION_CMD="winnowmap --version"

REPETITIVE_K15="${REPETITIVE_K15:-${OUTPUT_DIR}/../../datasets/cache/repetitive_k15.txt}"
if [[ ! -f "$REPETITIVE_K15" ]]; then
  echo "[$TOOL] ERROR: repetitive k-mer file not found at $REPETITIVE_K15" >&2
  echo "         Run benchmarks/datasets/make_winnowmap_kmers.sh first."   >&2
  exit 2
fi

ALIGN_CMD="winnowmap -W $REPETITIVE_K15 -t $THREADS -ax $PRESET $REF $READS_R1"

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
