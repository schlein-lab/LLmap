#!/usr/bin/env bash
# STAR runner — RNA-seq spliced alignment only (T5).
set -euo pipefail

TOOL="star"
TOOL_VERSION_CMD="STAR --version"

# STAR needs a pre-built genome dir, not a FASTA
GENOME_DIR="${STAR_GENOME_DIR:?STAR_GENOME_DIR not set}"
if [[ ! -d "$GENOME_DIR" || ! -f "$GENOME_DIR/SA" ]]; then
  echo "[$TOOL] STAR genome directory missing or incomplete at $GENOME_DIR" >&2
  exit 2
fi

# STAR writes its outputs into a directory; we collect just the SAM and
# emit it to stdout for the template's BAM-sort step.
TMP_STAR="${OUTPUT_DIR}/star_tmp"
mkdir -p "$TMP_STAR"

ALIGN_CMD="STAR --runThreadN $THREADS --genomeDir $GENOME_DIR --readFilesIn $READS_R1 \
  --readFilesCommand 'zcat -f' --outSAMtype SAM --outSAMunmapped Within \
  --outFileNamePrefix ${TMP_STAR}/ --outStd SAM"

export TOOL TOOL_VERSION_CMD ALIGN_CMD
exec bash "$(dirname "$0")/_template.sh"
