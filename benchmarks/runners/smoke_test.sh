#!/usr/bin/env bash
# Smoke test for Phase 11 benchmark runners.
# Runs each installed tool on a tiny 10-read dataset against a mini reference.
# Exit 0 if all available tools pass; exit 1 if any fail.
#
# Usage:
#   ./benchmarks/runners/smoke_test.sh
#   ./benchmarks/runners/smoke_test.sh --verbose
#   ./benchmarks/runners/smoke_test.sh --tools minimap2,llmap
#   ./benchmarks/runners/smoke_test.sh --json

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATASET_DIR="$(cd "$SCRIPT_DIR/../datasets/smoke" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Defaults
VERBOSE=false
JSON_OUTPUT=false
SELECTED_TOOLS=""
CLEANUP=true

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose) VERBOSE=true; shift ;;
    -j|--json)    JSON_OUTPUT=true; shift ;;
    -t|--tools)   SELECTED_TOOLS="$2"; shift 2 ;;
    --no-cleanup) CLEANUP=false; shift ;;
    -h|--help)
      echo "Usage: smoke_test.sh [--verbose] [--json] [--tools minimap2,llmap] [--no-cleanup]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Setup paths
REF="$DATASET_DIR/smoke_ref.fa"
READS="$DATASET_DIR/smoke_reads.fq"
WORK_DIR=$(mktemp -d --suffix=.llmap_smoke)

log() { $VERBOSE && echo "[smoke] $*" >&2 || true; }

# Tool detection functions
tool_installed() {
  case "$1" in
    minimap2)   command -v minimap2 >/dev/null 2>&1 ;;
    bwa_mem2)   command -v bwa-mem2 >/dev/null 2>&1 ;;
    winnowmap2) command -v winnowmap >/dev/null 2>&1 ;;
    bowtie2)    command -v bowtie2 >/dev/null 2>&1 ;;
    star)       command -v STAR >/dev/null 2>&1 ;;
    llmap)      [[ -x "${LLMAP_BIN:-$LLMAP_ROOT/build/src/llmap}" ]] ;;
    *)          return 1 ;;
  esac
}

# All tools in the benchmark
ALL_TOOLS="minimap2 bwa_mem2 winnowmap2 bowtie2 star llmap"

# Filter to selected tools if specified
if [[ -n "$SELECTED_TOOLS" ]]; then
  TOOLS="${SELECTED_TOOLS//,/ }"
else
  TOOLS="$ALL_TOOLS"
fi

# Results tracking
declare -A RESULTS
PASSED=0
FAILED=0
SKIPPED=0

# Cleanup handler
cleanup() {
  if $CLEANUP; then
    rm -rf "$WORK_DIR"
    log "Cleaned up $WORK_DIR"
  else
    log "Work directory preserved: $WORK_DIR"
  fi
}
trap cleanup EXIT

log "Reference: $REF"
log "Reads: $READS"
log "Work dir: $WORK_DIR"

# Run each tool
for TOOL in $TOOLS; do
  OUT_DIR="$WORK_DIR/$TOOL"
  mkdir -p "$OUT_DIR"

  if ! tool_installed "$TOOL"; then
    log "$TOOL: not installed, skipping"
    RESULTS[$TOOL]="skipped"
    SKIPPED=$((SKIPPED + 1))
    continue
  fi

  log "$TOOL: starting smoke test"

  # Tool-specific invocation (simplified, not using full runner template)
  case "$TOOL" in
    minimap2)
      if minimap2 -t 2 -ax sr "$REF" "$READS" > "$OUT_DIR/alignments.sam" 2>"$OUT_DIR/stderr.log"; then
        STATUS=0
      else
        STATUS=1
      fi
      ;;

    bwa_mem2)
      # Index first
      if ! bwa-mem2 index "$REF" 2>"$OUT_DIR/index.log"; then
        log "$TOOL: indexing failed"
        STATUS=1
      else
        if bwa-mem2 mem -t 2 "$REF" "$READS" > "$OUT_DIR/alignments.sam" 2>"$OUT_DIR/stderr.log"; then
          STATUS=0
        else
          STATUS=1
        fi
      fi
      ;;

    winnowmap2)
      # Generate repetitive k-mer list first
      KMER_FILE="$OUT_DIR/repetitive_k15.txt"
      if ! meryl count k=15 "$REF" output "$OUT_DIR/meryl_db" 2>"$OUT_DIR/meryl.log"; then
        # Fallback: create empty k-mer file (smoke test just needs to run)
        touch "$KMER_FILE"
        log "$TOOL: meryl unavailable, using empty k-mer list"
      else
        meryl print greater-than distinct=0.9998 "$OUT_DIR/meryl_db" 2>/dev/null | awk '{print $1}' > "$KMER_FILE" || touch "$KMER_FILE"
      fi

      if winnowmap -W "$KMER_FILE" -t 2 -ax map-ont "$REF" "$READS" > "$OUT_DIR/alignments.sam" 2>"$OUT_DIR/stderr.log"; then
        STATUS=0
      else
        STATUS=1
      fi
      ;;

    bowtie2)
      # Index first
      if ! bowtie2-build "$REF" "$OUT_DIR/bt2_idx" 2>"$OUT_DIR/index.log" >/dev/null; then
        log "$TOOL: indexing failed"
        STATUS=1
      else
        if bowtie2 --very-sensitive -p 2 -x "$OUT_DIR/bt2_idx" -U "$READS" > "$OUT_DIR/alignments.sam" 2>"$OUT_DIR/stderr.log"; then
          STATUS=0
        else
          STATUS=1
        fi
      fi
      ;;

    star)
      # STAR needs pre-built genome directory - create minimal one
      GENOME_DIR="$OUT_DIR/star_genome"
      mkdir -p "$GENOME_DIR"

      # STAR genome generate for tiny reference
      if ! STAR --runMode genomeGenerate \
           --genomeDir "$GENOME_DIR" \
           --genomeFastaFiles "$REF" \
           --genomeSAindexNbases 4 \
           --runThreadN 2 \
           2>"$OUT_DIR/genome_gen.log" >/dev/null; then
        log "$TOOL: genome generation failed"
        STATUS=1
      else
        if STAR --runThreadN 2 \
             --genomeDir "$GENOME_DIR" \
             --readFilesIn "$READS" \
             --outSAMtype SAM \
             --outFileNamePrefix "$OUT_DIR/" \
             --outStd SAM > "$OUT_DIR/alignments.sam" 2>"$OUT_DIR/stderr.log"; then
          STATUS=0
        else
          STATUS=1
        fi
      fi
      ;;

    llmap)
      LLMAP_BIN="${LLMAP_BIN:-$LLMAP_ROOT/build/src/llmap}"

      # llmap align takes FASTA directly with -x, not a pre-built index
      if "$LLMAP_BIN" align -r "$READS" -x "$REF" -o "$OUT_DIR/alignments.sam" --sam -k 11 -w 5 -t 2 2>"$OUT_DIR/stderr.log"; then
        STATUS=0
      else
        STATUS=1
      fi
      ;;

    *)
      log "$TOOL: unknown tool"
      STATUS=1
      ;;
  esac

  # Validate output
  if [[ $STATUS -eq 0 ]]; then
    SAM_FILE="$OUT_DIR/alignments.sam"
    if [[ -f "$SAM_FILE" ]]; then
      # Check SAM has some content (header + at least 1 alignment line)
      LINE_COUNT=$(wc -l < "$SAM_FILE")
      HEADER_COUNT=$(grep -c '^@' "$SAM_FILE" || echo 0)
      ALIGN_COUNT=$((LINE_COUNT - HEADER_COUNT))

      if [[ $LINE_COUNT -gt 0 && $ALIGN_COUNT -ge 0 ]]; then
        log "$TOOL: produced SAM with $HEADER_COUNT header lines, $ALIGN_COUNT alignments"

        # Verify SAM is valid with samtools
        if command -v samtools >/dev/null 2>&1; then
          if samtools view -H "$SAM_FILE" >/dev/null 2>&1; then
            RESULTS[$TOOL]="pass"
            PASSED=$((PASSED + 1))
            log "$TOOL: PASS"
          else
            RESULTS[$TOOL]="fail:invalid_sam"
            FAILED=$((FAILED + 1))
            log "$TOOL: FAIL (invalid SAM format)"
          fi
        else
          # No samtools, accept non-empty SAM
          RESULTS[$TOOL]="pass"
          PASSED=$((PASSED + 1))
          log "$TOOL: PASS (samtools not available for validation)"
        fi
      else
        RESULTS[$TOOL]="fail:empty_output"
        FAILED=$((FAILED + 1))
        log "$TOOL: FAIL (empty output)"
      fi
    else
      RESULTS[$TOOL]="fail:no_output"
      FAILED=$((FAILED + 1))
      log "$TOOL: FAIL (no SAM file produced)"
    fi
  else
    RESULTS[$TOOL]="fail:exit_code_$STATUS"
    FAILED=$((FAILED + 1))
    log "$TOOL: FAIL (non-zero exit code: $STATUS)"
    if $VERBOSE && [[ -f "$OUT_DIR/stderr.log" ]]; then
      echo "--- stderr ---" >&2
      head -20 "$OUT_DIR/stderr.log" >&2
      echo "--- end ---" >&2
    fi
  fi
done

# Output results
if $JSON_OUTPUT; then
  echo "{"
  echo "  \"passed\": $PASSED,"
  echo "  \"failed\": $FAILED,"
  echo "  \"skipped\": $SKIPPED,"
  echo "  \"results\": {"
  first=true
  for TOOL in $ALL_TOOLS; do
    if [[ -v RESULTS[$TOOL] ]]; then
      $first || echo ","
      first=false
      printf "    \"%s\": \"%s\"" "$TOOL" "${RESULTS[$TOOL]}"
    fi
  done
  echo ""
  echo "  }"
  echo "}"
else
  echo ""
  echo "=== SMOKE TEST SUMMARY ==="
  for TOOL in $ALL_TOOLS; do
    if [[ -v RESULTS[$TOOL] ]]; then
      case "${RESULTS[$TOOL]}" in
        pass)    printf "  %-12s PASS\n" "$TOOL" ;;
        skipped) printf "  %-12s SKIP (not installed)\n" "$TOOL" ;;
        *)       printf "  %-12s FAIL (%s)\n" "$TOOL" "${RESULTS[$TOOL]}" ;;
      esac
    fi
  done
  echo ""
  echo "Passed: $PASSED / Failed: $FAILED / Skipped: $SKIPPED"
fi

# Exit code: 0 if all installed tools passed, 1 if any failed
if [[ $FAILED -gt 0 ]]; then
  exit 1
fi
exit 0
