#!/usr/bin/env bash
# Runner template for the Phase 11 benchmark.
# Each per-tool runner sources this and overrides ALIGN_CMD.
#
# Contract:
#   Inputs (env vars):
#     TOOL          short name of the tool, e.g. "minimap2"
#     TASK          task id (T1..T6)
#     REPLICATE     replicate index (0..2)
#     REF           reference FASTA path
#     READS_R1      reads (single-end) or R1 (paired-end)
#     READS_R2      R2 (paired-end only; empty otherwise)
#     PRESET        tool-specific preset (e.g. "map-hifi", "sr", "splice")
#     THREADS       thread count
#     OUTPUT_DIR    where to write outputs
#     SEED          random seed for this replicate
#
#   Outputs (written to OUTPUT_DIR):
#     alignments.bam        sorted+indexed BAM
#     manifest.json         tool/version/command/dataset SHAs/seed/host
#     resources.json        wallclock/cpu/RAM/output_size from /usr/bin/time -v
#     done.flag             touched on success (idempotency marker)
#
# The runner does NOT compute mapping metrics — that's metrics/compute.py.

set -euo pipefail

: "${TOOL:?TOOL not set}"
: "${TASK:?TASK not set}"
: "${REPLICATE:?REPLICATE not set}"
: "${REF:?REF not set}"
: "${READS_R1:?READS_R1 not set}"
: "${PRESET:?PRESET not set}"
: "${THREADS:=16}"
: "${OUTPUT_DIR:?OUTPUT_DIR not set}"
: "${SEED:=42}"

mkdir -p "$OUTPUT_DIR"
TIME_LOG="$OUTPUT_DIR/time.log"
COMMAND_LOG="$OUTPUT_DIR/command.log"

# Skip if already done
if [[ -f "$OUTPUT_DIR/done.flag" ]]; then
  echo "[$TOOL/$TASK/$REPLICATE] already done, skipping" >&2
  exit 0
fi

# --- override these in the per-tool runner ---
ALIGN_CMD="${ALIGN_CMD:?per-tool runner must set ALIGN_CMD}"
TOOL_VERSION_CMD="${TOOL_VERSION_CMD:?per-tool runner must set TOOL_VERSION_CMD}"

# Resolve and record tool version
TOOL_VERSION="$($TOOL_VERSION_CMD 2>&1 | head -3 | tr '\n' ' ')"

# Hash inputs
sha256() { sha256sum "$1" | awk '{print $1}'; }
REF_SHA="$(sha256 "$REF")"
R1_SHA="$(sha256 "$READS_R1")"
R2_SHA=""
[[ -n "${READS_R2:-}" && -f "${READS_R2:-}" ]] && R2_SHA="$(sha256 "$READS_R2")"

# Record command
echo "$ALIGN_CMD" > "$COMMAND_LOG"

# Run with /usr/bin/time -v to capture resource usage
echo "[$TOOL/$TASK/$REPLICATE] starting at $(date -Iseconds)" >&2
/usr/bin/time -v -o "$TIME_LOG" bash -c "$ALIGN_CMD" \
  > "$OUTPUT_DIR/alignments.sam" \
  2> "$OUTPUT_DIR/stderr.log"

# Sort + index BAM
samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/alignments.bam" "$OUTPUT_DIR/alignments.sam"
samtools index "$OUTPUT_DIR/alignments.bam"
rm -f "$OUTPUT_DIR/alignments.sam"

OUTPUT_SIZE="$(stat -c %s "$OUTPUT_DIR/alignments.bam")"

# Parse /usr/bin/time -v output into resources.json
python3 - "$TIME_LOG" "$OUTPUT_SIZE" > "$OUTPUT_DIR/resources.json" <<'PY'
import json, re, sys
time_log, output_size = sys.argv[1], int(sys.argv[2])
with open(time_log) as f:
    text = f.read()
def grab(pat, cast=float):
    m = re.search(pat, text)
    return cast(m.group(1)) if m else None
wall_str = grab(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', str)
def parse_wall(s):
    parts = s.split(':')
    parts = [float(p) for p in parts]
    if len(parts) == 2:
        return parts[0]*60 + parts[1]
    if len(parts) == 3:
        return parts[0]*3600 + parts[1]*60 + parts[2]
    return None
out = {
    "wallclock_seconds": parse_wall(wall_str) if wall_str else None,
    "user_cpu_seconds": grab(r'User time \(seconds\): ([0-9.]+)'),
    "system_cpu_seconds": grab(r'System time \(seconds\): ([0-9.]+)'),
    "peak_rss_bytes": int((grab(r'Maximum resident set size \(kbytes\): (\d+)') or 0) * 1024),
    "exit_status": int(grab(r'Exit status: (\d+)', int) or -1),
    "output_file_bytes": output_size,
}
json.dump(out, sys.stdout, indent=2)
PY

# Manifest
python3 - > "$OUTPUT_DIR/manifest.json" <<PY
import json, os, socket, time
manifest = {
    "tool": "$TOOL",
    "tool_version": "$TOOL_VERSION".strip(),
    "task": "$TASK",
    "replicate": int("$REPLICATE"),
    "seed": int("$SEED"),
    "threads": int("$THREADS"),
    "preset": "$PRESET",
    "host": socket.gethostname(),
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    "command": open("$COMMAND_LOG").read().strip(),
    "inputs": {
        "ref":      {"path": "$REF",      "sha256": "$REF_SHA"},
        "reads_r1": {"path": "$READS_R1", "sha256": "$R1_SHA"},
        "reads_r2": {"path": "$READS_R2", "sha256": "$R2_SHA"} if "$READS_R2" else None,
    },
}
print(json.dumps(manifest, indent=2))
PY

touch "$OUTPUT_DIR/done.flag"
echo "[$TOOL/$TASK/$REPLICATE] done at $(date -Iseconds)" >&2
