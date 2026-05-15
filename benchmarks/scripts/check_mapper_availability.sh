#!/usr/bin/env bash
# check_mapper_availability.sh
#
# Probe the host for the 4 cross-species benchmark mappers and cache the
# result to benchmarks/.mapper_status.tsv. Idempotent: re-running just
# refreshes the cache.
#
# Each mapper is checked in this order:
#   1. native binary on $PATH
#   2. apptainer container at benchmarks/containers/<mapper>.sif (if apptainer present)
#   3. otherwise UNAVAILABLE
#
# Output columns (tab-separated):
#   mapper  status  source  version  command_prefix
#
# Where:
#   status         = OK | UNAVAILABLE
#   source         = native | apptainer | -
#   command_prefix = "" for native, "apptainer exec ... " for apptainer
#
# Other scripts source this file via SCOREBOARD-style lookup or just read
# the TSV.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
LLMAP_ROOT="$(dirname "$BENCH_ROOT")"
STATUS_FILE="$BENCH_ROOT/.mapper_status.tsv"
CONTAINER_DIR="$LLMAP_ROOT/containers"

LLMAP_BIN_DEFAULT="$LLMAP_ROOT/build/src/llmap"
LLMAP_BIN="${LLMAP_BIN:-$LLMAP_BIN_DEFAULT}"

QUIET=false
FORCE=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --quiet) QUIET=true; shift ;;
    --force) FORCE=true; shift ;;
    -h|--help)
      sed -n '1,30p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
  esac
done

log() { $QUIET || echo "[mapper-check] $*" >&2; }

have_apptainer=false
if command -v apptainer >/dev/null 2>&1; then
  have_apptainer=true
fi

# probe_one <mapper> <native_check_cmd> <version_cmd_native>
# Returns four fields via stdout: status<TAB>source<TAB>version<TAB>prefix
probe_one() {
  local mapper="$1" native_bin="$2" version_cmd="$3"
  local sif="$CONTAINER_DIR/${mapper}.sif"
  if [[ "$mapper" == "llmap" ]]; then
    # llmap is special: explicit binary path
    if [[ -x "$LLMAP_BIN" ]]; then
      local v
      v="$("$LLMAP_BIN" --version 2>&1 | head -1 | tr '\t' ' ' || echo unknown)"
      printf 'OK\tnative\t%s\t%s\n' "$v" "$LLMAP_BIN"
      return
    fi
  else
    if command -v "$native_bin" >/dev/null 2>&1; then
      local v
      v="$(eval "$version_cmd" 2>&1 | head -1 | tr '\t' ' ' || echo unknown)"
      printf 'OK\tnative\t%s\t%s\n' "$v" "$native_bin"
      return
    fi
  fi
  if $have_apptainer && [[ -f "$sif" ]]; then
    local v
    v="$(apptainer exec "$sif" "$native_bin" --version 2>&1 | head -1 | tr '\t' ' ' || echo unknown)"
    printf 'OK\tapptainer\t%s\tapptainer exec %s %s\n' "$v" "$sif" "$native_bin"
    return
  fi
  printf 'UNAVAILABLE\t-\t-\t-\n'
}

# Cache short-circuit
if [[ -f "$STATUS_FILE" && "$FORCE" != "true" ]]; then
  log "cache hit: $STATUS_FILE (use --force to refresh)"
  cat "$STATUS_FILE"
  exit 0
fi

log "probing mappers on $(hostname)"

tmp="$(mktemp)"
{
  printf 'mapper\tstatus\tsource\tversion\tcommand_prefix\n'
  while IFS='|' read -r mapper bin vcmd; do
    [[ -z "$mapper" ]] && continue
    fields="$(probe_one "$mapper" "$bin" "$vcmd")"
    printf '%s\t%s\n' "$mapper" "$fields"
  done <<EOF
llmap|llmap|llmap --version
minimap2|minimap2|minimap2 --version
bwa-mem2|bwa-mem2|bwa-mem2 version
winnowmap|winnowmap|winnowmap --version
EOF
} > "$tmp"

mv "$tmp" "$STATUS_FILE"

if ! $QUIET; then
  log "wrote $STATUS_FILE"
  column -t -s $'\t' "$STATUS_FILE" >&2 || cat "$STATUS_FILE" >&2
fi

cat "$STATUS_FILE"
