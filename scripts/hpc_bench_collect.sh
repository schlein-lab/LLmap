#!/usr/bin/env bash
# the HPC cluster_bench_collect.sh - pull the HPC cluster bench results back, rebuild scoreboard + site.
#
# Steps:
#   1. rsync /beegfs/u/<user>/llmap_bench/benchmarks/reports/   -> local
#      rsync /beegfs/u/<user>/llmap_bench/benchmarks/scoreboard.tsv (if exists)
#   2. run benchmarks/scripts/aggregate_scoreboard.py to regenerate the local
#      scoreboard from the freshly-synced per-cell summary.json files
#   3. run web/losslessmap.com/build.py to rebuild the static site
#
# Usage:
#   the HPC cluster_bench_collect.sh                # full sync + scoreboard + site rebuild
#   the HPC cluster_bench_collect.sh --no-site      # skip site rebuild
#   the HPC cluster_bench_collect.sh --dry-run      # rsync --dry-run, skip downstream

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

REMOTE_HOST="<hpc-login>"
REMOTE_REPORTS="/beegfs/u/<user>/llmap_bench/benchmarks/reports/"
REMOTE_SCOREBOARD="/beegfs/u/<user>/llmap_bench/benchmarks/scoreboard.tsv"
LOCAL_REPORTS="$LLMAP_ROOT/benchmarks/reports/"
LOCAL_SCOREBOARD="$LLMAP_ROOT/benchmarks/scoreboard.tsv"

DRY_RUN=false
NO_SITE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --no-site) NO_SITE=true; shift ;;
    -h|--help) sed -n '1,15p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
  esac
done

say() { echo "[the HPC cluster_bench_collect] $*" >&2; }

RSYNC_FLAGS=(-avz --human-readable)
$DRY_RUN && RSYNC_FLAGS+=(--dry-run -i)

say "step 1a: rsync reports/ from $REMOTE_HOST"
mkdir -p "$LOCAL_REPORTS"
rsync "${RSYNC_FLAGS[@]}" "$REMOTE_HOST:$REMOTE_REPORTS" "$LOCAL_REPORTS"

say "step 1b: rsync remote scoreboard.tsv (if present)"
rsync "${RSYNC_FLAGS[@]}" "$REMOTE_HOST:$REMOTE_SCOREBOARD" "$LOCAL_SCOREBOARD" 2>/dev/null \
  || say "  (no remote scoreboard.tsv yet, will regenerate locally)"

if $DRY_RUN; then
  say "dry-run mode: skipping scoreboard / site rebuild"
  exit 0
fi

# ---------- step 2: regenerate scoreboard locally from synced summaries ----------
AGG="$LLMAP_ROOT/benchmarks/scripts/aggregate_scoreboard.py"
if [[ -x "$AGG" ]]; then
  say "step 2: regenerate scoreboard via $AGG"
  python3 "$AGG" || say "WARN: aggregate_scoreboard.py exited non-zero"
else
  say "WARN: $AGG missing or not executable; skipping aggregation"
fi

# ---------- step 3: rebuild static site ----------
if $NO_SITE; then
  say "step 3: --no-site set, skipping site rebuild"
  exit 0
fi

SITE_BUILD="$LLMAP_ROOT/web/losslessmap.com/build.py"
if [[ -f "$SITE_BUILD" ]]; then
  say "step 3: rebuild losslessmap.com via $SITE_BUILD"
  ( cd "$LLMAP_ROOT/web/losslessmap.com" && python3 build.py ) \
    || say "WARN: site rebuild failed"
else
  say "WARN: $SITE_BUILD missing; skipping site rebuild"
fi

say "done. results in $LOCAL_REPORTS, scoreboard at $LOCAL_SCOREBOARD"
