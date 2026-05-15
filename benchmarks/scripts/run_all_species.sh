#!/usr/bin/env bash
# run_all_species.sh — loop over all registered organisms x informative tiers.
#
# Tiers swept: 1, 2, 5, 6, 10 (the most-informative subset).
# Organisms come from knowledge/organisms/ (one dir per organism), minus
# meta-buckets like 'test_fictive'. Use --organisms to override.
#
# Usage:
#   run_all_species.sh [--organisms human,mouse,...] [--tiers 1,2,5,6,10]
#                      [--mappers llmap,minimap2,bwa-mem2,winnowmap]
#                      [--threads N] [--continue-on-error]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
LLMAP_ROOT="$(dirname "$BENCH_ROOT")"
ORG_DIR="$LLMAP_ROOT/knowledge/organisms"
SCOREBOARD="$BENCH_ROOT/scoreboard.tsv"

DEFAULT_TIERS="1,2,5,6,10"
DEFAULT_MAPPERS="llmap,minimap2,bwa-mem2,winnowmap"
THREADS=4
CONTINUE=true   # default: skip failing tiers but keep going

ORGANISMS=""
TIERS="$DEFAULT_TIERS"
MAPPERS="$DEFAULT_MAPPERS"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --organisms) ORGANISMS="$2"; shift 2 ;;
    --tiers)     TIERS="$2"; shift 2 ;;
    --mappers)   MAPPERS="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --abort-on-error) CONTINUE=false; shift ;;
    --continue-on-error) CONTINUE=true; shift ;;
    -h|--help) sed -n '1,20p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown arg $1" >&2; exit 2 ;;
  esac
done

# Run mapper availability probe once up front
bash "$SCRIPT_DIR/check_mapper_availability.sh" --force | column -t -s $'\t' >&2 || true

if [[ -z "$ORGANISMS" ]]; then
  # Auto-discover organisms; drop test/meta buckets
  mapfile -t ORG_ARR < <(find "$ORG_DIR" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' \
    | grep -v -E '^(test_fictive)$' | sort)
  ORGANISMS="$(IFS=','; echo "${ORG_ARR[*]}")"
fi

IFS=',' read -ra ORG_LIST <<< "$ORGANISMS"
IFS=',' read -ra TIER_LIST <<< "$TIERS"

echo "[run_all_species] organisms = ${ORG_LIST[*]}"  >&2
echo "[run_all_species] tiers     = ${TIER_LIST[*]}" >&2
echo "[run_all_species] mappers   = $MAPPERS"        >&2
echo "[run_all_species] scoreboard= $SCOREBOARD"     >&2

ok=0; fail=0
for org in "${ORG_LIST[@]}"; do
  for tier in "${TIER_LIST[@]}"; do
    echo "------------------------------------------------------------" >&2
    echo "[run_all_species] $org tier=$tier" >&2
    if bash "$SCRIPT_DIR/run_species_bench.sh" \
         --organism "$org" --tier "$tier" \
         --mappers "$MAPPERS" --threads "$THREADS"; then
      ok=$((ok+1))
    else
      fail=$((fail+1))
      echo "[run_all_species] FAILED $org tier=$tier" >&2
      $CONTINUE || exit 1
    fi
  done
done

echo "------------------------------------------------------------" >&2
echo "[run_all_species] DONE  ok=$ok fail=$fail  -> $SCOREBOARD" >&2
