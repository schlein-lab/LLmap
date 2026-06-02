#!/usr/bin/env bash
# the HPC cluster_bench_status.sh - lightweight status snapshot of the the HPC cluster bench.
#
# Shows:
#   1. squeue rows for our user (array-expanded)
#   2. count of done.flag / the HPC cluster_task.json files (per-cell completion)
#   3. tail of the most recent slurm_logs/bench_*.out and .err
#
# Usage:
#   the HPC cluster_bench_status.sh                 # default: 30-line squeue + last 20 log lines
#   the HPC cluster_bench_status.sh --full-queue    # whole squeue, not just head
#   the HPC cluster_bench_status.sh --tail N        # last N lines of latest log (default 20)

set -uo pipefail

REMOTE_HOST="<hpc-login>"
REMOTE_ROOT="/beegfs/u/<user>/llmap_bench"
REMOTE_USER="<user>"

FULL_QUEUE=false
TAIL_N=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --full-queue) FULL_QUEUE=true; shift ;;
    --tail) TAIL_N="$2"; shift 2 ;;
    -h|--help) sed -n '1,15p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
  esac
done

echo "============================================================"
echo " the HPC cluster bench status  ($(date -Iseconds))"
echo "============================================================"
echo

echo "[1] SLURM queue (user=$REMOTE_USER):"
if $FULL_QUEUE; then
  ssh "$REMOTE_HOST" "squeue -u $REMOTE_USER --array -o '%.18i %.9P %.20j %.8T %.10M %.10l %.6D %R'" \
    || echo "  (squeue failed)"
else
  ssh "$REMOTE_HOST" "squeue -u $REMOTE_USER --array -o '%.18i %.9P %.20j %.8T %.10M %.10l %.6D %R' | head -30" \
    || echo "  (squeue failed)"
fi
echo

echo "[2] Per-cell completion (done.flag):"
ssh "$REMOTE_HOST" "
  cd $REMOTE_ROOT/benchmarks/reports 2>/dev/null || { echo '  (no reports dir yet)'; exit 0; }
  total=\$(find . -mindepth 2 -maxdepth 2 -type d | wc -l)
  done=\$(find . -name done.flag -type f | wc -l)
  failed=\$(find . -name the HPC cluster_task.json -type f -exec grep -l '\"exit\": [^0]' {} \; 2>/dev/null | wc -l)
  echo \"  done=\$done  failed=\$failed  cells_seen=\$total  expected=85\"
" || true
echo

echo "[3] Latest SLURM log (tail -n $TAIL_N):"
ssh "$REMOTE_HOST" "
  cd $REMOTE_ROOT/slurm_logs 2>/dev/null || { echo '  (no slurm_logs dir yet)'; exit 0; }
  latest_out=\$(ls -t bench_*.out 2>/dev/null | head -1)
  latest_err=\$(ls -t bench_*.err 2>/dev/null | head -1)
  if [[ -n \"\$latest_out\" ]]; then
    echo \"  --- \$latest_out ---\"
    tail -n $TAIL_N \"\$latest_out\"
  else
    echo '  (no .out files yet)'
  fi
  if [[ -n \"\$latest_err\" && -s \"\$latest_err\" ]]; then
    echo
    echo \"  --- \$latest_err (non-empty) ---\"
    tail -n $TAIL_N \"\$latest_err\"
  fi
" || true
echo
