#!/usr/bin/env bash
# the HPC cluster_bench_submit.sh - sync LLmap bench to the HPC cluster and submit SLURM array.
#
# Usage:
#   the HPC cluster_bench_submit.sh                 # rsync + sbatch
#   the HPC cluster_bench_submit.sh --dry-run       # print every step, run nothing
#   the HPC cluster_bench_submit.sh --no-submit     # rsync only, skip sbatch
#   the HPC cluster_bench_submit.sh --test-only     # rsync + sbatch --test-only
#
# Companion scripts:
#   the HPC cluster_bench_status.sh   - squeue + tail latest log
#   the HPC cluster_bench_collect.sh  - rsync results back, rebuild scoreboard + site

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

REMOTE_HOST="<hpc-login>"
REMOTE_USER="<user>"
REMOTE_ROOT="/beegfs/u/<user>/llmap_bench"

# the HPC cluster: login node lacks munge; sbatch must run on the SLURM frontend.
# rsync goes via login (BeeGFS is shared), sbatch via front1 (ProxyJump'd).
REMOTE_SUBMIT_HOST="<hpc-frontend>"
REMOTE_SBATCH="/syssw/slurm/current/bin/sbatch"

DRY_RUN=false
NO_SUBMIT=false
TEST_ONLY=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)   DRY_RUN=true; shift ;;
    --no-submit) NO_SUBMIT=true; shift ;;
    --test-only) TEST_ONLY=true; shift ;;
    -h|--help)
      sed -n '1,15p' "$0"
      exit 0 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
  esac
done

say() { echo "[the HPC cluster_bench_submit] $*" >&2; }

run() {
  # echo + execute (or just echo in dry-run)
  echo "  \$ $*"
  if ! $DRY_RUN; then
    eval "$@"
  fi
}

# ---------- preflight ----------
if ! $DRY_RUN && ! ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE_HOST" true 2>/dev/null; then
  echo "ERROR: cannot ssh $REMOTE_HOST (no batchmode key/agent). Try: ssh $REMOTE_HOST" >&2
  exit 1
fi

say "local root  = $LLMAP_ROOT"
say "remote host = $REMOTE_HOST"
say "remote root = $REMOTE_ROOT"
say "dry-run     = $DRY_RUN"

# ---------- step 1: ensure remote dirs ----------
say "step 1: ensure remote workspace exists"
run "ssh '$REMOTE_HOST' 'mkdir -p $REMOTE_ROOT/slurm_logs $REMOTE_ROOT/benchmarks/reports'"

# ---------- step 2: rsync each tree we need ----------
RSYNC_OPTS=(-az --delete --human-readable
            --exclude '__pycache__' --exclude '*.pyc' --exclude '.venv'
            --exclude '.git' --exclude 'build/CMakeFiles'
            --exclude 'benchmarks/reports'    # don't clobber remote results
            --exclude 'benchmarks/datasets/cache' # huge, generated on remote
            --exclude '.locks')

if $DRY_RUN; then
  RSYNC_OPTS+=(--dry-run -i)
fi

for sub in benchmarks knowledge build src scripts CMakeLists.txt; do
  if [[ -e "$LLMAP_ROOT/$sub" ]]; then
    say "step 2: rsync $sub/ -> $REMOTE_HOST:$REMOTE_ROOT/$sub"
    run "rsync ${RSYNC_OPTS[*]} '$LLMAP_ROOT/$sub' '$REMOTE_HOST:$REMOTE_ROOT/'"
  fi
done

# Explicitly ensure the SLURM script + manifest land where the job expects them.
say "step 2b: rsync scripts/the HPC cluster_bench.* (explicit, in case scripts/ was skipped)"
run "rsync -avz '$LLMAP_ROOT/scripts/the HPC cluster_bench.slurm' '$LLMAP_ROOT/scripts/the HPC cluster_bench_array_manifest.tsv' '$REMOTE_HOST:$REMOTE_ROOT/scripts/'"

# ---------- step 3: submit ----------
if $NO_SUBMIT; then
  say "step 3: --no-submit set, stopping after rsync"
  exit 0
fi

if $TEST_ONLY; then
  say "step 3: sbatch --test-only (validates SLURM syntax, does not run)"
  run "ssh '$REMOTE_SUBMIT_HOST' 'cd $REMOTE_ROOT && $REMOTE_SBATCH --test-only scripts/the HPC cluster_bench.slurm'"
  exit 0
fi

if $DRY_RUN; then
  say "step 3: (dry-run) would submit:"
  echo "  \$ ssh $REMOTE_SUBMIT_HOST 'cd $REMOTE_ROOT && $REMOTE_SBATCH scripts/the HPC cluster_bench.slurm'"
  exit 0
fi

say "step 3: submitting SLURM array job"
SUBMIT_OUT="$(ssh "$REMOTE_SUBMIT_HOST" "cd $REMOTE_ROOT && $REMOTE_SBATCH scripts/the HPC cluster_bench.slurm" 2>&1)"
echo "$SUBMIT_OUT"
JOB_ID="$(echo "$SUBMIT_OUT" | awk '/Submitted batch job/ {print $4; exit}')"

if [[ -z "$JOB_ID" ]]; then
  echo "ERROR: could not parse job id from sbatch output" >&2
  exit 1
fi

say "submitted job id: $JOB_ID  (array 1-85)"
cat <<EOF >&2

monitor with:
  bash $SCRIPT_DIR/the HPC cluster_bench_status.sh
  ssh $REMOTE_HOST 'squeue -j $JOB_ID --array | head -30'
  ssh $REMOTE_HOST 'tail -f $REMOTE_ROOT/slurm_logs/bench_${JOB_ID}_1.out'

collect results when done:
  bash $SCRIPT_DIR/the HPC cluster_bench_collect.sh
EOF
