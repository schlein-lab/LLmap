#!/usr/bin/env bash
# build_llmap_container_the HPC cluster.sh — produce llmap.sif and ship it to the HPC cluster.
#
# the HPC cluster's apptainer module explicitly states:
#   "you cannot directly create containers here due to restrictive permissions
#    in the multiuser environment. ... your local Linux box."
# So we ALWAYS build locally (apptainer or singularity-in-docker) and then
# rsync the .sif to /beegfs/u/<user>/llmap_bench/containers/.
#
# Usage:
#   bash scripts/build_llmap_container_the HPC cluster.sh                 # build + ship
#   bash scripts/build_llmap_container_the HPC cluster.sh --ship-only     # skip build, just rsync existing .sif
#   bash scripts/build_llmap_container_the HPC cluster.sh --no-validate   # don't run the remote --version check
#
# Companion: scripts/build_llmap_container.sh actually does the local build.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLMAP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOCAL_SIF="$LLMAP_ROOT/containers/llmap.sif"

REMOTE_RSYNC_HOST="<hpc-login>"
REMOTE_SUBMIT_HOST="<hpc-frontend>"
REMOTE_USER="<user>"
REMOTE_BENCH="/beegfs/u/${REMOTE_USER}/llmap_bench"
REMOTE_CONT_DIR="${REMOTE_BENCH}/containers"
REMOTE_SIF="${REMOTE_CONT_DIR}/llmap.sif"

# Apptainer is installed via the RRZ module system; the env var pulls in
# the binary at /sw/env/system-gcc/apptainer/1.4.5/bin/apptainer. We probe
# that path directly because `module load` does not work in non-interactive
# ssh on front1.
REMOTE_APPTAINER="/sw/env/system-gcc/apptainer/1.4.5/bin/apptainer"

SHIP_ONLY=false
DO_VALIDATE=true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ship-only)   SHIP_ONLY=true; shift ;;
    --no-validate) DO_VALIDATE=false; shift ;;
    -h|--help)     sed -n '1,20p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
  esac
done

say() { echo "[build_llmap_container_the HPC cluster] $*" >&2; }

# ---- step 1: ensure a local .sif exists ------------------------------------
if ! $SHIP_ONLY; then
  say "step 1: build .sif locally"
  bash "$SCRIPT_DIR/build_llmap_container.sh"
fi

if [[ ! -f "$LOCAL_SIF" ]]; then
  echo "ERROR: $LOCAL_SIF not found. Build first (omit --ship-only)." >&2
  exit 1
fi

say "local .sif: $(ls -lh "$LOCAL_SIF" | awk '{print $5, $9}')"

# ---- step 2: preflight remote ssh ------------------------------------------
if ! ssh -o BatchMode=yes -o ConnectTimeout=10 "$REMOTE_RSYNC_HOST" true 2>/dev/null; then
  echo "ERROR: cannot ssh $REMOTE_RSYNC_HOST" >&2
  exit 1
fi

ssh "$REMOTE_RSYNC_HOST" "mkdir -p '$REMOTE_CONT_DIR'"

# ---- step 3: rsync .sif to BeeGFS ------------------------------------------
say "step 3: rsync .sif -> $REMOTE_RSYNC_HOST:$REMOTE_SIF"
rsync -avz --partial --progress "$LOCAL_SIF" "$REMOTE_RSYNC_HOST:$REMOTE_SIF"

# ---- step 4: validate on the HPC cluster --------------------------------------------
if $DO_VALIDATE; then
  say "step 4: validate via $REMOTE_SUBMIT_HOST (using $REMOTE_APPTAINER)"
  ssh "$REMOTE_SUBMIT_HOST" "test -x $REMOTE_APPTAINER && \
    $REMOTE_APPTAINER exec --bind /beegfs $REMOTE_SIF llmap --version" || {
    say "WARN: front1 validation failed (front1 may lack apptainer setuid)."
    say "      The .sif is on BeeGFS; SLURM compute nodes should be able to run it."
  }
fi

say "DONE. Remote .sif: $REMOTE_SIF"
say "Next: sbatch scripts/the HPC cluster_bench.slurm (after rsyncing scripts/ too)"
