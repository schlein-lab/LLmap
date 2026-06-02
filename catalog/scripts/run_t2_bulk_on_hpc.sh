#!/usr/bin/env bash
# Sync catalog/ to the HPC cluster and run build_t2_bulk.py on the login node.
#
# Per project policy 2026-06-02: NO local heavy compute. Parsing UCSC +
# T2T BEDs is light enough for <hpc-login> (< 60 min wall), so we do not
# need SLURM for this. If wall time grows past 60 min (extra assemblies,
# gnomAD-SV overlap, etc.) — switch to an sbatch wrapper.
#
# Usage:
#   bash catalog/scripts/run_t2_bulk_on_the HPC cluster.sh                 # default: grch38+chm13, v2026.Q2
#   bash catalog/scripts/run_t2_bulk_on_the HPC cluster.sh --dry-run
#   ASSEMBLIES="grch38 chm13 mm39" bash catalog/scripts/run_t2_bulk_on_the HPC cluster.sh

set -uo pipefail

REMOTE_HOST="<hpc-login>"
REMOTE_ROOT="/beegfs/u/<user>/llmap_catalog_build"

RELEASE="${RELEASE:-v2026.Q2}"
ASSEMBLIES="${ASSEMBLIES:-grch38 chm13}"

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CATALOG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

say()  { echo "[run_t2_bulk] $*" >&2; }
run()  { echo "  \$ $*"; $DRY_RUN || eval "$@"; }

say "release    = $RELEASE"
say "assemblies = $ASSEMBLIES"
say "remote     = $REMOTE_HOST:$REMOTE_ROOT"

# ---- preflight ssh
if ! $DRY_RUN && ! ssh -o BatchMode=yes -o ConnectTimeout=8 "$REMOTE_HOST" true 2>/dev/null; then
  echo "ERROR: cannot ssh $REMOTE_HOST (VPN down?)" >&2
  exit 1
fi

# ---- sync catalog/ (script + schema; no large artifacts)
run "ssh '$REMOTE_HOST' 'mkdir -p $REMOTE_ROOT/catalog/bulk'"
run "rsync -az --human-readable --delete \
    --exclude 'bulk/' --exclude '__pycache__' \
    '$CATALOG_ROOT/' '$REMOTE_HOST:$REMOTE_ROOT/catalog/'"

# ---- run the build on <hpc-login>
# python3 is in miniforge3 base on the HPC cluster (per the HPC cluster_build_slurm_ops memory)
REMOTE_CMD="source /home/<user>/miniforge3/etc/profile.d/conda.sh && conda activate base && \
    cd $REMOTE_ROOT && \
    python3 catalog/scripts/build_t2_bulk.py \
        --release $RELEASE \
        --assemblies $ASSEMBLIES \
        --out catalog/bulk/"

say "running build on $REMOTE_HOST"
run "ssh '$REMOTE_HOST' '$REMOTE_CMD'"

# ---- pull summaries back locally so we can inspect without re-ssh
say "pulling summaries back"
run "rsync -az '$REMOTE_HOST:$REMOTE_ROOT/catalog/bulk/*.summary.json' \
              '$REMOTE_HOST:$REMOTE_ROOT/catalog/bulk/*.release_summary.json' \
              '$CATALOG_ROOT/bulk/' 2>/dev/null || true"

say "done — JSONL stays on the HPC cluster ($REMOTE_ROOT/catalog/bulk/), summaries pulled locally"
