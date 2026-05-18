#!/usr/bin/env bash
# the HPC cluster_bench.slurm - SLURM array job for LLmap cross-species bench on the HPC cluster.
#
# One array task = one (organism, tier) cell from the manifest.
# Resources are per-task defaults; tier 5/6/10 tasks override via scontrol
# update because we cannot have per-task --mem in a single array directive.
# Workaround: we request the *maximum* (48G / 8h / 8 cpus) so heavy cells fit;
# light cells just leave most of it unused. SLURM does not bill unused mem.
#
# Submitted from the local box by the HPC cluster_bench_submit.sh.
#
# the HPC cluster specifics (the previous module-load chain did not work):
#   - `module load <foo>` for gcc/cmake/python/samtools/minimap2 is NOT
#     available on the compute nodes. Only `apptainer` is.
#   - Everything llmap needs (gcc-13 / libstdc++ for GLIBCXX_3.4.32,
#     onnxruntime, minimap2, samtools) is bundled in llmap.sif.
#   - Fallback to /home/<user>/miniforge3/bin for samtools/minimap2 if
#     the container ever fails for some reason.

#SBATCH --job-name=llmap_bench
#SBATCH --partition=std
#SBATCH --account=your_slurm_account
#SBATCH --array=1-85
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --chdir=/beegfs/u/<user>/llmap_bench
#SBATCH --output=/beegfs/u/<user>/llmap_bench/slurm_logs/bench_%A_%a.out
#SBATCH --error=/beegfs/u/<user>/llmap_bench/slurm_logs/bench_%A_%a.err
# (no --mem: the HPC cluster disallows; memory is proportional to cpus-per-task)

set -uo pipefail

# ---------- locate workspace ----------
BENCH_ROOT="/beegfs/u/<user>/llmap_bench"
cd "$BENCH_ROOT" || { echo "FATAL: cannot cd $BENCH_ROOT" >&2; exit 1; }

mkdir -p slurm_logs

MANIFEST="$BENCH_ROOT/scripts/the HPC cluster_bench_array_manifest.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo "FATAL: manifest not found: $MANIFEST" >&2
  exit 1
fi

# ---------- pick our row ----------
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
ROW="$(awk -F'\t' -v t="$TASK_ID" 'NR>1 && $1==t {print; exit}' "$MANIFEST")"
if [[ -z "$ROW" ]]; then
  echo "FATAL: no manifest row for task_id=$TASK_ID" >&2
  exit 1
fi
IFS=$'\t' read -r MROW_ID ORGANISM TIER MEM_GB TLIMIT <<< "$ROW"

echo "[the HPC cluster_bench] task=$TASK_ID  organism=$ORGANISM  tier=$TIER  mem=${MEM_GB}G  time=$TLIMIT"
echo "[the HPC cluster_bench] host=$(hostname)  date=$(date -Iseconds)"
echo "[the HPC cluster_bench] job=${SLURM_JOB_ID:-?}  array=${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID:-?}"

# ---------- environment ----------
# the HPC cluster: use the natively-built llmap binary + miniforge3 conda for
# samtools/minimap2. No container, no module system, no apptainer needed.
source /home/<user>/miniforge3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate base 2>/dev/null || true

export PATH="/home/<user>/miniforge3/bin:$PATH"
export TMPDIR="/beegfs/u/<user>/tmp"
mkdir -p "$TMPDIR"

export LLMAP_BIN="$BENCH_ROOT/build_the HPC cluster/src/llmap"
export MINIMAP2_BIN="/home/<user>/miniforge3/bin/minimap2"
export SAMTOOLS_BIN="/home/<user>/miniforge3/bin/samtools"

if [[ ! -x "$LLMAP_BIN" ]]; then
  echo "FATAL: llmap binary not found at $LLMAP_BIN" >&2
  echo "       Rebuild on the HPC cluster: cmake -S . -B build_the HPC cluster -DCMAKE_BUILD_TYPE=Release -DLLMAP_BUILD_TESTS=OFF" >&2
  echo "                          cmake --build build_the HPC cluster -j 4 --target llmap" >&2
  exit 1
fi

echo "[the HPC cluster_bench] llmap:    $LLMAP_BIN ($("$LLMAP_BIN" --version 2>&1 | head -1))"
echo "[the HPC cluster_bench] minimap2: $MINIMAP2_BIN ($("$MINIMAP2_BIN" --version 2>&1 | head -1))"
echo "[the HPC cluster_bench] samtools: $SAMTOOLS_BIN ($("$SAMTOOLS_BIN" --version 2>&1 | head -1))"

# ---------- run the bench ----------
cd "$BENCH_ROOT/benchmarks/scripts" || exit 1

# Refresh the mapper-availability cache on this node.
bash check_mapper_availability.sh --force --quiet || true

THREADS="${SLURM_CPUS_PER_TASK:-4}"

set +e
bash run_species_bench.sh \
    --organism "$ORGANISM" \
    --tier "$TIER" \
    --threads "$THREADS"
RC=$?
set -e

# Mirror outputs into the per-cell reports dir on BeeGFS (run_species_bench
# already writes to $BENCH_ROOT/benchmarks/reports/<organism>/tier<N>/).
CELL_DIR="$BENCH_ROOT/benchmarks/reports/$ORGANISM/tier${TIER}"
mkdir -p "$CELL_DIR"
echo "{\"task_id\": $TASK_ID, \"organism\": \"$ORGANISM\", \"tier\": $TIER, \"exit\": $RC, \"job\": \"${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID:-?}\", \"host\": \"$(hostname)\", \"end\": \"$(date -Iseconds)\"}" \
    > "$CELL_DIR/the HPC cluster_task.json"

if [[ $RC -eq 0 ]]; then
  touch "$CELL_DIR/done.flag"
  echo "[the HPC cluster_bench] OK  $ORGANISM tier=$TIER"
else
  echo "[the HPC cluster_bench] FAIL $ORGANISM tier=$TIER rc=$RC"
fi

exit $RC
