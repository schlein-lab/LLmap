#!/usr/bin/env bash
# run_species_bench.sh — 4-mapper cross-species benchmark runner.
#
# Usage:
#   run_species_bench.sh --organism <name> --tier <N> \
#       [--mappers llmap,minimap2,bwa-mem2,winnowmap] \
#       [--threads N] [--seed N] [--force]
#
# Pipeline:
#   1. gen_synth_reads.py --organism X --tier N -> reads.fastq + truth.tsv
#   2. For each requested mapper that is available on this host:
#        - build index (cached)
#        - map reads -> SAM -> sorted+indexed BAM
#        - write timing.json
#   3. analyze_bench.py -> accuracy.json per mapper + comparison.tsv + summary.json
#   4. Append one summary line per (mapper x tier) to benchmarks/scoreboard.tsv
#
# Mapper availability is probed via check_mapper_availability.sh; missing
# mappers are skipped with a WARN line, never abort the harness.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(dirname "$SCRIPT_DIR")"
LLMAP_ROOT="$(dirname "$BENCH_ROOT")"
REPORTS_ROOT="$BENCH_ROOT/reports"
SCOREBOARD="$BENCH_ROOT/scoreboard.tsv"
STATUS_FILE="$BENCH_ROOT/.mapper_status.tsv"

LLMAP_BIN="${LLMAP_BIN:-$LLMAP_ROOT/build/src/llmap}"

# ---- args --------------------------------------------------------------------
ORGANISM=""
TIER=""
MAPPERS="llmap,minimap2,bwa-mem2,winnowmap"
THREADS=4
SEED=42
FORCE=false

usage() {
  sed -n '1,30p' "$0"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --organism)  ORGANISM="$2"; shift 2 ;;
    --tier)      TIER="$2"; shift 2 ;;
    --mappers)   MAPPERS="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --seed)      SEED="$2"; shift 2 ;;
    --force)     FORCE=true; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "ERROR: unknown arg $1" >&2; usage; exit 2 ;;
  esac
done

[[ -z "$ORGANISM" ]] && { echo "ERROR: --organism required" >&2; exit 2; }
[[ -z "$TIER" ]] && { echo "ERROR: --tier required" >&2; exit 2; }

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not on PATH" >&2
  exit 2
fi

# --- concurrency lock ---------------------------------------------------------
# Prevent two concurrent invocations from writing the same per-tier outputs.
# Without this, parallel runs (e.g. two `run_all_species.sh` started by mistake)
# clobber each other's aln.sam and the OOM killer takes one out mid-write,
# leaving 0-byte aln.sam + 0-byte .time files. Empirically this was the
# root cause of "WARN: llmap produced empty SAM" on tier 5/6.
LOCK_DIR="${BENCH_ROOT}/.locks"
mkdir -p "$LOCK_DIR"
LOCK_FILE="$LOCK_DIR/${ORGANISM}_tier${TIER}.lock"
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  echo "ERROR: another run_species_bench.sh is already running for $ORGANISM tier=$TIER (lock $LOCK_FILE)" >&2
  echo "       refusing to clobber its outputs; aborting." >&2
  exit 4
fi
trap 'flock -u 9 2>/dev/null || true; rm -f "$LOCK_FILE" 2>/dev/null || true' EXIT

# ---- mapper availability ----------------------------------------------------
bash "$SCRIPT_DIR/check_mapper_availability.sh" --quiet >/dev/null

# Returns "OK<TAB>source<TAB>version<TAB>cmdprefix" if mapper is OK, else "".
lookup_mapper() {
  local mapper="$1"
  awk -F'\t' -v m="$mapper" 'NR>1 && $1==m && $2=="OK" {print $2"\t"$3"\t"$4"\t"$5; exit}' "$STATUS_FILE"
}

# ---- output layout ----------------------------------------------------------
TIER_DIR="$REPORTS_ROOT/$ORGANISM/tier${TIER}"
mkdir -p "$TIER_DIR"

# ---- step 1: generate reads -------------------------------------------------
# The generator owns reference creation (it lives under benchmarks/datasets/cache/...).
# It writes reads.fastq OR reads.R1.fastq+reads.R2.fastq, truth.tsv, manifest.json.
NEED_GEN=true
if [[ -f "$TIER_DIR/manifest.json" && "$FORCE" != "true" ]]; then
  if [[ -f "$TIER_DIR/truth.tsv" ]] && { [[ -f "$TIER_DIR/reads.fastq" ]] || [[ -f "$TIER_DIR/reads.R1.fastq" ]]; }; then
    NEED_GEN=false
    echo "[run_species_bench] cached generator output in $TIER_DIR" >&2
  fi
fi
if $NEED_GEN; then
  LLMAP_BIN="$LLMAP_BIN" python3 "$SCRIPT_DIR/gen_synth_reads.py" \
    --organism "$ORGANISM" --tier "$TIER" \
    --output-dir "$TIER_DIR" --seed "$SEED"
fi

TRUTH="$TIER_DIR/truth.tsv"
[[ -f "$TRUTH" ]] || { echo "ERROR: missing $TRUTH after generator" >&2; exit 3; }

# Reference path lives in manifest.json (generator-owned, shared across tiers)
REF="$(python3 -c "import json,sys;print(json.load(open(sys.argv[1])).get('reference_fa',''))" "$TIER_DIR/manifest.json")"
[[ -n "$REF" && -f "$REF" ]] || { echo "ERROR: cannot resolve reference_fa from $TIER_DIR/manifest.json" >&2; exit 3; }

# Detect single-end vs paired reads
if [[ -f "$TIER_DIR/reads.fastq" ]]; then
  READS="$TIER_DIR/reads.fastq"
  READS_R2=""
  PAIRED=false
elif [[ -f "$TIER_DIR/reads.R1.fastq" && -f "$TIER_DIR/reads.R2.fastq" ]]; then
  READS="$TIER_DIR/reads.R1.fastq"
  READS_R2="$TIER_DIR/reads.R2.fastq"
  PAIRED=true
else
  echo "ERROR: no reads.fastq or reads.R[12].fastq in $TIER_DIR" >&2
  exit 3
fi
echo "[run_species_bench] ref=$REF reads=$READS paired=$PAIRED" >&2

# Reference helpers (build once per ref)
ensure_fai() {
  [[ -f "${REF}.fai" ]] || samtools faidx "$REF"
}

# Detect read type for preset selection
detect_long_reads() {
  awk 'NR%4==2 {sum+=length($0); n++} n>=200 {exit} END {if (n>0) printf "%d", sum/n; else printf "0"}' "$READS"
}
MEAN_LEN="$(detect_long_reads)"
IS_LONG_READ=true
if $PAIRED || [[ "$MEAN_LEN" -lt 500 ]]; then
  IS_LONG_READ=false
fi
echo "[run_species_bench] mean read length ~${MEAN_LEN}bp -> long_read=$IS_LONG_READ" >&2

# ---- step 2: per-mapper runs -------------------------------------------------

run_llmap() {
  local out="$TIER_DIR/llmap/rep0"
  mkdir -p "$out"
  local bam="$out/alignments.bam"
  local timing="$out/timing.json"
  if [[ -f "$bam" && "$FORCE" != "true" ]]; then
    echo "[llmap] cached $bam"; return 0
  fi
  local idx="$TIER_DIR/.idx/llmap.idx"
  mkdir -p "$(dirname "$idx")"
  if [[ ! -f "$idx" ]]; then
    echo "[llmap] building index" >&2
    "$LLMAP_BIN" index -r "$REF" -o "$idx" -k 19 -w 19 1>&2 || \
      { echo "WARN: llmap index failed, skipping" >&2; return 0; }
  fi
  local sam="$out/aln.sam"
  local preset="map-hifi"
  $IS_LONG_READ || preset="sr"
  echo "[llmap] aligning preset=$preset -> $bam" >&2
  local t0; t0=$(date +%s.%N)
  local llmap_rc=0
  /usr/bin/time -v -o "$out/.time" \
     "$LLMAP_BIN" align -x "$preset" --reference "$REF" -r "$READS" \
                        -o "$sam" --threads "$THREADS" --sam 2>"$out/stderr.log" \
     || llmap_rc=$?
  local t1; t1=$(date +%s.%N)
  if [[ $llmap_rc -ne 0 ]]; then
    echo "WARN: llmap align failed rc=$llmap_rc (see $out/stderr.log)" >&2
    # Persist a structured marker so analyze_bench/aggregator can see this.
    printf '{"mapper":"llmap","failure_mode":"align_rc_%s","tier_dir":"%s"}\n' \
      "$llmap_rc" "$TIER_DIR" > "$out/failure.json"
    return 0
  fi
  if [[ ! -s "$sam" ]]; then
    echo "WARN: llmap produced empty SAM (preset=$preset reads=$READS ref=$REF)" >&2
    printf '{"mapper":"llmap","failure_mode":"empty_sam","preset":"%s","reads":"%s","ref":"%s"}\n' \
      "$preset" "$READS" "$REF" > "$out/failure.json"
    return 0
  fi
  awk -F'\t' 'NF >= 11 || /^@/' "$sam" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  rm -f "$sam"
  emit_timing "llmap" "$timing" "$t0" "$t1" "$out/.time" "$bam"
}

run_minimap2() {
  local out="$TIER_DIR/minimap2/rep0"
  mkdir -p "$out"
  local bam="$out/alignments.bam"
  local timing="$out/timing.json"
  if [[ -f "$bam" && "$FORCE" != "true" ]]; then
    echo "[minimap2] cached $bam"; return 0
  fi
  local info; info="$(lookup_mapper minimap2)"
  [[ -z "$info" ]] && { echo "WARN: minimap2 not installed, skipping" >&2; return 0; }
  local cmd_prefix; cmd_prefix="$(echo "$info" | awk -F'\t' '$1=="OK" && $2=="apptainer" {print $4}')"
  local preset="map-hifi"
  $IS_LONG_READ || preset="sr"
  echo "[minimap2] preset=$preset -> $bam" >&2
  local sam="$out/aln.sam"
  local extra=""
  $PAIRED && extra="'$READS_R2'"
  local mbin="minimap2"
  [[ -n "$cmd_prefix" ]] && mbin="$cmd_prefix"
  local t0; t0=$(date +%s.%N)
  /usr/bin/time -v -o "$out/.time" bash -c \
    "$mbin -t $THREADS -ax $preset '$REF' '$READS' $extra" > "$sam" 2>>"$out/stderr.log"
  local t1; t1=$(date +%s.%N)
  awk -F'\t' 'NF >= 11 || /^@/' "$sam" | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  rm -f "$sam"
  emit_timing "minimap2" "$timing" "$t0" "$t1" "$out/.time" "$bam"
}

run_bwa_mem2() {
  local out="$TIER_DIR/bwa-mem2/rep0"
  mkdir -p "$out"
  local bam="$out/alignments.bam"
  local timing="$out/timing.json"
  if [[ -f "$bam" && "$FORCE" != "true" ]]; then
    echo "[bwa-mem2] cached $bam"; return 0
  fi
  if $IS_LONG_READ; then
    echo "WARN: bwa-mem2 is short-read only; tier reads look long (~${MEAN_LEN}bp) -- skipping" >&2
    return 0
  fi
  local info; info="$(lookup_mapper bwa-mem2)"
  [[ -z "$info" ]] && { echo "WARN: bwa-mem2 not installed, skipping" >&2; return 0; }
  if [[ ! -f "${REF}.bwt.2bit.64" ]]; then
    echo "[bwa-mem2] building index" >&2
    bwa-mem2 index "$REF" >&2
  fi
  local sam="$out/aln.sam"
  local extra=""
  $PAIRED && extra="'$READS_R2'"
  echo "[bwa-mem2] aligning -> $bam" >&2
  local t0; t0=$(date +%s.%N)
  /usr/bin/time -v -o "$out/.time" bash -c \
    "bwa-mem2 mem -t $THREADS '$REF' '$READS' $extra" > "$sam" 2>>"$out/stderr.log"
  local t1; t1=$(date +%s.%N)
  awk -F'\t' 'NF >= 11 || /^@/' "$sam" | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  rm -f "$sam"
  emit_timing "bwa-mem2" "$timing" "$t0" "$t1" "$out/.time" "$bam"
}

run_winnowmap() {
  local out="$TIER_DIR/winnowmap/rep0"
  mkdir -p "$out"
  local bam="$out/alignments.bam"
  local timing="$out/timing.json"
  if [[ -f "$bam" && "$FORCE" != "true" ]]; then
    echo "[winnowmap] cached $bam"; return 0
  fi
  if ! $IS_LONG_READ; then
    echo "WARN: winnowmap is long-read only; skipping" >&2
    return 0
  fi
  local info; info="$(lookup_mapper winnowmap)"
  [[ -z "$info" ]] && { echo "WARN: winnowmap not installed, skipping" >&2; return 0; }
  # repetitive_kmers cached next to ref
  local rep_k="$TIER_DIR/.idx/repetitive_k15.txt"
  mkdir -p "$(dirname "$rep_k")"
  if [[ ! -f "$rep_k" ]]; then
    if command -v meryl >/dev/null 2>&1; then
      echo "[winnowmap] generating repetitive k-mer list with meryl" >&2
      ( cd "$(dirname "$rep_k")" \
        && meryl count k=15 output merylDB "$REF" >/dev/null 2>&1 \
        && meryl print greater-than distinct=0.9998 merylDB > "$(basename "$rep_k")" )
    else
      echo "WARN: meryl missing; writing empty repetitive_kmers.txt (winnowmap will run without masking)" >&2
      : > "$rep_k"
    fi
  fi
  local sam="$out/aln.sam"
  echo "[winnowmap] aligning -> $bam" >&2
  local t0; t0=$(date +%s.%N)
  /usr/bin/time -v -o "$out/.time" bash -c \
    "winnowmap -W '$rep_k' -t $THREADS -ax map-pb '$REF' '$READS'" > "$sam" 2>>"$out/stderr.log"
  local t1; t1=$(date +%s.%N)
  awk -F'\t' 'NF >= 11 || /^@/' "$sam" | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  rm -f "$sam"
  emit_timing "winnowmap" "$timing" "$t0" "$t1" "$out/.time" "$bam"
}

emit_timing() {
  local mapper="$1" out_json="$2" t0="$3" t1="$4" time_log="$5" bam="$6"
  python3 - "$mapper" "$out_json" "$t0" "$t1" "$time_log" "$bam" <<'PY'
import json, os, re, sys
mapper, out_json, t0, t1, time_log, bam = sys.argv[1:7]
wall = max(0.0, float(t1) - float(t0))
peak_rss = 0
user_s = sys_s = None
if os.path.isfile(time_log):
    txt = open(time_log).read()
    m = re.search(r'Maximum resident set size \(kbytes\): (\d+)', txt)
    if m: peak_rss = int(m.group(1)) * 1024
    m = re.search(r'User time \(seconds\): ([0-9.]+)', txt)
    if m: user_s = float(m.group(1))
    m = re.search(r'System time \(seconds\): ([0-9.]+)', txt)
    if m: sys_s = float(m.group(1))
out = {
    "mapper": mapper,
    "wallclock_seconds": wall,
    "user_cpu_seconds": user_s,
    "system_cpu_seconds": sys_s,
    "peak_rss_bytes": peak_rss,
    "bam_bytes": os.path.getsize(bam) if os.path.isfile(bam) else 0,
}
with open(out_json, "w") as f:
    json.dump(out, f, indent=2)
PY
}

ensure_fai

# Dispatch the requested mappers
IFS=',' read -ra MAPPER_LIST <<< "$MAPPERS"
for m in "${MAPPER_LIST[@]}"; do
  case "$m" in
    llmap)      run_llmap     || echo "WARN: llmap run failed" >&2 ;;
    minimap2)   run_minimap2  || echo "WARN: minimap2 run failed" >&2 ;;
    bwa-mem2)   run_bwa_mem2  || echo "WARN: bwa-mem2 run failed" >&2 ;;
    winnowmap)  run_winnowmap || echo "WARN: winnowmap run failed" >&2 ;;
    *) echo "WARN: unknown mapper '$m', skipping" >&2 ;;
  esac
done

# ---- step 3: accuracy --------------------------------------------------------
python3 "$SCRIPT_DIR/analyze_bench.py" --tier-dir "$TIER_DIR" --mappers "$MAPPERS"

# ---- step 4: scoreboard ------------------------------------------------------
if [[ ! -f "$SCOREBOARD" ]]; then
  printf 'organism\ttier\tmapper\tn_total\tn_mapped\tmapping_rate\tf1_within_100bp\tf1_within_1kb\twallclock_seconds\tpeak_rss_bytes\ttimestamp\n' > "$SCOREBOARD"
fi

python3 - "$TIER_DIR" "$ORGANISM" "$TIER" "$SCOREBOARD" <<'PY'
import json, os, sys, time, pathlib
tier_dir, organism, tier, scoreboard = sys.argv[1:5]
summary_path = pathlib.Path(tier_dir) / "summary.json"
if not summary_path.is_file():
    sys.exit(0)
summary = json.loads(summary_path.read_text())
ts = time.strftime("%Y-%m-%dT%H:%M:%S%z")
with open(scoreboard, "a") as f:
    for mapper, acc in summary.get("per_mapper", {}).items():
        ov = acc["overall"]
        timing_path = pathlib.Path(tier_dir) / mapper / "rep0" / "timing.json"
        wall = rss = 0
        if timing_path.is_file():
            t = json.loads(timing_path.read_text())
            wall = t.get("wallclock_seconds", 0) or 0
            rss = t.get("peak_rss_bytes", 0) or 0
        f.write(f"{organism}\t{tier}\t{mapper}\t{ov['n_total']}\t{ov['n_mapped']}\t"
                f"{ov['mapping_rate']:.6f}\t{ov['f1_within_100bp']:.6f}\t"
                f"{ov['f1_within_1kb']:.6f}\t{wall:.3f}\t{rss}\t{ts}\n")
print(f"[run_species_bench] appended scoreboard rows for {organism}/tier{tier}", file=sys.stderr)
PY

echo "[run_species_bench] DONE  $ORGANISM tier=$TIER  -> $TIER_DIR" >&2
