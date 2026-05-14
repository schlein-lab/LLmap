#!/usr/bin/env bash
# Verify that all required tools are reachable at pinned versions.
# Writes benchmarks/reports/tool_versions.json.
set -euo pipefail

OUT="${1:-$(dirname "$0")/../reports/tool_versions.json}"
mkdir -p "$(dirname "$OUT")"

record() {
  local key="$1" cmd="$2" expect="$3"
  local got
  got="$(eval "$cmd" 2>&1 | head -1 | tr -d '\r')" || got="MISSING"
  local sha=""
  local bin
  bin="$(command -v "${cmd%% *}" 2>/dev/null || true)"
  [[ -n "$bin" && -f "$bin" ]] && sha="$(sha256sum "$bin" | awk '{print $1}')"
  echo "  \"$key\": {\"command\": \"$cmd\", \"reported\": \"$got\", \"expected\": \"$expect\", \"binary\": \"$bin\", \"sha256\": \"$sha\"},"
}

{
  echo "{"
  record minimap2  "minimap2 --version"   "2.28"
  record bwa_mem2  "bwa-mem2 version"     "2.2.1"
  record winnowmap "winnowmap --version"  "2.03"
  record star      "STAR --version"       "2.7.11b"
  record bowtie2   "bowtie2 --version"    "2.5.4"
  record samtools  "samtools --version"   "1.21"
  record llmap     "${LLMAP_BIN:-llmap} --version"  "1.0.0"
  echo "  \"_check_run_at\": \"$(date -Iseconds)\""
  echo "}"
} > "$OUT"

echo "wrote $OUT"
