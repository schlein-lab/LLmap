#!/usr/bin/env bash
# Verify that all required tools are reachable at pinned versions.
# Writes benchmarks/reports/tool_versions.json and exits non-zero if any mismatch.
#
# Usage:
#   ./check_versions.sh [--json-only] [output.json]
#
# Exit codes:
#   0 = all tools present and versions match
#   1 = at least one tool missing or version mismatch
#   2 = script error
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-${SCRIPT_DIR}/../reports/tool_versions.json}"
JSON_ONLY=false

if [[ "${1:-}" == "--json-only" ]]; then
  JSON_ONLY=true
  OUT="${2:-${SCRIPT_DIR}/../reports/tool_versions.json}"
fi

mkdir -p "$(dirname "$OUT")"

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m'

declare -a RESULTS=()
MISMATCH_COUNT=0
MISSING_COUNT=0

# Extract version number from tool output
# Args: raw_output expected_pattern
extract_version() {
  local raw="$1"
  # Common patterns: "2.28-r...", "version 2.2.1", "STAR_2.7.11b", "bowtie2-2.5.4"
  echo "$raw" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?[a-z]?' | head -1
}

# Check a single tool
# Args: name cmd expected_version
check_tool() {
  local name="$1"
  local cmd="$2"
  local expected="$3"

  local raw_output=""
  local got=""
  local bin=""
  local sha=""
  local status="ok"

  # Try to run the command
  if raw_output="$(eval "$cmd" 2>&1)"; then
    got="$(extract_version "$raw_output")"
    bin="$(command -v "${cmd%% *}" 2>/dev/null || true)"
    [[ -n "$bin" && -f "$bin" ]] && sha="$(sha256sum "$bin" | awk '{print $1}')"

    if [[ -z "$got" ]]; then
      status="parse_error"
      got="$raw_output"
    elif [[ "$got" != "$expected" ]]; then
      status="mismatch"
      ((MISMATCH_COUNT++)) || true
    fi
  else
    status="missing"
    got="NOT_FOUND"
    ((MISSING_COUNT++)) || true
  fi

  # Store result for JSON
  RESULTS+=("$(printf '    "%s": {"command": "%s", "reported": "%s", "expected": "%s", "binary": "%s", "sha256": "%s", "status": "%s"}' \
    "$name" "$cmd" "$got" "$expected" "${bin:-}" "${sha:-}" "$status")")

  # Print status unless JSON-only
  if [[ "$JSON_ONLY" == "false" ]]; then
    case "$status" in
      ok)
        printf "  ${GREEN}[OK]${NC} %-12s %s\n" "$name" "$got"
        ;;
      mismatch)
        printf "  ${YELLOW}[MISMATCH]${NC} %-12s got %s, expected %s\n" "$name" "$got" "$expected"
        ;;
      missing)
        printf "  ${RED}[MISSING]${NC} %-12s command not found: %s\n" "$name" "$cmd"
        ;;
      parse_error)
        printf "  ${YELLOW}[PARSE]${NC} %-12s could not parse version from: %s\n" "$name" "${got:0:50}..."
        ;;
    esac
  fi
}

# Header
if [[ "$JSON_ONLY" == "false" ]]; then
  echo "Checking benchmark tool versions..."
  echo ""
fi

# Check all tools from tools.yaml
check_tool "minimap2"  "minimap2 --version"                 "2.28"
check_tool "bwa_mem2"  "bwa-mem2 version 2>&1"              "2.2.1"
check_tool "winnowmap" "winnowmap --version 2>&1"          "2.03"
check_tool "star"      "STAR --version 2>&1"               "2.7.11b"
check_tool "bowtie2"   "bowtie2 --version 2>&1"            "2.5.4"
check_tool "samtools"  "samtools --version 2>&1"           "1.21"
check_tool "seqkit"    "seqkit version 2>&1"               "2.8.2"
check_tool "llmap"     "${LLMAP_BIN:-llmap} --version 2>&1" "1.0.0"

# Write JSON output
{
  echo "{"
  IFS=","
  for i in "${!RESULTS[@]}"; do
    if [[ $i -lt $((${#RESULTS[@]} - 1)) ]]; then
      echo "${RESULTS[$i]},"
    else
      echo "${RESULTS[$i]}"
    fi
  done
  echo "  ,\"_check_run_at\": \"$(date -Iseconds)\""
  echo "  ,\"_mismatch_count\": $MISMATCH_COUNT"
  echo "  ,\"_missing_count\": $MISSING_COUNT"
  echo "}"
} > "$OUT"

# Summary
if [[ "$JSON_ONLY" == "false" ]]; then
  echo ""
  echo "Results written to: $OUT"
  echo ""
fi

# Exit code
TOTAL_ERRORS=$((MISMATCH_COUNT + MISSING_COUNT))
if [[ $TOTAL_ERRORS -gt 0 ]]; then
  if [[ "$JSON_ONLY" == "false" ]]; then
    echo "${RED}FAIL${NC}: $MISSING_COUNT missing, $MISMATCH_COUNT version mismatch(es)"
    echo ""
    echo "Run 'benchmarks/runners/install_tools.sh' to install missing tools."
  fi
  exit 1
fi

if [[ "$JSON_ONLY" == "false" ]]; then
  echo "${GREEN}PASS${NC}: all tools present at expected versions"
fi
exit 0
