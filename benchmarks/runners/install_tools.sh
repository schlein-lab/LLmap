#!/usr/bin/env bash
# Install benchmark tools via conda or module load.
#
# Usage:
#   ./install_tools.sh [--conda-only | --module-only] [--env <name>]
#
# Options:
#   --conda-only   Only install via conda (skip module checks)
#   --module-only  Only load via modules (skip conda)
#   --env <name>   Conda environment name (default: llmap-bench)
#   --dry-run      Print commands without executing
#   --help         Show this help
#
# On Hummel-2, prefer module load; conda is the fallback for local machines.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="llmap-bench"
MODE="auto"
DRY_RUN=false

usage() {
  grep '^#' "$0" | grep -v '^#!/' | sed 's/^# //'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --conda-only)  MODE="conda"; shift ;;
    --module-only) MODE="module"; shift ;;
    --env)         ENV_NAME="$2"; shift 2 ;;
    --dry-run)     DRY_RUN=true; shift ;;
    --help|-h)     usage ;;
    *)             echo "Unknown option: $1"; exit 2 ;;
  esac
done

run_cmd() {
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[dry-run] $*"
  else
    echo "+ $*"
    eval "$@"
  fi
}

# Detect environment
HAS_MODULE=false
HAS_CONDA=false

if command -v module &>/dev/null; then
  HAS_MODULE=true
fi

if command -v conda &>/dev/null || command -v mamba &>/dev/null; then
  HAS_CONDA=true
fi

echo "=== LLmap Benchmark Tool Installation ==="
echo ""
echo "Environment: ${ENV_NAME}"
echo "Mode: ${MODE}"
echo "Module system: ${HAS_MODULE}"
echo "Conda: ${HAS_CONDA}"
echo ""

# Module-based installation (Hummel-2 and similar HPC clusters)
install_via_modules() {
  echo "=== Loading tools via module system ==="

  local modules=(
    "minimap2/2.28"
    "bwa-mem2/2.2.1"
    "winnowmap/2.03"
    "star/2.7.11b"
    "bowtie2/2.5.4"
    "samtools/1.21"
    "htslib/1.21"
  )

  for mod in "${modules[@]}"; do
    if module avail "$mod" 2>&1 | grep -q "$mod"; then
      run_cmd "module load $mod"
    else
      echo "  [SKIP] Module not available: $mod"
    fi
  done

  echo ""
  echo "Note: Module loads are session-local. Add to your ~/.bashrc or job script:"
  echo ""
  for mod in "${modules[@]}"; do
    echo "  module load $mod"
  done
}

# Conda-based installation
install_via_conda() {
  echo "=== Installing tools via conda ==="

  local CONDA_CMD="conda"
  if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba (faster solver)"
  fi

  # Create environment if it doesn't exist
  if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment: ${ENV_NAME}"
    run_cmd "$CONDA_CMD create -n ${ENV_NAME} -y python=3.11"
  fi

  echo "Installing tools into: ${ENV_NAME}"

  # Install bioconda tools
  local tools=(
    "bioconda::minimap2=2.28"
    "bioconda::bwa-mem2=2.2.1"
    "bioconda::winnowmap=2.03"
    "bioconda::star=2.7.11b"
    "bioconda::bowtie2=2.5.4"
    "bioconda::samtools=1.21"
    "bioconda::htslib=1.21"
    "bioconda::seqkit=2.8.2"
  )

  run_cmd "$CONDA_CMD install -n ${ENV_NAME} -c bioconda -c conda-forge -y ${tools[*]}"

  echo ""
  echo "Installation complete. Activate with:"
  echo "  conda activate ${ENV_NAME}"
}

# LLmap installation reminder
install_llmap_note() {
  echo ""
  echo "=== LLmap ==="
  echo "LLmap is built from source. Ensure it's in your PATH:"
  echo ""
  echo "  export LLMAP_ROOT=/path/to/llmap"
  echo '  export LLMAP_BIN="${LLMAP_ROOT}/build/llmap"'
  echo '  export PATH="${LLMAP_ROOT}/build:${PATH}"'
  echo ""
  echo "To build LLmap:"
  echo "  cd \$LLMAP_ROOT"
  echo "  mkdir -p build && cd build"
  echo "  cmake -DLLMAP_ENABLE_CUDA=ON .."
  echo "  cmake --build . -j\$(nproc)"
}

# Main logic
case "$MODE" in
  auto)
    if [[ "$HAS_MODULE" == "true" ]]; then
      install_via_modules
    elif [[ "$HAS_CONDA" == "true" ]]; then
      install_via_conda
    else
      echo "ERROR: Neither module nor conda available."
      echo "Install miniconda: https://docs.conda.io/en/latest/miniconda.html"
      exit 1
    fi
    ;;
  module)
    if [[ "$HAS_MODULE" != "true" ]]; then
      echo "ERROR: module command not available"
      exit 1
    fi
    install_via_modules
    ;;
  conda)
    if [[ "$HAS_CONDA" != "true" ]]; then
      echo "ERROR: conda/mamba not available"
      exit 1
    fi
    install_via_conda
    ;;
esac

install_llmap_note

echo ""
echo "=== Verification ==="
echo "Run version check:"
echo "  ${SCRIPT_DIR}/check_versions.sh"
echo ""
