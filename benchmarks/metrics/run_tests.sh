#!/usr/bin/env bash
# Run unit tests for benchmark metrics modules.
#
# Usage:
#   ./run_tests.sh           # run all tests
#   ./run_tests.sh -v        # verbose output
#   ./run_tests.sh compute   # run only compute tests
#   ./run_tests.sh concordance # run only concordance tests

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

VERBOSE=""
TARGET="all"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE="-v"
            shift
            ;;
        compute|concordance|all)
            TARGET="$1"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [-v|--verbose] [compute|concordance|all]"
            echo ""
            echo "Options:"
            echo "  -v, --verbose    Verbose test output"
            echo "  compute          Run compute.py tests only"
            echo "  concordance      Run concordance.py tests only"
            echo "  all              Run all tests (default)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Check dependencies
if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 not found" >&2
    exit 1
fi

if ! python3 -c "import pytest" 2>/dev/null; then
    echo "ERROR: pytest not installed. Run: pip install pytest" >&2
    exit 1
fi

if ! python3 -c "import pysam" 2>/dev/null; then
    echo "ERROR: pysam not installed. Run: pip install pysam" >&2
    exit 1
fi

if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools not found" >&2
    exit 1
fi

echo "=== Benchmark Metrics Unit Tests ==="
echo ""

case "$TARGET" in
    compute)
        echo "Running compute.py tests..."
        python3 -m pytest test_compute.py $VERBOSE
        ;;
    concordance)
        echo "Running concordance.py tests..."
        python3 -m pytest test_concordance.py $VERBOSE
        ;;
    all)
        echo "Running all tests..."
        python3 -m pytest test_compute.py test_concordance.py $VERBOSE
        ;;
esac

echo ""
echo "=== All tests passed ==="
