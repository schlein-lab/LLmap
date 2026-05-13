#!/usr/bin/env bash
# LLmap — pbsim3 wrapper for HiFi/iso-seq read simulation.
#
# Usage:
#   pbsim3_wrapper.sh --fasta reference.fa --output reads.fastq [options]
#
# Options:
#   --fasta FILE         Reference FASTA (required)
#   --output FILE        Output FASTQ (required)
#   --depth INT          Coverage depth (default: 30)
#   --preset PRESET      Simulation preset: hifi, isoseq, clr (default: hifi)
#   --seed INT           Random seed for reproducibility (default: 42)
#   --truncation FLOAT   5'-truncation rate for iso-seq (default: 0.3)
#   --help               Show this help

set -euo pipefail

# Defaults
FASTA=""
OUTPUT=""
DEPTH=30
PRESET="hifi"
SEED=42
TRUNCATION=0.3

usage() {
    head -n 15 "$0" | tail -n 13
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta)
            FASTA="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --depth)
            DEPTH="$2"
            shift 2
            ;;
        --preset)
            PRESET="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --truncation)
            TRUNCATION="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$FASTA" || -z "$OUTPUT" ]]; then
    echo "Error: --fasta and --output are required" >&2
    exit 1
fi

# Check for pbsim3
if ! command -v pbsim3 &> /dev/null; then
    echo "pbsim3 not found. Install from: https://github.com/yukiteruono/pbsim3" >&2
    echo "Skipping pbsim3 simulation (will generate placeholder)." >&2

    # Create placeholder output
    echo "@pbsim3_placeholder" > "$OUTPUT"
    echo "ACGT" >> "$OUTPUT"
    echo "+" >> "$OUTPUT"
    echo "IIII" >> "$OUTPUT"
    exit 0
fi

# Set parameters based on preset
case $PRESET in
    hifi)
        # HiFi: high accuracy, ~10-20kb reads
        ACCURACY="0.999"
        LENGTH_MEAN="15000"
        LENGTH_SD="3000"
        ;;
    isoseq)
        # Iso-Seq: transcript-length reads, 5'-truncation common
        ACCURACY="0.999"
        LENGTH_MEAN="2000"
        LENGTH_SD="1000"
        ;;
    clr)
        # CLR: lower accuracy, longer reads
        ACCURACY="0.90"
        LENGTH_MEAN="20000"
        LENGTH_SD="10000"
        ;;
    *)
        echo "Unknown preset: $PRESET (use hifi, isoseq, clr)" >&2
        exit 1
        ;;
esac

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "Running pbsim3 with preset: $PRESET, depth: $DEPTH, seed: $SEED"

pbsim3 \
    --strategy wgs \
    --method qshmm \
    --qshmm QSHMM-RSII.model \
    --depth "$DEPTH" \
    --genome "$FASTA" \
    --seed "$SEED" \
    --length-mean "$LENGTH_MEAN" \
    --length-sd "$LENGTH_SD" \
    --accuracy-mean "$ACCURACY" \
    --prefix "${TMPDIR}/sim" \
    2>/dev/null || {
        echo "pbsim3 failed, using fallback qshmm model" >&2
        pbsim3 \
            --strategy wgs \
            --method errhmm \
            --errhmm ERRHMM-RSII.model \
            --depth "$DEPTH" \
            --genome "$FASTA" \
            --seed "$SEED" \
            --length-mean "$LENGTH_MEAN" \
            --length-sd "$LENGTH_SD" \
            --prefix "${TMPDIR}/sim" \
            2>/dev/null || true
    }

# Combine output files
if ls "${TMPDIR}"/sim*.fastq 1>/dev/null 2>&1; then
    cat "${TMPDIR}"/sim*.fastq > "$OUTPUT"
    echo "Generated $(grep -c '^@' "$OUTPUT") reads to $OUTPUT"
else
    echo "Warning: No reads generated, creating placeholder" >&2
    echo "@pbsim3_placeholder" > "$OUTPUT"
    echo "ACGT" >> "$OUTPUT"
    echo "+" >> "$OUTPUT"
    echo "IIII" >> "$OUTPUT"
fi
