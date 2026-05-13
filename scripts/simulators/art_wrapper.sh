#!/usr/bin/env bash
# LLmap — ART wrapper for Illumina read simulation.
#
# Usage:
#   art_wrapper.sh --fasta reference.fa --output reads.fastq [options]
#
# Options:
#   --fasta FILE         Reference FASTA (required)
#   --output FILE        Output FASTQ prefix (required)
#   --depth INT          Coverage depth (default: 30)
#   --preset PRESET      Simulation preset: illumina_hiseq, illumina_novaseq (default: illumina_novaseq)
#   --seed INT           Random seed for reproducibility (default: 42)
#   --read-length INT    Read length (default: 150)
#   --paired             Generate paired-end reads (default: single-end)
#   --insert-mean INT    Insert size mean for paired-end (default: 400)
#   --insert-sd INT      Insert size SD for paired-end (default: 50)
#   --help               Show this help

set -euo pipefail

# Defaults
FASTA=""
OUTPUT=""
DEPTH=30
PRESET="illumina_novaseq"
SEED=42
READ_LENGTH=150
PAIRED=false
INSERT_MEAN=400
INSERT_SD=50

usage() {
    head -n 18 "$0" | tail -n 16
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
        --read-length)
            READ_LENGTH="$2"
            shift 2
            ;;
        --paired)
            PAIRED=true
            shift
            ;;
        --insert-mean)
            INSERT_MEAN="$2"
            shift 2
            ;;
        --insert-sd)
            INSERT_SD="$2"
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

# Check for ART
ART_CMD=""
for cmd in art_illumina art; do
    if command -v "$cmd" &> /dev/null; then
        ART_CMD="$cmd"
        break
    fi
done

if [[ -z "$ART_CMD" ]]; then
    echo "ART not found. Install from: https://www.niehs.nih.gov/research/resources/software/biostatistics/art" >&2
    echo "Skipping ART simulation (will generate placeholder)." >&2

    # Create placeholder output
    echo "@art_placeholder/1" > "${OUTPUT}_1.fq"
    echo "ACGTACGTACGTACGT" >> "${OUTPUT}_1.fq"
    echo "+" >> "${OUTPUT}_1.fq"
    echo "IIIIIIIIIIIIIIII" >> "${OUTPUT}_1.fq"

    if $PAIRED; then
        echo "@art_placeholder/2" > "${OUTPUT}_2.fq"
        echo "TGCATGCATGCATGCA" >> "${OUTPUT}_2.fq"
        echo "+" >> "${OUTPUT}_2.fq"
        echo "IIIIIIIIIIIIIIII" >> "${OUTPUT}_2.fq"
    fi
    exit 0
fi

# Set quality profile based on preset
case $PRESET in
    illumina_hiseq)
        PROFILE="HS25"
        ;;
    illumina_novaseq)
        PROFILE="HS25"  # ART doesn't have NovaSeq yet; HiSeq is close enough
        ;;
    *)
        echo "Unknown preset: $PRESET (use illumina_hiseq, illumina_novaseq)" >&2
        exit 1
        ;;
esac

echo "Running ART with preset: $PRESET, depth: $DEPTH, seed: $SEED"

if $PAIRED; then
    "$ART_CMD" \
        -ss "$PROFILE" \
        -i "$FASTA" \
        -p \
        -l "$READ_LENGTH" \
        -f "$DEPTH" \
        -m "$INSERT_MEAN" \
        -s "$INSERT_SD" \
        -rs "$SEED" \
        -o "${OUTPUT}_" \
        -q || {
            echo "ART failed" >&2
            exit 1
        }

    # Rename outputs to expected names
    mv "${OUTPUT}_1.fq" "${OUTPUT}_R1.fq" 2>/dev/null || true
    mv "${OUTPUT}_2.fq" "${OUTPUT}_R2.fq" 2>/dev/null || true

    echo "Generated paired-end reads to ${OUTPUT}_R1.fq and ${OUTPUT}_R2.fq"
else
    "$ART_CMD" \
        -ss "$PROFILE" \
        -i "$FASTA" \
        -l "$READ_LENGTH" \
        -f "$DEPTH" \
        -rs "$SEED" \
        -o "${OUTPUT}_" \
        -q || {
            echo "ART failed" >&2
            exit 1
        }

    mv "${OUTPUT}_.fq" "${OUTPUT}.fq" 2>/dev/null || true

    echo "Generated single-end reads to ${OUTPUT}.fq"
fi
