#!/usr/bin/env python3
"""Unit tests for compute.py metrics extraction.

Creates mini BAM files in a temp directory and validates that compute.py
correctly extracts mapping metrics from them.
"""

import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Generator

import pysam
import pytest


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def reference_path(temp_dir: Path) -> Path:
    """Create a minimal reference FASTA."""
    ref = temp_dir / "ref.fa"
    ref.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
    return ref


@pytest.fixture
def indexed_reference(reference_path: Path, temp_dir: Path) -> Path:
    """Create a reference with FAI index for BAM creation."""
    subprocess.run(["samtools", "faidx", str(reference_path)], check=True)
    return reference_path


def create_bam(
    bam_path: Path,
    ref_path: Path,
    alignments: list[dict],
) -> Path:
    """Create a BAM file with the given alignments.

    Each alignment is a dict with keys:
        name: read name
        flag: SAM flag (0=mapped forward, 4=unmapped, 256=secondary, 2048=supplementary)
        ref_name: reference name (or "*" for unmapped)
        pos: 0-based position
        mapq: mapping quality
        cigar: CIGAR string
        seq: sequence
        qual: quality string (same length as seq)
    """
    header = pysam.AlignmentHeader.from_references(
        ["chr1"], [40]
    )
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for aln in alignments:
            a = pysam.AlignedSegment()
            a.query_name = aln["name"]
            a.flag = aln["flag"]
            if aln["flag"] & 4:  # unmapped
                a.reference_id = -1
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigarstring = None
            else:
                a.reference_id = 0  # chr1
                a.reference_start = aln["pos"]
                a.mapping_quality = aln["mapq"]
                a.cigarstring = aln["cigar"]
            a.query_sequence = aln["seq"]
            a.query_qualities = pysam.qualitystring_to_array(aln["qual"])
            out.write(a)

    # Sort and index
    sorted_bam = bam_path.with_suffix(".sorted.bam")
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return sorted_bam


def run_compute(
    bam_path: Path,
    out_dir: Path,
    total: int,
    truth_path: Path | None = None,
    tolerance_bp: int = 10,
) -> dict:
    """Run compute.py and return the mapping_summary.json contents."""
    script = Path(__file__).parent / "compute.py"
    cmd = [
        "python3", str(script),
        "--bam", str(bam_path),
        "--out-dir", str(out_dir),
        "--total", str(total),
        "--tolerance-bp", str(tolerance_bp),
    ]
    if truth_path:
        cmd.extend(["--truth", str(truth_path)])

    subprocess.run(cmd, check=True)
    return json.loads((out_dir / "mapping_summary.json").read_text())


class TestSummarize:
    """Tests for the summarize function via compute.py main."""

    def test_all_mapped(self, temp_dir: Path) -> None:
        """All reads mapped with high MAPQ."""
        alignments = [
            {"name": f"read{i}", "flag": 0, "pos": i * 2, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"}
            for i in range(5)
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)
        out = temp_dir / "out"

        summary = run_compute(bam, out, total=5)

        assert summary["total_input_reads"] == 5
        assert summary["primary_mapped"] == 5
        assert summary["unmapped"] == 0
        assert summary["secondary_count"] == 0
        assert summary["supplementary_count"] == 0
        assert summary["mapping_rate"] == 1.0
        assert summary["drop_rate"] == 0.0

    def test_all_unmapped(self, temp_dir: Path) -> None:
        """All reads unmapped."""
        alignments = [
            {"name": f"read{i}", "flag": 4, "pos": 0, "mapq": 0,
             "cigar": "*", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"}
            for i in range(3)
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)
        out = temp_dir / "out"

        summary = run_compute(bam, out, total=3)

        assert summary["total_input_reads"] == 3
        assert summary["primary_mapped"] == 0
        assert summary["unmapped"] == 3
        assert summary["mapping_rate"] == 0.0
        assert summary["drop_rate"] == 1.0

    def test_mixed_alignments(self, temp_dir: Path) -> None:
        """Mix of primary, secondary, supplementary, unmapped."""
        alignments = [
            # 2 primary mapped
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read2", "flag": 0, "pos": 10, "mapq": 30,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            # 1 secondary
            {"name": "read1", "flag": 256, "pos": 20, "mapq": 0,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            # 1 supplementary
            {"name": "read2", "flag": 2048, "pos": 5, "mapq": 10,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            # 1 unmapped
            {"name": "read3", "flag": 4, "pos": 0, "mapq": 0,
             "cigar": "*", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)
        out = temp_dir / "out"

        summary = run_compute(bam, out, total=3)

        assert summary["primary_mapped"] == 2
        assert summary["secondary_count"] == 1
        assert summary["supplementary_count"] == 1
        assert summary["unmapped"] == 1
        # mapping_rate is primary/total
        assert summary["mapping_rate"] == pytest.approx(2 / 3, rel=1e-6)

    def test_mapq_histogram(self, temp_dir: Path) -> None:
        """MAPQ histogram is correctly generated."""
        alignments = [
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read2", "flag": 0, "pos": 10, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read3", "flag": 0, "pos": 20, "mapq": 30,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read4", "flag": 0, "pos": 5, "mapq": 0,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)
        out = temp_dir / "out"

        run_compute(bam, out, total=4)
        hist = json.loads((out / "mapq_histogram.json").read_text())

        assert hist["60"] == 2
        assert hist["30"] == 1
        assert hist["0"] == 1

    def test_empty_bam(self, temp_dir: Path) -> None:
        """Empty BAM file handled gracefully."""
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", [])
        out = temp_dir / "out"

        summary = run_compute(bam, out, total=0)

        assert summary["primary_mapped"] == 0
        assert summary["unmapped"] == 0


class TestGroundTruth:
    """Tests for ground truth evaluation."""

    def test_perfect_truth(self, temp_dir: Path) -> None:
        """All reads map to correct positions."""
        alignments = [
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read2", "flag": 0, "pos": 10, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            {"name": "read3", "flag": 0, "pos": 20, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text(
            "read1\tchr1\t0\n"
            "read2\tchr1\t10\n"
            "read3\tchr1\t20\n"
        )

        out = temp_dir / "out"
        run_compute(bam, out, total=3, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 3
        assert gt["false_positive"] == 0
        assert gt["false_negative"] == 0
        assert gt["recall"] == 1.0
        assert gt["precision"] == 1.0
        assert gt["f1"] == 1.0

    def test_position_tolerance(self, temp_dir: Path) -> None:
        """Reads within tolerance are true positives, outside are FP.

        recall = TP / (TP + FN), precision = TP / (TP + FP)
        Here: read1 at pos=5 vs truth=0, |5-0|=5 <= 10 -> TP
              read2 at pos=25 vs truth=10, |25-10|=15 > 10 -> FP
        FN = 0 (both truth reads have primary alignments)
        So: recall = 1/(1+0) = 1.0, precision = 1/(1+1) = 0.5
        """
        alignments = [
            # Mapped 5bp away from truth (within default 10bp tolerance)
            {"name": "read1", "flag": 0, "pos": 5, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
            # Mapped 15bp away (outside tolerance)
            {"name": "read2", "flag": 0, "pos": 25, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text(
            "read1\tchr1\t0\n"
            "read2\tchr1\t10\n"
        )

        out = temp_dir / "out"
        run_compute(bam, out, total=2, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 1
        assert gt["false_positive"] == 1
        assert gt["false_negative"] == 0
        assert gt["recall"] == 1.0  # TP/(TP+FN) = 1/(1+0)
        assert gt["precision"] == 0.5  # TP/(TP+FP) = 1/(1+1)

    def test_wrong_chromosome(self, temp_dir: Path) -> None:
        """Read on wrong chromosome is false positive."""
        # We only have chr1 in our test setup, so simulate by truth saying chr2
        alignments = [
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text("read1\tchr2\t0\n")

        out = temp_dir / "out"
        run_compute(bam, out, total=1, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 0
        assert gt["false_positive"] == 1

    def test_unmapped_in_truth(self, temp_dir: Path) -> None:
        """Read in truth but unmapped is false negative."""
        alignments = [
            {"name": "read1", "flag": 4, "pos": 0, "mapq": 0,
             "cigar": "*", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text("read1\tchr1\t0\n")

        out = temp_dir / "out"
        run_compute(bam, out, total=1, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 0
        assert gt["false_negative"] == 1

    def test_missing_from_bam(self, temp_dir: Path) -> None:
        """Read in truth but missing from BAM is false negative."""
        alignments = [
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text(
            "read1\tchr1\t0\n"
            "read2\tchr1\t10\n"  # This read is not in BAM
        )

        out = temp_dir / "out"
        run_compute(bam, out, total=2, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 1
        assert gt["false_negative"] == 1

    def test_custom_tolerance(self, temp_dir: Path) -> None:
        """Custom tolerance parameter works."""
        alignments = [
            {"name": "read1", "flag": 0, "pos": 5, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text("read1\tchr1\t0\n")

        out = temp_dir / "out"
        # With tolerance 3, pos=5 is outside (|5-0|=5 > 3)
        run_compute(bam, out, total=1, truth_path=truth_path, tolerance_bp=3)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 0
        assert gt["false_positive"] == 1
        assert gt["tolerance_bp"] == 3

    def test_truth_with_comments(self, temp_dir: Path) -> None:
        """Truth file with comments and blank lines."""
        alignments = [
            {"name": "read1", "flag": 0, "pos": 0, "mapq": 60,
             "cigar": "10M", "seq": "ACGTACGTAC", "qual": "IIIIIIIIII"},
        ]
        bam = create_bam(temp_dir / "test.bam", temp_dir / "ref.fa", alignments)

        truth_path = temp_dir / "truth.tsv"
        truth_path.write_text(
            "# This is a comment\n"
            "\n"
            "read1\tchr1\t0\n"
            "   \n"
        )

        out = temp_dir / "out"
        run_compute(bam, out, total=1, truth_path=truth_path)

        gt = json.loads((out / "ground_truth.json").read_text())
        assert gt["true_positive"] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
