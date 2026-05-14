#!/usr/bin/env python3
"""Unit tests for concordance.py pairwise BAM comparison.

Creates pairs of mini BAM files and validates classification of agreement
between two aligners.
"""

import json
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


def create_bam(
    bam_path: Path,
    alignments: list[dict],
    refs: list[tuple[str, int]] | None = None,
) -> Path:
    """Create a BAM file with the given alignments.

    Each alignment is a dict with keys:
        name: read name
        flag: SAM flag (0=mapped forward, 4=unmapped, 256=secondary, 2048=supplementary)
        ref_name: reference name (or "*" for unmapped)
        ref_id: reference index (0-based)
        pos: 0-based position
        mapq: mapping quality
        cigar: CIGAR string
        seq: sequence
        qual: quality string
    """
    if refs is None:
        refs = [("chr1", 1000), ("chr2", 1000)]

    header = pysam.AlignmentHeader.from_references(
        [r[0] for r in refs], [r[1] for r in refs]
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
                a.reference_id = aln.get("ref_id", 0)
                a.reference_start = aln["pos"]
                a.mapping_quality = aln.get("mapq", 60)
                a.cigarstring = aln.get("cigar", "10M")
            a.query_sequence = aln.get("seq", "ACGTACGTAC")
            a.query_qualities = pysam.qualitystring_to_array(
                aln.get("qual", "IIIIIIIIII")
            )
            out.write(a)

    sorted_bam = bam_path.with_suffix(".sorted.bam")
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return sorted_bam


def run_concordance(
    bam_a: Path,
    bam_b: Path,
    out_dir: Path,
    name_a: str = "toolA",
    name_b: str = "toolB",
    tolerance_bp: int = 10,
    max_examples: int = 100,
) -> dict:
    """Run concordance.py and return the concordance.json contents."""
    script = Path(__file__).parent / "concordance.py"
    cmd = [
        "python3", str(script),
        "--bam-a", str(bam_a),
        "--bam-b", str(bam_b),
        "--name-a", name_a,
        "--name-b", name_b,
        "--out-dir", str(out_dir),
        "--tolerance-bp", str(tolerance_bp),
        "--max-examples", str(max_examples),
    ]
    subprocess.run(cmd, check=True)
    return json.loads((out_dir / "concordance.json").read_text())


class TestClassification:
    """Tests for read classification logic."""

    def test_same_position(self, temp_dir: Path) -> None:
        """Both aligners map to the same position."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["same_position"] == 1
        assert result["total_reads"] == 1

    def test_same_position_within_tolerance(self, temp_dir: Path) -> None:
        """Positions within tolerance are same_position."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 105},  # 5bp away
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out, tolerance_bp=10)

        assert result["counts"]["same_position"] == 1

    def test_same_chrom_diff_pos(self, temp_dir: Path) -> None:
        """Same chromosome but positions differ beyond tolerance."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 200},  # 100bp away
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out, tolerance_bp=10)

        assert result["counts"]["same_chrom_diff_pos"] == 1

    def test_diff_chrom(self, temp_dir: Path) -> None:
        """Different chromosomes."""
        alignments_a = [
            {"name": "read1", "flag": 0, "ref_id": 0, "pos": 100},  # chr1
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "ref_id": 1, "pos": 100},  # chr2
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["diff_chrom"] == 1

    def test_one_a_only(self, temp_dir: Path) -> None:
        """Mapped by A only (B unmapped)."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 4, "pos": 0},  # unmapped
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["oneA"] == 1

    def test_one_b_only(self, temp_dir: Path) -> None:
        """Mapped by B only (A unmapped)."""
        alignments_a = [
            {"name": "read1", "flag": 4, "pos": 0},  # unmapped
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["oneB"] == 1

    def test_both_unmapped(self, temp_dir: Path) -> None:
        """Both aligners unmapped."""
        alignments_a = [
            {"name": "read1", "flag": 4, "pos": 0},
        ]
        alignments_b = [
            {"name": "read1", "flag": 4, "pos": 0},
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["both_unmapped"] == 1

    def test_read_only_in_a(self, temp_dir: Path) -> None:
        """Read present only in BAM A."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
            {"name": "read2", "flag": 0, "pos": 200},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},
            # read2 missing
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["same_position"] == 1
        assert result["counts"]["oneA"] == 1

    def test_read_only_in_b(self, temp_dir: Path) -> None:
        """Read present only in BAM B."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},
            {"name": "read2", "flag": 0, "pos": 200},
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["same_position"] == 1
        assert result["counts"]["oneB"] == 1


class TestMixedScenarios:
    """Tests with multiple reads of different classifications."""

    def test_mixed_all_categories(self, temp_dir: Path) -> None:
        """Multiple reads across all classification categories."""
        alignments_a = [
            {"name": "same_pos", "flag": 0, "pos": 100, "ref_id": 0},
            {"name": "same_chrom", "flag": 0, "pos": 100, "ref_id": 0},
            {"name": "diff_chrom", "flag": 0, "pos": 100, "ref_id": 0},
            {"name": "only_a", "flag": 0, "pos": 100, "ref_id": 0},
            {"name": "only_b", "flag": 4, "pos": 0, "ref_id": 0},
            {"name": "both_unmap", "flag": 4, "pos": 0, "ref_id": 0},
        ]
        alignments_b = [
            {"name": "same_pos", "flag": 0, "pos": 105, "ref_id": 0},  # within tol
            {"name": "same_chrom", "flag": 0, "pos": 500, "ref_id": 0},  # outside tol
            {"name": "diff_chrom", "flag": 0, "pos": 100, "ref_id": 1},  # chr2
            {"name": "only_a", "flag": 4, "pos": 0, "ref_id": 0},  # unmapped
            {"name": "only_b", "flag": 0, "pos": 100, "ref_id": 0},  # mapped
            {"name": "both_unmap", "flag": 4, "pos": 0, "ref_id": 0},  # unmapped
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out, tolerance_bp=10)

        assert result["counts"]["same_position"] == 1
        assert result["counts"]["same_chrom_diff_pos"] == 1
        assert result["counts"]["diff_chrom"] == 1
        assert result["counts"]["oneA"] == 1
        assert result["counts"]["oneB"] == 1
        assert result["counts"]["both_unmapped"] == 1
        assert result["total_reads"] == 6

    def test_percentages_correct(self, temp_dir: Path) -> None:
        """Percentages sum to 100."""
        alignments_a = [
            {"name": f"read{i}", "flag": 0, "pos": i * 10}
            for i in range(10)
        ]
        alignments_b = [
            {"name": f"read{i}", "flag": 0, "pos": i * 10}
            for i in range(10)
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        total_pct = sum(result["percentages"].values())
        assert total_pct == pytest.approx(100.0, rel=1e-3)


class TestDivergenceExamples:
    """Tests for divergence_examples.tsv output."""

    def test_examples_written(self, temp_dir: Path) -> None:
        """Divergent reads are written to examples file."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},
            {"name": "read2", "flag": 0, "pos": 100},
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},  # same
            {"name": "read2", "flag": 0, "pos": 500},  # different
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        run_concordance(bam_a, bam_b, out)

        examples = (out / "divergence_examples.tsv").read_text()
        lines = examples.strip().split("\n")

        # Header + 1 divergent read (same_position is not divergent)
        assert len(lines) == 2
        assert "read2" in lines[1]
        assert "same_chrom_diff_pos" in lines[1]

    def test_max_examples_limit(self, temp_dir: Path) -> None:
        """Examples are limited to max_examples."""
        alignments_a = [
            {"name": f"read{i}", "flag": 0, "pos": 100}
            for i in range(50)
        ]
        alignments_b = [
            {"name": f"read{i}", "flag": 0, "pos": 500}  # all divergent
            for i in range(50)
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        run_concordance(bam_a, bam_b, out, max_examples=10)

        examples = (out / "divergence_examples.tsv").read_text()
        lines = examples.strip().split("\n")

        # Header + max 10 examples
        assert len(lines) == 11


class TestSecondarySupplementary:
    """Tests for handling of secondary and supplementary alignments."""

    def test_secondary_ignored(self, temp_dir: Path) -> None:
        """Secondary alignments are ignored in concordance."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},  # primary
            {"name": "read1", "flag": 256, "pos": 200},  # secondary
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},  # primary
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        # Only primary counted
        assert result["counts"]["same_position"] == 1
        assert result["total_reads"] == 1

    def test_supplementary_ignored(self, temp_dir: Path) -> None:
        """Supplementary alignments are ignored in concordance."""
        alignments_a = [
            {"name": "read1", "flag": 0, "pos": 100},  # primary
            {"name": "read1", "flag": 2048, "pos": 300},  # supplementary
        ]
        alignments_b = [
            {"name": "read1", "flag": 0, "pos": 100},  # primary
        ]

        bam_a = create_bam(temp_dir / "a.bam", alignments_a)
        bam_b = create_bam(temp_dir / "b.bam", alignments_b)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["same_position"] == 1
        assert result["total_reads"] == 1


class TestMetadata:
    """Tests for output metadata."""

    def test_names_in_output(self, temp_dir: Path) -> None:
        """Tool names appear in output."""
        alignments = [{"name": "read1", "flag": 0, "pos": 100}]

        bam_a = create_bam(temp_dir / "a.bam", alignments)
        bam_b = create_bam(temp_dir / "b.bam", alignments)
        out = temp_dir / "out"

        result = run_concordance(
            bam_a, bam_b, out, name_a="minimap2", name_b="llmap"
        )

        assert result["name_a"] == "minimap2"
        assert result["name_b"] == "llmap"

    def test_tolerance_in_output(self, temp_dir: Path) -> None:
        """Tolerance appears in output."""
        alignments = [{"name": "read1", "flag": 0, "pos": 100}]

        bam_a = create_bam(temp_dir / "a.bam", alignments)
        bam_b = create_bam(temp_dir / "b.bam", alignments)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out, tolerance_bp=25)

        assert result["tolerance_bp"] == 25


class TestEdgeCases:
    """Tests for edge cases."""

    def test_empty_bams(self, temp_dir: Path) -> None:
        """Empty BAMs produce zero counts.

        Note: total_reads is set to 1 when empty to avoid division by zero
        in percentage calculations, but all counts should be 0.
        """
        bam_a = create_bam(temp_dir / "a.bam", [])
        bam_b = create_bam(temp_dir / "b.bam", [])
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        # total_reads uses `or 1` to avoid div-by-zero
        assert result["total_reads"] == 1
        assert all(v == 0 for v in result["counts"].values())

    def test_single_read_perfect(self, temp_dir: Path) -> None:
        """Single read perfectly concordant."""
        alignments = [{"name": "single", "flag": 0, "pos": 42}]

        bam_a = create_bam(temp_dir / "a.bam", alignments)
        bam_b = create_bam(temp_dir / "b.bam", alignments)
        out = temp_dir / "out"

        result = run_concordance(bam_a, bam_b, out)

        assert result["counts"]["same_position"] == 1
        assert result["percentages"]["same_position"] == 100.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
