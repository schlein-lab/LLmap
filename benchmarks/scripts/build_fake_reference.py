#!/usr/bin/env python3
"""Stitch a fake reference FASTA from an organism's specific_loci JSONs.

This is fine for benchmarking the *mapper*: we just need self-consistent
coordinates between the reads we draw and the reference we map against.
The sequence content does not have to match any real genome.

Per locus:
  - read the coordinates entry for the requested assembly_key (falls back
    to the first available key if missing)
  - allocate a contig named after that chr (concatenating loci that share
    a chr onto the same contig, leaving gaps between blocks)
  - fill the locus block with deterministically-seeded random DNA
  - fill inter-locus gaps with random DNA of length min(10kb, gap)

Output:
  <output>/reference.fa
  <output>/reference.fa.fai-style.tsv  (informal, name/length/loci_count)
  <output>/loci_index.tsv               (name, chr, start, end, source_path)

CLI:
  python build_fake_reference.py --organism mouse \
      --registry benchmarks/species_registry.json \
      --root /home/<user>/llmap-local \
      --output benchmarks/datasets/cache/mouse_mm39
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

ALPHABET = "ACGT"
PAD_BETWEEN_LOCI = 10_000  # generated bp between blocks on same contig
MIN_CONTIG_LEN = 1_000     # avoid pathologically tiny contigs


def iter_loci_jsons(loci_dir: Path) -> Iterable[Path]:
    if not loci_dir.exists():
        return []
    return sorted(p for p in loci_dir.rglob("*.json") if p.is_file())


def pick_coords(record: dict, assembly_key: str) -> Optional[Tuple[str, int, int, str]]:
    coords = record.get("coordinates")
    if not isinstance(coords, dict) or not coords:
        return None
    chosen = coords.get(assembly_key)
    used_key = assembly_key
    if chosen is None:
        # fall back to first valid entry
        for k, v in coords.items():
            if isinstance(v, dict) and "chr" in v and "start" in v and "end" in v:
                chosen = v
                used_key = k
                break
    if not isinstance(chosen, dict):
        return None
    chrom = str(chosen.get("chr") or "")
    try:
        start = int(chosen.get("start", 0))
        end = int(chosen.get("end", 0))
    except (TypeError, ValueError):
        return None
    if not chrom or end <= start:
        return None
    return chrom, start, end, used_key


def random_seq(rng: random.Random, n: int) -> str:
    # bytearray + random.choices is reasonably fast and stable across versions
    return "".join(rng.choices(ALPHABET, k=n))


def write_fasta(path: Path, contigs: Dict[str, str], line_width: int = 80) -> None:
    with path.open("w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width])
                fh.write("\n")


def collect_loci(
    loci_dir: Path,
    assembly_key: str,
    sub_organisms: Optional[List[str]] = None,
    registry_root: Optional[Path] = None,
) -> List[Tuple[str, int, int, str, Path]]:
    """Return list of (chr, start, end, locus_name, source_path) tuples.

    For organisms with sub_organisms (great_apes), recurse into each subdir.
    """
    out: List[Tuple[str, int, int, str, Path]] = []
    if sub_organisms and registry_root is not None:
        for sub in sub_organisms:
            sub_dir = registry_root / loci_dir / sub / "specific_loci"
            if not sub_dir.exists():
                sub_dir = registry_root / loci_dir / sub
            for jp in iter_loci_jsons(sub_dir):
                _maybe_add(out, jp, assembly_key, prefix=sub)
    else:
        for jp in iter_loci_jsons(loci_dir):
            _maybe_add(out, jp, assembly_key)
    return out


def _maybe_add(
    out: List[Tuple[str, int, int, str, Path]],
    jp: Path,
    assembly_key: str,
    prefix: Optional[str] = None,
) -> None:
    try:
        rec = json.loads(jp.read_text())
    except (OSError, json.JSONDecodeError):
        return
    pc = pick_coords(rec, assembly_key)
    if pc is None:
        return
    chrom, start, end, _used = pc
    name = rec.get("name") or jp.stem
    if prefix:
        chrom = f"{prefix}_{chrom}"
        name = f"{prefix}__{name}"
    out.append((chrom, start, end, name, jp))


def build_contigs(
    loci: List[Tuple[str, int, int, str, Path]],
    seed: int,
) -> Tuple[Dict[str, str], List[Tuple[str, str, int, int, str]]]:
    """Group loci by chr, place onto a single contig per chr.

    Returns (contigs, index) where index entries are
    (locus_name, chr, mapped_start, mapped_end, source_path_str).
    """
    by_chr: Dict[str, List[Tuple[int, int, str, Path]]] = defaultdict(list)
    for chrom, start, end, name, src in loci:
        by_chr[chrom].append((start, end, name, src))

    contigs: Dict[str, str] = {}
    index: List[Tuple[str, str, int, int, str]] = []

    for chrom, items in sorted(by_chr.items()):
        items.sort(key=lambda x: (x[0], x[1]))
        rng = random.Random(f"{seed}:{chrom}")
        cursor = 0
        buf: List[str] = []
        last_orig_end: Optional[int] = None
        for orig_start, orig_end, name, src in items:
            block_len = max(1, orig_end - orig_start)
            # leading pad
            if last_orig_end is None:
                pad = min(PAD_BETWEEN_LOCI, max(1000, block_len // 4))
            else:
                # honor original gap, clamped
                gap = max(0, orig_start - last_orig_end)
                pad = min(PAD_BETWEEN_LOCI, max(500, gap if gap > 0 else 500))
            buf.append(random_seq(rng, pad))
            cursor += pad
            # locus block
            buf.append(random_seq(rng, block_len))
            mapped_start = cursor
            mapped_end = cursor + block_len
            cursor = mapped_end
            index.append((name, chrom, mapped_start, mapped_end, str(src)))
            last_orig_end = orig_end
        # trailing pad to round up
        if cursor < MIN_CONTIG_LEN:
            buf.append(random_seq(rng, MIN_CONTIG_LEN - cursor))
            cursor = MIN_CONTIG_LEN
        else:
            buf.append(random_seq(rng, 2000))
            cursor += 2000
        contigs[chrom] = "".join(buf)
    return contigs, index


def main() -> int:
    ap = argparse.ArgumentParser(description="Build fake reference from specific_loci JSONs")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--registry", required=True, help="species_registry.json path")
    ap.add_argument("--root", required=True, help="llmap-local repo root (paths in registry are resolved against this)")
    ap.add_argument("--output", required=True, help="output directory (created if missing)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    registry = json.loads(Path(args.registry).read_text())
    if args.organism not in registry:
        sys.stderr.write(f"organism {args.organism!r} not in registry\n")
        return 2
    cfg = registry[args.organism]

    loci_dir = root / cfg["loci_dir"]
    assembly_key = cfg["assembly_key"]
    sub_organisms = cfg.get("sub_organisms")

    if not loci_dir.exists():
        sys.stderr.write(f"loci_dir does not exist: {loci_dir}\n")
        return 3

    loci = collect_loci(loci_dir, assembly_key, sub_organisms=sub_organisms, registry_root=root)
    if not loci:
        sys.stderr.write(
            f"no usable specific_loci JSONs found for {args.organism} "
            f"(loci_dir={loci_dir}, assembly_key={assembly_key})\n"
        )
        return 4

    contigs, index = build_contigs(loci, seed=args.seed)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    fa_path = out_dir / "reference.fa"
    write_fasta(fa_path, contigs)

    # small index files
    with (out_dir / "loci_index.tsv").open("w") as fh:
        fh.write("locus_name\tcontig\tmapped_start\tmapped_end\tsource_path\n")
        for row in index:
            fh.write("\t".join(str(c) for c in row) + "\n")
    with (out_dir / "reference.contigs.tsv").open("w") as fh:
        fh.write("contig\tlength\tloci_count\n")
        per_chr_count: Dict[str, int] = defaultdict(int)
        for _, contig, _, _, _ in index:
            per_chr_count[contig] += 1
        for name, seq in contigs.items():
            fh.write(f"{name}\t{len(seq)}\t{per_chr_count.get(name, 0)}\n")

    print(
        json.dumps(
            {
                "organism": args.organism,
                "assembly_key": assembly_key,
                "loci_used": len(index),
                "contigs": len(contigs),
                "total_bp": sum(len(s) for s in contigs.values()),
                "reference_fa": str(fa_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
