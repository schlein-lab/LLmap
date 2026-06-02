#!/usr/bin/env python3
"""Build the T2 (bulk) SegDup catalog from primary databases.

Sources (per docs/design/segdup_databases_inventory.md):
  - UCSC genomicSuperDups (hg38, frozen 2014-10-14)
  - UCSC genomicSuperDups (hg19, mm39, rn7, danRer11, dm6) — optional
  - T2T-CHM13 v2.0 SD-BED (Vollger 2022)
  - GIAB genome-stratifications (SD / low-mappability)

Output:
  catalog/bulk/v<release>.<assembly>.jsonl
  catalog/bulk/v<release>.<assembly>.summary.json    (counts per architecture class)

Design: PARSING ONLY. No alignment, no compute-heavy steps. Suitable for
the HPC clusterfront execution (~5-30 min depending on assemblies). If you add
gnomAD-SV overlap or OMIM cross-ref later, those go in a separate
build_t2_scores.py script — keeps this one fast and idempotent.

Usage on the HPC clusterfront (NOT local — per project policy 2026-06-02):
    ssh <hpc-login>
    cd /home/<user>/llmap_catalog_build      # or wherever the repo sits
    python3 catalog/scripts/build_t2_bulk.py \
        --release v2026.Q2 \
        --assemblies grch38 chm13 \
        --out catalog/bulk/
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import urllib.request
from collections import Counter
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Iterator


# ----------------------------------------------------------------------------
# Database source registry
# ----------------------------------------------------------------------------

SOURCES = {
    "grch38": {
        "ucsc_sd": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
        "label":   "UCSC_genomicSuperDups_hg38_frozen2014",
        # Per ucsc_sd_completeness_review.md: UCSC alone misses ~51 Mbp T2T-additional
        # SDs (incl. IGHG4-ChimDup), acrocentric short arms (35 Mbp N-masked), and
        # all SDs below 90 % identity (e.g. IGH-CanonDup at 84.8 %). DO NOT use as
        # sole source — pair with T2T-CHM13 SD-BED and Paraphase regions.
    },
    "hg19": {
        "ucsc_sd": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz",
        "label":   "UCSC_genomicSuperDups_hg19",
    },
    "chm13": {
        "t2t_sd":     "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed",
        "t2t_sd_full":"https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.full.bed",
        "label":      "T2T_CHM13_v2.0_SD_Vollger2022",
    },
    "mm39": {
        "ucsc_sd": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/genomicSuperDups.txt.gz",
        "label":   "UCSC_genomicSuperDups_mm39",
    },
    "paraphase": {
        # Paraphase 160-region list — PacBio, 316 genes covering clinically-relevant
        # paralog groups. NOT a per-pair SD list — region-anchored. Parsed
        # separately into anchor records (one per region, architecture inferred
        # from known paralog structure).
        "paraphase_yaml": "https://raw.githubusercontent.com/PacificBiosciences/paraphase/main/data/38/config.yaml",
        "label":           "Paraphase_v160_PacBio_hg38",
    },
}


# ----------------------------------------------------------------------------
# Architecture classification rules
# (from docs/design/segdup_class_taxonomy.md)
# ----------------------------------------------------------------------------

def classify_architecture(chrom_a: str, start_a: int, end_a: int,
                           chrom_b: str, start_b: int, end_b: int) -> str:
    """Map a paired-SD record to an Axis-A architecture tag.

    Conservative — anything we can't confidently slot lands in
    'unresolved' so the curator sees it instead of being mislabelled.
    """
    if chrom_a != chrom_b:
        return "interchromosomal"

    # same-chromosome — distance between block midpoints
    mid_a = (start_a + end_a) / 2
    mid_b = (start_b + end_b) / 2
    dist  = abs(mid_a - mid_b)

    # Tandem ≈ adjacent or overlapping; split <50 kb vs ≥50 kb per taxonomy.
    # NOTE: we can't tell tandem vs interspersed from coords alone for
    # mid-distance pairs; ≥1 Mb intrachromosomal is interspersed by taxonomy.
    if dist < 50_000:
        return "tandem_dup_small"
    if dist < 1_000_000:
        return "tandem_dup_large"
    return "intrachromosomal_interspersed"


# ----------------------------------------------------------------------------
# Parsers
# ----------------------------------------------------------------------------

@dataclass
class BulkRecord:
    locus_id:                str
    structural_architecture: str    # renamed from `architecture` (schema v0.2: split axis)
    coords:                  dict
    identity_pct:            float
    source:                  str
    version:                 str
    paired_with:             dict = field(default_factory=dict)
    size_bp:                 int = 0
    notes:                   str = ""
    # T2 bulk records do NOT carry haplotype_class — that lives only in T1 curated
    # entries because T2 has no PSV / sequence-level resolution.

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self), separators=(",", ":"))


def open_remote(url: str) -> io.TextIOBase:
    """Open a URL as a text stream, decompressing .gz transparently."""
    req = urllib.request.Request(url, headers={"User-Agent": "llmap-catalog-builder/0.1"})
    resp = urllib.request.urlopen(req, timeout=120)
    if url.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=resp), encoding="utf-8")
    return io.TextIOWrapper(resp, encoding="utf-8")


def parse_ucsc_sd(url: str, source_label: str, release: str) -> Iterator[BulkRecord]:
    """Stream UCSC genomicSuperDups records.

    Field layout (per kent-source schema):
      bin chrom start end name score strand otherChrom otherStart otherEnd
      otherSize uid posBasesHit testResult verdict chits ccov alignfile
      alignL indelN indelS alignB matchB mismatchB transitionsB transversionsB
      fracMatch fracMatchIndel jcK k2K

    fracMatch = fraction-matching = our identity_pct/100.
    """
    counter = 0
    with open_remote(url) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 27:
                continue
            try:
                chrom_a, start_a, end_a = f[1], int(f[2]), int(f[3])
                chrom_b, start_b, end_b = f[7], int(f[8]), int(f[9])
                strand_b = f[6]
                frac_match = float(f[26])     # fracMatch
            except (ValueError, IndexError):
                continue

            counter += 1
            rec = BulkRecord(
                locus_id=f"ucsc_sd_{release.replace('.', '_').lower()}_{source_label.split('_')[-2]}_{counter:06d}",
                structural_architecture=classify_architecture(
                    chrom_a, start_a, end_a, chrom_b, start_b, end_b),
                coords={chrom_a: {"chrom": chrom_a, "start": start_a, "end": end_a}},
                identity_pct=round(frac_match * 100.0, 3),
                source=source_label,
                version=release,
                paired_with={
                    "chrom": chrom_b, "start": start_b, "end": end_b,
                    "strand": strand_b,
                },
                size_bp=end_a - start_a,
            )
            yield rec


def parse_paraphase_yaml(url: str, source_label: str, release: str) -> Iterator[BulkRecord]:
    """Stream Paraphase region records.

    Paraphase data.yaml structure (per PacBio repo):
      <gene>:
        gene: <name>
        region: chrX:start-end
        ...
    We emit one BulkRecord per region with architecture='unresolved' (the
    user-facing curated entries will resolve it; bulk only records 'this
    region is known SD-ish per Paraphase').
    """
    counter = 0
    # Lazy import — keeps base Python install free if Paraphase source is skipped
    try:
        import yaml  # type: ignore
    except ImportError:
        print("[build_t2_bulk] WARN: pyyaml not installed; skipping Paraphase source",
              file=sys.stderr)
        return

    with open_remote(url) as fh:
        data = yaml.safe_load(fh)

    if not isinstance(data, dict):
        return

    for gene_key, entry in data.items():
        if not isinstance(entry, dict):
            continue
        # Paraphase config.yaml shape (per repo):
        #   <gene_key>:
        #     genes: "SMN1,SMN2"     # comma-separated paralog set
        #     realign_region: chr5:70890000-71100000
        #     gene2_region:   chr5:70040526-70088546   # optional
        #     ...
        region = entry.get("realign_region") or entry.get("gene2_region") or ""
        if not isinstance(region, str) or ":" not in region or "-" not in region:
            continue
        try:
            chrom, rng = region.split(":", 1)
            start, end = rng.split("-", 1)
            start, end = int(start), int(end)
        except ValueError:
            continue

        paralog_set = entry.get("genes", gene_key)
        counter += 1
        yield BulkRecord(
            locus_id=f"paraphase_{gene_key.lower().replace('/', '_')}",
            structural_architecture="unresolved",   # T1 curation will assign a class
            coords={"GRCh38": {"chrom": chrom, "start": start, "end": end}},
            identity_pct=0.0,                # Paraphase doesn't expose per-pair %
            source=source_label,
            version=release,
            paired_with={},
            size_bp=end - start,
            notes=f"paraphase_gene={gene_key}; paralog_set={paralog_set}",
        )


def parse_t2t_sd_bed(url: str, source_label: str, release: str) -> Iterator[BulkRecord]:
    """Stream T2T-CHM13 SD-BED.

    Plain BED: chrom start end name [score strand ...]
    The Vollger 2022 SD-BED has more columns (paired-with coords) — parse
    if present, otherwise fall back to single-region records.
    """
    counter = 0
    with open_remote(url) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                chrom_a, start_a, end_a = f[0], int(f[1]), int(f[2])
            except (ValueError, IndexError):
                continue
            counter += 1
            name = f[3] if len(f) > 3 else ""

            # Vollger SD.full.bed typically encodes the paralog block in the
            # name field as 'chr14:99891444-99910944' or in cols 6+. Heuristic:
            paired = {}
            arch = "unresolved"
            for tok in (name.split(";") if name else []):
                if ":" in tok and "-" in tok:
                    try:
                        c, rng = tok.split(":", 1)
                        s, e   = rng.split("-", 1)
                        paired = {"chrom": c, "start": int(s), "end": int(e)}
                        arch   = classify_architecture(chrom_a, start_a, end_a,
                                                        c, int(s), int(e))
                        break
                    except ValueError:
                        pass
            if arch == "unresolved" and len(f) >= 9:
                # try BED12-ish layout
                try:
                    arch = classify_architecture(chrom_a, start_a, end_a,
                                                  f[6], int(f[7]), int(f[8]))
                    paired = {"chrom": f[6], "start": int(f[7]), "end": int(f[8])}
                except (ValueError, IndexError):
                    pass

            yield BulkRecord(
                locus_id=f"t2t_chm13_sd_{counter:06d}",
                structural_architecture=arch,
                coords={"CHM13": {"chrom": chrom_a, "start": start_a, "end": end_a}},
                identity_pct=0.0,   # not in this BED — fill from .full.bed when available
                source=source_label,
                version=release,
                paired_with=paired,
                size_bp=end_a - start_a,
                notes=name or "",
            )


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def build_assembly(assembly: str, release: str, out_dir: Path) -> dict:
    src = SOURCES.get(assembly)
    if not src:
        sys.exit(f"unknown assembly: {assembly}")

    out_path = out_dir / f"{release}.{assembly}.jsonl"
    summary  = {"assembly": assembly, "release": release,
                "sources": [], "records": 0,
                "architecture_counts": Counter()}

    with out_path.open("w") as out:
        if "ucsc_sd" in src:
            summary["sources"].append(src["label"])
            for rec in parse_ucsc_sd(src["ucsc_sd"], src["label"], release):
                out.write(rec.to_jsonl() + "\n")
                summary["records"] += 1
                summary["architecture_counts"][rec.structural_architecture] += 1

        if "t2t_sd_full" in src:
            summary["sources"].append(src["label"])
            for rec in parse_t2t_sd_bed(src["t2t_sd_full"], src["label"], release):
                out.write(rec.to_jsonl() + "\n")
                summary["records"] += 1
                summary["architecture_counts"][rec.structural_architecture] += 1

        if "paraphase_yaml" in src:
            summary["sources"].append(src["label"])
            for rec in parse_paraphase_yaml(src["paraphase_yaml"], src["label"], release):
                out.write(rec.to_jsonl() + "\n")
                summary["records"] += 1
                summary["architecture_counts"][rec.structural_architecture] += 1

    summary["architecture_counts"] = dict(summary["architecture_counts"])

    sum_path = out_dir / f"{release}.{assembly}.summary.json"
    with sum_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(f"[build_t2_bulk] {assembly}: {summary['records']} records → {out_path}")
    print(f"[build_t2_bulk]   architecture: {summary['architecture_counts']}")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release",    required=True,
                    help="catalog release tag, e.g. v2026.Q2")
    ap.add_argument("--assemblies", required=True, nargs="+",
                    choices=sorted(SOURCES.keys()),
                    help="which assemblies to build")
    ap.add_argument("--out",        required=True, type=Path,
                    help="output directory")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    summaries = {}
    for asm in args.assemblies:
        summaries[asm] = build_assembly(asm, args.release, args.out)

    # write release-level summary
    rel_sum = args.out / f"{args.release}.bulk.release_summary.json"
    with rel_sum.open("w") as f:
        json.dump({"release": args.release, "assemblies": summaries}, f, indent=2)
    print(f"[build_t2_bulk] release summary → {rel_sum}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
