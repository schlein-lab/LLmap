#!/usr/bin/env python3
"""Per-read genotyping of HPRC iso-seq reads against the IGHG4-ChimDup T1
catalog entry.

This is the actual silence-test that hap1-vs-hap2 mapping counts do NOT
provide: for each read aligned to hap1, look up the CH1 positions 69 and
252 and assign the read to canonical (ref allele) or duplicate (alt
allele). Counts give the within-haplotype expression ratio of the two
tandem-dup copies.

Reads the T1 entry to get the SNP positions in chromosomal coordinates
(GRCh38), then projects them onto each per-sample hap1 contig via the
sample's GFF annotation (IGHG4 gene → CDS exon 1 → position offset).

This script is a POST-processor for the SAMs produced by hprc_igh_bench.
It does NOT re-run any mapping. Pure parsing.

Usage on the HPC clusterfront (or <hpc-login>):
    python3 per_read_genotype_chimdup.py \\
        --catalog-entry catalog/curated/IGHG4_chimdup_tandem.json \\
        --bench-root /beegfs/u/<user>/llmap_bench/hprc_igh \\
        --manifest /beegfs/u/<user>/llmap_bench/scripts/hprc_igh_manifest.tsv \\
        --out /beegfs/u/<user>/llmap_bench/hprc_igh/per_read_genotype.tsv
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


# ----------------------------------------------------------------------------
# Catalog parsing
# ----------------------------------------------------------------------------

@dataclass
class DiscriminatingSnp:
    """One CH1 SNP from the T1 catalog entry."""
    id: str
    cds_offset: int            # 1-based mRNA-position within IGHG4 CDS
    grch38_pos: int            # 1-based chromosomal position on chr14, GRCh38
    ref: str
    alt: str
    klass: str                 # DUP_fixed / ORIG_polymorph / ...


def load_discriminating_snps(catalog_path: Path) -> list[DiscriminatingSnp]:
    with catalog_path.open() as f:
        entry = json.load(f)
    out = []
    for s in entry["diagnostic_features"]["discriminating_snps"]:
        out.append(DiscriminatingSnp(
            id=s["id"],
            cds_offset=int(s["position"]["cds_offset"]),
            grch38_pos=int(s["position"]["grch38_chr14"]),
            ref=s["ref"], alt=s["alt"],
            klass=s["class"],
        ))
    return out


# ----------------------------------------------------------------------------
# Per-sample GFF projection
# ----------------------------------------------------------------------------

def open_maybe_gz(p: Path):
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return p.open("r")


def find_ighg4_exon1_starts(gff_path: Path) -> list[tuple[str, int, int, str, str]]:
    """Return list of (contig, exon1_start, exon1_end, strand, gene_id)
    for every IGHG4 CH1 / exon-number-1 in this hap's annotation.

    HPRC R2 GFFs (CAT-source) use:
      - feature type 'transcript' (not 'mRNA')
      - attribute key 'gene_name=IGHG4'
      - attribute 'exon_number=1' on the CH1 exon record
      - strand '-' for IGHG4

    A single hap may have multiple IGHG4 gene records (tandem duplicate
    → e.g. HG00329 hap1 has G030290 + G030293). We return one entry per
    gene_id, preferring the IGHG4-201 / IGHG4-202 transcript whose
    exon_number=1 record matches the 295 bp CH1 length pattern.
    """
    ighg4_gene_ids: set[str] = set()
    transcript_to_gene: dict[str, str] = {}
    # per gene_id, collect exon-number=1 records
    by_gene: dict[str, list[tuple[str, int, int, str]]] = {}

    with open_maybe_gz(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)

            if ftype == "gene":
                if attrs.get("gene_name") == "IGHG4":
                    ighg4_gene_ids.add(attrs.get("ID", ""))

            elif ftype == "transcript":
                if attrs.get("gene_name") == "IGHG4":
                    parent = attrs.get("Parent", "")
                    transcript_to_gene[attrs.get("ID", "")] = parent

            elif ftype == "exon":
                if attrs.get("gene_name") != "IGHG4":
                    continue
                if attrs.get("exon_number") != "1":
                    continue
                # walk through possibly-multiple parents (mRNAs) to find the gene
                for parent in attrs.get("Parent", "").split(","):
                    gene_id = transcript_to_gene.get(parent)
                    if not gene_id or gene_id not in ighg4_gene_ids:
                        continue
                    by_gene.setdefault(gene_id, []).append(
                        (f[0], int(f[3]), int(f[4]), f[6]))

    results: list[tuple[str, int, int, str, str]] = []
    for gene_id, exons in by_gene.items():
        # Prefer the exon matching the canonical 295 bp CH1 length; else first.
        chosen = next((e for e in exons if (e[2] - e[1] + 1) == 295), exons[0])
        results.append((*chosen, gene_id))
    return results


# ----------------------------------------------------------------------------
# SAM read genotyping
# ----------------------------------------------------------------------------

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def base_at_ref_pos(read_pos: int, cigar: str, seq: str, query_pos_target: int) -> str | None:
    """Walk CIGAR to find the read-base aligned to the given reference position
    (1-based, relative to mapping start).

    read_pos = SAM POS (1-based on the contig)
    query_pos_target = the contig position we want to genotype
    Returns the read base, or None if the position is deleted/clipped.
    """
    ref = read_pos
    qry = 0
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op in ("M", "=", "X"):
            if ref <= query_pos_target < ref + length:
                offset = query_pos_target - ref
                if qry + offset < len(seq):
                    return seq[qry + offset]
                return None
            ref += length
            qry += length
        elif op in ("I", "S"):
            qry += length
        elif op in ("D", "N"):
            if ref <= query_pos_target < ref + length:
                return None  # deleted in read
            ref += length
        # H, P don't advance either
    return None


def read_sam_extract_offset(sam_path: Path, expected_contig: str) -> tuple[str, int] | None:
    """Read @SQ headers, find the rname matching the gff contig (allowing for
    the LLmap extract suffix ':start-end'), and return (rname_in_sam, offset_bp).

    HPRC extract step emits SAM with RNAME = '<contig>:<start>-<end>' where
    <start> is 1-based.  Position in alignments is 1-based relative to that
    sub-region, so to convert a GFF (orig-contig) coord to SAM-space we
    subtract `<start> - 1`.
    """
    if not sam_path.exists():
        return None
    with sam_path.open() as f:
        for line in f:
            if not line.startswith("@"):
                break
            if not line.startswith("@SQ"):
                continue
            sn = next((t.split(":", 1)[1] for t in line.split("\t")
                       if t.startswith("SN:")), None)
            if not sn:
                continue
            # split on the LAST ':' so embedded ':' in contig names stays intact
            if ":" in sn:
                base, region = sn.rsplit(":", 1)
                if base == expected_contig and "-" in region:
                    try:
                        start = int(region.split("-", 1)[0])
                        return sn, start - 1
                    except ValueError:
                        continue
            if sn == expected_contig:
                return sn, 0
    return None


def iter_sam(sam_path: Path, sam_rname: str) -> Iterator[dict]:
    """Yield mapped primary alignments for `sam_rname` (already resolved)."""
    if not sam_path.exists():
        return
    with sam_path.open() as f:
        for line in f:
            if line.startswith("@"):
                continue
            f0 = line.rstrip("\n").split("\t")
            if len(f0) < 11:
                continue
            flag = int(f0[1])
            if flag & 0x4:        # unmapped
                continue
            if flag & 0x100:      # secondary
                continue
            if flag & 0x800:      # supplementary
                continue
            if f0[2] != sam_rname:
                continue
            yield {
                "qname": f0[0], "flag": flag, "rname": f0[2],
                "pos": int(f0[3]), "mapq": int(f0[4]),
                "cigar": f0[5], "seq": f0[9],
            }


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def genotype_sample(sample: str, bench_root: Path, snps: list[DiscriminatingSnp]) -> list[dict]:
    """For one HPRC sample, walk hap1 + hap2 SAMs and tally per-SNP read
    base distributions at the projected positions."""
    sample_dir = bench_root / "samples" / sample
    if not sample_dir.is_dir():
        return []

    rows: list[dict] = []
    for hap in (1, 2):
        sam = sample_dir / f"llmap.hap{hap}.sam"
        coords_tsv = sample_dir / f"hap{hap}.igh.coords.tsv"
        if not sam.exists() or not coords_tsv.exists():
            continue

        # The extract step writes one row per IGH gene with contig/start/end/strand.
        # We need IGHG4 CDS-exon-1 starts; the GFF holds those.
        # GFF location: standard HPRC R2 layout
        #   assemblies_r2/<SAMPLE>/annotation/<SAMPLE>_hap<N>.gff3.gz
        # Fall back: parse from extract.log if standard path doesn't exist.
        gff_path = Path(
            f"/beegfs/u/<user>/<group>/shared/references/human_pangenome/"
            f"assemblies_r2/{sample}/annotation/{sample}_hap{hap}.gff3.gz"
        )
        if not gff_path.exists():
            log = sample_dir / "extract.log"
            if log.exists():
                # walk log lines, find the one matching this hap
                in_hap = False
                resolved = None
                for line in log.read_text().splitlines():
                    if line.startswith(f"[asm] hap={hap}"):
                        in_hap = True
                        continue
                    if line.startswith("[asm] hap="):
                        in_hap = False
                    if in_hap and line.startswith("[asm] gff="):
                        resolved = line.split("=", 1)[1].strip()
                        break
                if resolved and Path(resolved).exists():
                    gff_path = Path(resolved)
                else:
                    continue
            else:
                continue

        exon1_starts = find_ighg4_exon1_starts(Path(gff_path))
        if not exon1_starts:
            continue

        # Per IGHG4 copy on this hap (gene_id distinguishes tandem duplicates):
        for copy_idx, (contig, start, end, strand, gene_id) in enumerate(exon1_starts):
            # Resolve SAM RNAME + offset for this contig (LLmap extract uses
            # '<contig>:<extract_start>-<extract_end>' as RNAME; positions are
            # 1-based relative to extract_start).
            res = read_sam_extract_offset(sam, contig)
            if res is None:
                continue
            sam_rname, offset_bp = res

            for snp in snps:
                # CH1 starts at the 5' end of the CH1 exon on the gene strand
                if strand == "+":
                    target_pos_orig = start + (snp.cds_offset - 1)
                else:
                    target_pos_orig = end - (snp.cds_offset - 1)
                # Map to SAM-relative 1-based coord
                target_pos = target_pos_orig - offset_bp

                # Tally bases across all primary mapped reads at target_pos
                base_counts: dict[str, int] = {}
                for aln in iter_sam(sam, sam_rname):
                    b = base_at_ref_pos(aln["pos"], aln["cigar"],
                                         aln["seq"], target_pos)
                    if b:
                        base_counts[b.upper()] = base_counts.get(b.upper(), 0) + 1

                # Classify the call
                total = sum(base_counts.values())
                ref_n = base_counts.get(snp.ref, 0)
                alt_n = base_counts.get(snp.alt, 0)
                if total == 0:
                    call = "no_coverage"
                elif ref_n > alt_n * 3:
                    call = "canonical"
                elif alt_n > ref_n * 3:
                    call = "dup_fixed_pattern" if snp.klass == "DUP_fixed" else "orig_polymorph_alt"
                else:
                    call = "mixed"

                rows.append({
                    "sample": sample, "hap": hap,
                    "ighg4_copy_idx": copy_idx,
                    "gene_id": gene_id,
                    "contig": contig,
                    "snp_id": snp.id, "snp_class": snp.klass,
                    "target_pos": target_pos,
                    "ref": snp.ref, "alt": snp.alt,
                    "ref_count": ref_n, "alt_count": alt_n, "total": total,
                    "base_counts": json.dumps(base_counts, sort_keys=True),
                    "call": call,
                })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog-entry", required=True, type=Path,
                    help="path to catalog/curated/IGHG4_chimdup_tandem.json")
    ap.add_argument("--bench-root", required=True, type=Path,
                    help="root of HPRC IGH bench (contains samples/<S>/)")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="hprc_igh_manifest.tsv")
    ap.add_argument("--out", required=True, type=Path,
                    help="output TSV path")
    args = ap.parse_args()

    snps = load_discriminating_snps(args.catalog_entry)
    print(f"[per_read_genotype] {len(snps)} discriminating SNPs loaded "
          f"from {args.catalog_entry.name}")

    samples: list[str] = []
    with args.manifest.open() as f:
        next(f)  # header
        for line in f:
            f0 = line.rstrip("\n").split("\t")
            if len(f0) >= 2:
                samples.append(f0[1])
    print(f"[per_read_genotype] {len(samples)} samples in manifest")

    with args.out.open("w") as out:
        header = ["sample", "hap", "ighg4_copy_idx", "gene_id", "contig",
                  "snp_id", "snp_class", "target_pos",
                  "ref", "alt", "ref_count", "alt_count", "total",
                  "base_counts", "call"]
        out.write("\t".join(header) + "\n")

        n_with_data = 0
        for s in samples:
            rows = genotype_sample(s, args.bench_root, snps)
            if rows:
                n_with_data += 1
            for r in rows:
                out.write("\t".join(str(r.get(h, "")) for h in header) + "\n")

    print(f"[per_read_genotype] {n_with_data}/{len(samples)} samples produced rows → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
