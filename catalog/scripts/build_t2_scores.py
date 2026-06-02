#!/usr/bin/env python3
"""Score T2 bulk records and produce a T1-promotion candidate list.

Promotion-score rubric (per docs/design/segdup_catalog_spec.md §1.1):
  - OMIM / ClinGen disease gene overlap:        +5
  - Antibody / immune-receptor locus:           auto-T1 (no threshold)
  - gnomAD-SV allele frequency > 1 % anywhere:  +2
  - HPRC pangenome haplotype-variability flag:  +2
  - GIAB low-mappability flag:                  +1

Threshold for promotion candidacy: score ≥ 5 (manual review required).

This is the SECOND-pass enrichment step. The first-pass T2 records from
build_t2_bulk.py carry only coordinates + architecture. This script
joins them with auxiliary annotations (downloaded once, cached locally)
and emits:
  - catalog/bulk/v<release>.<assembly>.scored.jsonl  (T2 + score field)
  - catalog/bulk/v<release>.promotion_candidates.tsv  (score ≥ 5, sorted)

Sources of enrichment (each optional — skipped with warning if unavailable):
  - GENCODE gene BED → OMIM / ClinGen overlap
  - gnomAD-SV v4 sites VCF → AF
  - HPRC R2 polymorphic-SD list (Jeong 2025 suppl, if available)
  - GIAB stratifications BED (low-mappability subset)

Designed to run on <hpc-login>. Pure parsing/joining, no compute.

Usage:
    python3 build_t2_scores.py \\
        --release v2026.Q2 \\
        --in catalog/bulk/v2026.Q2.grch38.jsonl \\
        --assembly grch38 \\
        --omim-clingen-bed /beegfs/u/<user>/.../clingen_disease_genes.bed \\
        --gnomad-sv /beegfs/u/<user>/.../gnomad.v4.1.sv.sites.vcf.bgz \\
        --giab-stratifications /beegfs/u/<user>/.../GRCh38_lowmappabilityall.bed.gz \\
        --out-scored catalog/bulk/v2026.Q2.grch38.scored.jsonl \\
        --out-candidates catalog/bulk/v2026.Q2.promotion_candidates.tsv
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path


# ----------------------------------------------------------------------------
# Heuristic gene-name → category mapping
# (T1 auto-promote loci by gene-symbol regex match — for the auto-T1 rule)
# ----------------------------------------------------------------------------

# Immune-receptor / antibody locus prefixes (per IMGT / GENCODE biotype IG_*_gene)
AUTO_T1_GENE_PREFIXES = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

# MHC genes (block-level promotion target)
MHC_GENE_PREFIXES = ("HLA-", "MICA", "MICB", "TAP1", "TAP2", "HLA")


def gene_categories(gene_set: str) -> list[str]:
    """Return promotion categories triggered by a gene-name string."""
    cats = []
    upper = gene_set.upper()
    for token in upper.replace(",", " ").split():
        for pfx in AUTO_T1_GENE_PREFIXES:
            if token.startswith(pfx):
                cats.append("immune_receptor_auto_T1")
                break
        for pfx in MHC_GENE_PREFIXES:
            if token.startswith(pfx):
                cats.append("mhc_block")
                break
    return list(dict.fromkeys(cats))   # dedupe, preserve order


# ----------------------------------------------------------------------------
# Interval index for fast BED overlap
# ----------------------------------------------------------------------------

class IntervalIndex:
    """Per-chromosome sorted-list overlap test. Good enough for ~10⁵ intervals.

    For ~10⁷ T2 records × 10⁵ BED intervals this is O(N log M) per query.
    """

    def __init__(self) -> None:
        self.by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        self._sorted = False

    def add(self, chrom: str, start: int, end: int, label: str = "") -> None:
        self.by_chrom[chrom].append((start, end, label))

    def finalize(self) -> None:
        for k in self.by_chrom:
            self.by_chrom[k].sort()
        self._sorted = True

    def overlaps(self, chrom: str, start: int, end: int) -> list[str]:
        if not self._sorted:
            self.finalize()
        lst = self.by_chrom.get(chrom, [])
        hits = []
        # linear scan from first match — fine for genome-scale unless this is
        # called millions of times. Could be replaced with bisect + interval-tree.
        for s, e, lab in lst:
            if e < start:
                continue
            if s > end:
                break
            hits.append(lab)
        return hits


def load_bed(path: Path) -> IntervalIndex:
    idx = IntervalIndex()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                idx.add(f[0], int(f[1]), int(f[2]),
                        f[3] if len(f) > 3 else "")
            except ValueError:
                continue
    idx.finalize()
    return idx


# ----------------------------------------------------------------------------
# Scoring
# ----------------------------------------------------------------------------

def score_record(rec: dict,
                  assembly: str,
                  gene_bed: IntervalIndex | None,
                  gnomad_index: IntervalIndex | None,
                  giab_index: IntervalIndex | None) -> tuple[int, list[str], list[str]]:
    """Return (score, evidence_tags, gene_overlaps) for one T2 record."""
    score = 0
    tags: list[str] = []
    genes: list[str] = []

    coords = rec.get("coords", {})
    asm_key = "GRCh38" if assembly == "grch38" else "CHM13" if assembly == "chm13" else assembly
    if asm_key not in coords or coords[asm_key] is None:
        return 0, [], []
    c = coords[asm_key]
    chrom, start, end = c["chrom"], c["start"], c["end"]

    # Gene-symbol based auto-T1 (works without external BEDs)
    notes = rec.get("notes", "")
    if "paraphase_gene=" in notes or "paralog_set=" in notes:
        cats = gene_categories(notes)
        for cat in cats:
            tags.append(f"auto_T1:{cat}")
            score += 100   # well above threshold — paraphase already curated
        if cats:
            genes.append(notes)

    if gene_bed:
        hits = gene_bed.overlaps(chrom, start, end)
        if hits:
            score += 5
            tags.append("omim_clingen_overlap")
            genes.extend(hits)
            # extra auto-T1 for antibody/MHC overlaps
            for label in hits:
                cats = gene_categories(label)
                for cat in cats:
                    if f"auto_T1:{cat}" not in tags:
                        tags.append(f"auto_T1:{cat}")
                        score += 100

    if gnomad_index:
        if gnomad_index.overlaps(chrom, start, end):
            score += 2
            tags.append("gnomad_sv_overlap")

    if giab_index:
        if giab_index.overlaps(chrom, start, end):
            score += 1
            tags.append("giab_low_mappability")

    return score, tags, genes


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release",   required=True)
    ap.add_argument("--in",        dest="in_path", required=True, type=Path,
                    help="bulk JSONL from build_t2_bulk.py")
    ap.add_argument("--assembly",  required=True,
                    choices=["grch38", "chm13", "paraphase"])
    ap.add_argument("--omim-clingen-bed", type=Path, default=None,
                    help="BED of disease-gene loci (col 4 = gene name); optional")
    ap.add_argument("--gnomad-sv",        type=Path, default=None,
                    help="gnomAD-SV sites BED/VCF; optional")
    ap.add_argument("--giab-stratifications", type=Path, default=None,
                    help="GIAB low-mappability BED; optional")
    ap.add_argument("--out-scored",        required=True, type=Path)
    ap.add_argument("--out-candidates",    required=True, type=Path)
    ap.add_argument("--threshold",         type=int, default=5)
    args = ap.parse_args()

    gene_bed       = load_bed(args.omim_clingen_bed)       if args.omim_clingen_bed       else None
    gnomad_index   = load_bed(args.gnomad_sv)              if args.gnomad_sv              else None
    giab_index     = load_bed(args.giab_stratifications)   if args.giab_stratifications   else None

    print(f"[build_t2_scores] gene_bed={'loaded' if gene_bed else 'skipped'}, "
          f"gnomad={'loaded' if gnomad_index else 'skipped'}, "
          f"giab={'loaded' if giab_index else 'skipped'}", file=sys.stderr)

    n_scored = 0
    n_candidates = 0
    with args.in_path.open() as fin, \
         args.out_scored.open("w") as fout, \
         args.out_candidates.open("w") as cands:

        cands.write("locus_id\tarchitecture\tcoords\tscore\tevidence\tgenes\n")
        for line in fin:
            rec = json.loads(line)
            score, tags, genes = score_record(
                rec, args.assembly, gene_bed, gnomad_index, giab_index)
            rec["promotion_score"]    = score
            rec["promotion_evidence"] = tags
            rec["promotion_genes"]    = genes
            fout.write(json.dumps(rec, separators=(",", ":")) + "\n")
            n_scored += 1
            if score >= args.threshold:
                n_candidates += 1
                coords = rec.get("coords", {})
                asm_key = "GRCh38" if args.assembly == "grch38" else "CHM13"
                c = coords.get(asm_key) or {}
                cands.write("\t".join([
                    rec.get("locus_id", "?"),
                    rec.get("architecture", "?"),
                    f"{c.get('chrom', '?')}:{c.get('start','?')}-{c.get('end','?')}",
                    str(score),
                    ";".join(tags),
                    ";".join(genes),
                ]) + "\n")

    print(f"[build_t2_scores] {n_scored} scored, {n_candidates} candidates "
          f"(score ≥ {args.threshold}) → {args.out_candidates}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
