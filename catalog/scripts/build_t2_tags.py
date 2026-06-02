#!/usr/bin/env python3
"""Enrich T2 bulk records with evidence tags. NO numeric score.

Rationale (per [[feedback_avoid_artificial_scoring]]):
- numeric "+5 for OMIM, +2 for gnomAD" weights are arbitrary
- thresholds like "score >= 5" hide which evidence actually triggered promotion
- per-locus tag combinations are transparent and reproducible

What this script does:
- Annotates each T2 record with a flat list of evidence tags (boolean per tag)
- Records the actual evidence values (gnomAD-SV AF in percent, gene names,
  GIAB stratification label) — not converted to points
- Emits a sorted candidate list grouped by tag combination, so a curator
  can scan "clinically-relevant ∩ immune-receptor" before
  "polymorphic-only" candidates

Available evidence tags:
    immune_receptor     gene-symbol matches IGH/IGK/IGL/TR*
    mhc_block           gene-symbol matches HLA-*/MICA/MICB/TAP1/TAP2
    paraphase_curated   source = Paraphase config.yaml (already vetted by PacBio)
    clinically_relevant gene overlap with clinical-disease bed (ClinGen/OMIM)
    sv_polymorphic      overlap with gnomAD-SV sites (records AF if available)
    low_mappability     overlap with GIAB low-mappability stratification

Output formats:
  catalog/bulk/v<release>.<assembly>.tagged.jsonl   (T2 + tags + evidence detail)
  catalog/bulk/v<release>.candidates_by_tagset.tsv  (grouped sortable view)

Usage on <hpc-login>:
    python3 build_t2_tags.py \\
        --release v2026.Q2 \\
        --in catalog/bulk/v2026.Q2.grch38.jsonl \\
        --assembly grch38 \\
        --clinical-bed <optional path> \\
        --gnomad-sv-bed <optional path> \\
        --giab-lowmap-bed <optional path> \\
        --out-tagged catalog/bulk/v2026.Q2.grch38.tagged.jsonl \\
        --out-candidates catalog/bulk/v2026.Q2.candidates_by_tagset.tsv
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path


# Immune-receptor / antibody locus prefixes (per IMGT / GENCODE biotype IG_*_gene).
# Match must be longer than 3 chars to avoid TRAPPC*/IGLL*/IGHMBP2 false positives;
# real IG/TR symbols have a V/D/J/C segment letter at position 4 (e.g. IGHG4, IGKV1-5,
# TRAV12, TRBJ2-1) — except a handful of canonical single-segment locus aliases.
import re as _re
_IMMUNE_RE = _re.compile(
    r"^(IGH[VDJCAEGM]|IGK[VJC]|IGL[VJC]|TR[ABGD][VDJC])(\d|-|$)"
)

# MHC genes (block-level)
MHC_GENE_PREFIXES = ("HLA-", "HLA_", "MICA", "MICB", "TAP1", "TAP2")


def gene_symbols_from_record(rec: dict) -> list[str]:
    """Pull gene symbols out of `notes` (Paraphase carries them there)."""
    notes = rec.get("notes", "") or ""
    out: list[str] = []
    for token in notes.split(";"):
        t = token.strip()
        if t.startswith("paraphase_gene="):
            out.append(t.split("=", 1)[1])
        elif t.startswith("paralog_set="):
            out.extend(s.strip() for s in t.split("=", 1)[1].split(","))
    return out


# ----------------------------------------------------------------------------
# Lightweight BED index (per-chrom sorted list, linear scan within chrom)
# ----------------------------------------------------------------------------

class IntervalIndex:
    def __init__(self) -> None:
        self.by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        self._sorted = False

    def add(self, chrom: str, start: int, end: int, label: str = "") -> None:
        self.by_chrom[chrom].append((start, end, label))

    def finalize(self) -> None:
        for k in self.by_chrom:
            self.by_chrom[k].sort()
        self._sorted = True

    def hits(self, chrom: str, start: int, end: int) -> list[str]:
        if not self._sorted:
            self.finalize()
        out: list[str] = []
        for s, e, lab in self.by_chrom.get(chrom, []):
            if e < start:
                continue
            if s > end:
                break
            out.append(lab)
        return out


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
# Tag a single record
# ----------------------------------------------------------------------------

def tag_record(rec: dict,
                assembly_key: str,
                clinical_bed: IntervalIndex | None,
                gnomad_index: IntervalIndex | None,
                giab_index: IntervalIndex | None) -> tuple[list[str], dict]:
    """Return (sorted unique tag list, evidence dict)."""
    tags: list[str] = []
    evidence: dict = {}

    coords = rec.get("coords") or {}
    # Records may key coords by assembly ('GRCh38') or by chromosome ('chr1').
    # Pick the assembly-keyed entry if present; otherwise tag against every
    # chrom-keyed entry (some loci span multiple chromosomes / paired SDs).
    coord_entries: list[dict] = []
    if assembly_key in coords:
        coord_entries.append(coords[assembly_key])
    else:
        for k, v in coords.items():
            if isinstance(v, dict) and "chrom" in v and "start" in v and "end" in v:
                coord_entries.append(v)
    if not coord_entries:
        return tags, evidence

    # 1) Gene-symbol-driven tags (work without external beds)
    syms = gene_symbols_from_record(rec)
    for s in syms:
        up = s.upper()
        if _IMMUNE_RE.match(up):
            tags.append("immune_receptor")
            evidence.setdefault("immune_receptor_genes", []).append(s)
        if any(up.startswith(pfx) for pfx in MHC_GENE_PREFIXES):
            tags.append("mhc_block")
            evidence.setdefault("mhc_genes", []).append(s)

    # 2) Source-driven tag
    if "paraphase" in (rec.get("source") or "").lower():
        tags.append("paraphase_curated")

    # 3) BED-driven tags (loop over every coord entry; aggregate hits)
    for c in coord_entries:
        chrom, start, end = c["chrom"], c["start"], c["end"]
        if clinical_bed is not None:
            hits = clinical_bed.hits(chrom, start, end)
            if hits:
                tags.append("clinically_relevant")
                evidence.setdefault("clinical_genes", [])
                evidence["clinical_genes"].extend(h for h in hits if h)
                # immune-receptor / mhc auto-detect from clinical overlap too
                for h in hits:
                    up = (h or "").upper()
                    if _IMMUNE_RE.match(up):
                        tags.append("immune_receptor")
                    if any(up.startswith(pfx) for pfx in MHC_GENE_PREFIXES):
                        tags.append("mhc_block")
        if gnomad_index is not None:
            hits = gnomad_index.hits(chrom, start, end)
            if hits:
                tags.append("sv_polymorphic")
                evidence.setdefault("gnomad_sv_hits", [])
                # keep at most 5 per record (already limited below)
                for h in hits:
                    if len(evidence["gnomad_sv_hits"]) >= 5:
                        break
                    evidence["gnomad_sv_hits"].append(h)
        if giab_index is not None:
            hits = giab_index.hits(chrom, start, end)
            if hits:
                tags.append("low_mappability")
                evidence.setdefault("giab_labels", [])
                evidence["giab_labels"].extend(h for h in hits if h)

    # de-dup clinical_genes and giab_labels (limit to 5 for compactness)
    if "clinical_genes" in evidence:
        evidence["clinical_genes"] = sorted(set(evidence["clinical_genes"]))
    if "giab_labels" in evidence:
        evidence["giab_labels"] = sorted(set(evidence["giab_labels"]))[:5]

    return sorted(set(tags)), evidence


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release",   required=True)
    ap.add_argument("--in",        dest="in_path", required=True, type=Path)
    ap.add_argument("--assembly",  required=True,
                    choices=["grch38", "chm13", "paraphase"])
    ap.add_argument("--clinical-bed",   type=Path, default=None,
                    help="BED of clinical-disease loci (col 4 = gene name)")
    ap.add_argument("--gnomad-sv-bed",  type=Path, default=None,
                    help="gnomAD-SV sites as BED (col 4 = AF or ID)")
    ap.add_argument("--giab-lowmap-bed",type=Path, default=None,
                    help="GIAB low-mappability BED")
    ap.add_argument("--out-tagged",     required=True, type=Path)
    ap.add_argument("--out-candidates", required=True, type=Path)
    args = ap.parse_args()

    asm_key = {"grch38": "GRCh38", "chm13": "CHM13",
               "paraphase": "GRCh38"}[args.assembly]

    clinical = load_bed(args.clinical_bed)   if args.clinical_bed   else None
    gnomad   = load_bed(args.gnomad_sv_bed)  if args.gnomad_sv_bed  else None
    giab     = load_bed(args.giab_lowmap_bed)if args.giab_lowmap_bed else None

    print(f"[build_t2_tags] clinical={'on' if clinical else 'off'}, "
          f"gnomad={'on' if gnomad else 'off'}, "
          f"giab={'on' if giab else 'off'}", file=sys.stderr)

    by_tagset: dict[tuple, list[dict]] = defaultdict(list)
    n_total = 0
    with args.in_path.open() as fin, args.out_tagged.open("w") as fout:
        for line in fin:
            rec = json.loads(line)
            tags, evidence = tag_record(rec, asm_key, clinical, gnomad, giab)
            rec["evidence_tags"] = tags
            rec["evidence_detail"] = evidence
            fout.write(json.dumps(rec, separators=(",", ":")) + "\n")
            n_total += 1
            if tags:
                by_tagset[tuple(tags)].append(rec)

    # Emit candidates grouped + sorted: most-tagged combinations first
    rows_by_combo = sorted(by_tagset.items(),
                           key=lambda kv: (-len(kv[0]), kv[0]))
    with args.out_candidates.open("w") as out:
        out.write("tagset\tn_loci\tlocus_id\tarchitecture\tcoords\tgenes\n")
        for tags, recs in rows_by_combo:
            for rec in recs:
                coords = rec.get("coords") or {}
                c = coords.get(asm_key)
                if not c:
                    # fall back to first chrom-keyed entry
                    for v in coords.values():
                        if isinstance(v, dict) and "chrom" in v:
                            c = v
                            break
                c = c or {}
                coord_str = f"{c.get('chrom','?')}:{c.get('start','?')}-{c.get('end','?')}"
                genes = []
                ev = rec.get("evidence_detail") or {}
                genes.extend(ev.get("immune_receptor_genes", []))
                genes.extend(ev.get("mhc_genes", []))
                genes.extend(ev.get("clinical_genes", []))
                # Backward compatible: read structural_architecture, fall back to legacy 'architecture'
                arch = rec.get("structural_architecture") or rec.get("architecture", "?")
                out.write("\t".join([
                    "|".join(tags),
                    str(len(recs)),
                    rec.get("locus_id", "?"),
                    arch,
                    coord_str,
                    ",".join(sorted(set(genes)))[:200],
                ]) + "\n")

    print(f"[build_t2_tags] {n_total} records tagged "
          f"({sum(len(v) for v in by_tagset.values())} with at least one tag, "
          f"{len(by_tagset)} distinct tag combinations)", file=sys.stderr)
    print(f"[build_t2_tags] candidates by tagset → {args.out_candidates}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
