#!/usr/bin/env python3
"""Synthetic read generator for the LLmap cross-species benchmark harness.

Reads ``benchmarks/species_registry.json`` for organism profiles, then draws
reads from the organism's reference FASTA (auto-built via
``build_fake_reference.py`` if missing) according to a tier specification.

Outputs into ``--output-dir``:
  reads.fastq            (long-read tiers) or reads.R1.fastq + reads.R2.fastq (paired tiers)
  truth.tsv              read_id  source_contig  source_start  source_end  true_strand  read_class
  manifest.json          tier/organism/n_reads/total_bp/region_breakdown/seed/...

Tier scheme (per organism):
  1 :  1000 reads, mean  5 kb, HiFi-low-error
  2 :  1000 reads, mean 12 kb, HiFi
  3 :  1000 reads, mean  3 kb, ONT-like (5% error)
  4 : 10000 reads, mixed long-read profiles
  5 : same as 4 but reads concentrated on specific_loci regions (SD/centromere/telomere)
  6 : like 5 but inject SVs at NAHR-flank loci (DEL/DUP/INV)
  7 : 2x150 bp paired Illumina, low error
  8 : 2x100 bp paired short, ONT-style chimeras
  9 : 100% from synthetic repeat families
 10 : 10000-read full mix of all of the above

CLI:
  python gen_synth_reads.py --organism mouse --tier 5 --reads 10000 \
      --output-dir benchmarks/datasets/cache/mouse_t5
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Tier table
# ---------------------------------------------------------------------------

TIERS: Dict[int, Dict] = {
    1: {"reads": 1000,  "mean": 5000,  "sd": 1500, "min_len": 500,  "max_len": 15000,
        "profile": "hifi_low",   "paired": False, "region_bias": None,    "inject_sv": False, "force_repeat": False},
    2: {"reads": 1000,  "mean": 12000, "sd": 3000, "min_len": 500,  "max_len": 30000,
        "profile": "hifi",       "paired": False, "region_bias": None,    "inject_sv": False, "force_repeat": False},
    3: {"reads": 1000,  "mean": 3000,  "sd": 1000, "min_len": 300,  "max_len": 12000,
        "profile": "ont",        "paired": False, "region_bias": None,    "inject_sv": False, "force_repeat": False},
    4: {"reads": 10000, "mean": 8000,  "sd": 3000, "min_len": 300,  "max_len": 25000,
        "profile": "mixed_long", "paired": False, "region_bias": None,    "inject_sv": False, "force_repeat": False},
    5: {"reads": 10000, "mean": 8000,  "sd": 3000, "min_len": 300,  "max_len": 25000,
        "profile": "mixed_long", "paired": False, "region_bias": "hard",  "inject_sv": False, "force_repeat": False},
    6: {"reads": 10000, "mean": 8000,  "sd": 3000, "min_len": 300,  "max_len": 25000,
        "profile": "mixed_long", "paired": False, "region_bias": "nahr",  "inject_sv": True,  "force_repeat": False},
    7: {"reads": 10000, "mean": 150,   "sd": 0,    "min_len": 150,  "max_len": 150,
        "profile": "illumina",   "paired": True,  "insert": 500, "insert_sd": 50, "read_len": 150,
        "region_bias": None, "inject_sv": False, "force_repeat": False},
    8: {"reads": 10000, "mean": 100,   "sd": 0,    "min_len": 100,  "max_len": 100,
        "profile": "ont_short",  "paired": True,  "insert": 350, "insert_sd": 60, "read_len": 100,
        "chimera_rate": 0.05, "region_bias": None, "inject_sv": False, "force_repeat": False},
    9: {"reads": 10000, "mean": 2500,  "sd": 1000, "min_len": 300,  "max_len": 8000,
        "profile": "hifi",       "paired": False, "region_bias": None,   "inject_sv": False, "force_repeat": True},
    10: {"reads": 10000, "mean": 8000, "sd": 3000, "min_len": 200,  "max_len": 25000,
         "profile": "mixed_all", "paired": False, "region_bias": "mixed", "inject_sv": True,  "force_repeat": False},
}

# Profile -> (sub_rate, ins_rate, del_rate); None means "decide per read".
PROFILE_ERRORS: Dict[str, Tuple[Optional[float], Optional[float], Optional[float]]] = {
    "hifi_low":   (0.0005, 0.0001, 0.0001),
    "hifi":       (0.001,  0.0003, 0.0003),
    "ont":        (0.03,   0.01,   0.01),
    "ont_short":  (0.03,   0.01,   0.01),
    "illumina":   (0.001,  0.00005, 0.00005),
    "mixed_long": (None,   None,   None),
    "mixed_all":  (None,   None,   None),
}

NAHR_KEYWORDS = (
    "nahr", "sd_", "segmental_duplication", "flanking",
    "direct_orientation", "inverted_orientation", "duplicon",
)

ALPHABET = "ACGT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Contig:
    name: str
    seq: str

    def __len__(self) -> int:
        return len(self.seq)


@dataclass
class LocusRegion:
    contig: str
    start: int
    end: int
    klass: str
    name: str

    def length(self) -> int:
        return max(0, self.end - self.start)


# ---------------------------------------------------------------------------
# FASTA loader (no biopython dependency)
# ---------------------------------------------------------------------------

def load_fasta(path: Path) -> List[Contig]:
    contigs: List[Contig] = []
    name: Optional[str] = None
    parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    contigs.append(Contig(name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        contigs.append(Contig(name, "".join(parts).upper()))
    return contigs


# ---------------------------------------------------------------------------
# Region / class helpers
# ---------------------------------------------------------------------------

def classify_locus_name(name: str) -> str:
    n = name.lower()
    if "centromere" in n or "alpha_satellite" in n:
        return "centromere"
    if "telomere" in n or "subtelomeric" in n:
        return "telomere"
    if "segmental" in n or "_sd" in n or "sd_" in n:
        return "SD"
    if "vntr" in n or "trinuc" in n or "tandem" in n or "palindrome" in n or "homopolymer" in n:
        return "synth_stress"
    if "sine" in n or "line" in n or "ltr" in n or "repeat" in n or "low_complexity" in n:
        return "low_complexity"
    return "unique"


def load_regions(loci_index_path: Path) -> List[LocusRegion]:
    regions: List[LocusRegion] = []
    if not loci_index_path.exists():
        return regions
    with loci_index_path.open() as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            name, contig, s, e = parts[0], parts[1], parts[2], parts[3]
            try:
                start, end = int(s), int(e)
            except ValueError:
                continue
            regions.append(LocusRegion(contig, start, end, classify_locus_name(name), name))
    return regions


def filter_regions(regions: Sequence[LocusRegion], bias: str) -> List[LocusRegion]:
    if bias == "hard":
        return [r for r in regions if r.klass in ("SD", "centromere", "telomere", "low_complexity", "synth_stress")]
    if bias == "nahr":
        return [r for r in regions if r.klass == "SD" or any(k in r.name.lower() for k in NAHR_KEYWORDS)]
    if bias == "mixed":
        return list(regions)
    return list(regions)


def classify_by_overlap(contig: str, start: int, end: int, regions: List[LocusRegion]) -> str:
    best = "unique"
    for r in regions:
        if r.contig != contig:
            continue
        if r.end <= start or r.start >= end:
            continue
        if r.klass != "unique":
            return r.klass
        best = r.klass
    return best


# ---------------------------------------------------------------------------
# Read drawing
# ---------------------------------------------------------------------------

def pick_length(rng: random.Random, mean: int, sd: int, lo: int, hi: int) -> int:
    if sd <= 0:
        return max(lo, min(hi, mean))
    for _ in range(5):
        v = int(rng.gauss(mean, sd))
        if lo <= v <= hi:
            return v
    return max(lo, min(hi, mean))


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def errorize(rng: random.Random, seq: str, sub: float, ins: float, dele: float) -> str:
    if sub == 0.0 and ins == 0.0 and dele == 0.0:
        return seq
    out: List[str] = []
    for ch in seq:
        if rng.random() < dele:
            continue
        if rng.random() < sub:
            alt = rng.choice(ALPHABET.replace(ch, "")) if ch in ALPHABET else rng.choice(ALPHABET)
            out.append(alt)
        else:
            out.append(ch)
        if rng.random() < ins:
            out.append(rng.choice(ALPHABET))
    return "".join(out)


def per_read_error_rates(profile: str, rng: random.Random) -> Tuple[float, float, float]:
    if profile == "mixed_long":
        choices = ["hifi_low", "hifi", "ont"]
        weights = [0.3, 0.5, 0.2]
        return PROFILE_ERRORS[rng.choices(choices, weights=weights, k=1)[0]]  # type: ignore[return-value]
    if profile == "mixed_all":
        choices = ["hifi_low", "hifi", "ont", "illumina"]
        weights = [0.25, 0.4, 0.25, 0.1]
        return PROFILE_ERRORS[rng.choices(choices, weights=weights, k=1)[0]]  # type: ignore[return-value]
    return PROFILE_ERRORS[profile]  # type: ignore[return-value]


def random_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choices(ALPHABET, k=n))


# ---------------------------------------------------------------------------
# Region-weighted source picking
# ---------------------------------------------------------------------------

def draw_position_from_regions(
    rng: random.Random,
    regions: List[LocusRegion],
    contigs_by_name: Dict[str, Contig],
    read_len: int,
) -> Optional[Tuple[Contig, int, str]]:
    if not regions:
        return None
    weights = [max(1, r.length()) for r in regions]
    for _ in range(10):
        r = rng.choices(regions, weights=weights, k=1)[0]
        contig = contigs_by_name.get(r.contig)
        if contig is None:
            continue
        if len(contig) <= read_len:
            continue
        clo = max(0, r.start - read_len // 2)
        chi = max(clo + 1, min(len(contig) - read_len, r.end - read_len // 2))
        if chi <= clo:
            chi = clo + 1
        start = rng.randint(clo, chi)
        if start + read_len > len(contig):
            start = max(0, len(contig) - read_len)
        return contig, start, r.klass
    return None


def draw_position_uniform(
    rng: random.Random,
    contigs: List[Contig],
    read_len: int,
) -> Tuple[Contig, int, str]:
    weights = [max(1, len(c) - read_len) for c in contigs]
    for _ in range(10):
        c = rng.choices(contigs, weights=weights, k=1)[0]
        if len(c) <= read_len:
            continue
        start = rng.randint(0, len(c) - read_len)
        return c, start, "unique"
    c = contigs[0]
    return c, 0, "unique"


# ---------------------------------------------------------------------------
# SV injection
# ---------------------------------------------------------------------------

def inject_sv(rng: random.Random, seq: str) -> Tuple[str, str]:
    if len(seq) < 200:
        return seq, "none"
    sv_kind = rng.choice(["DEL", "DUP", "INV"])
    if sv_kind == "DEL":
        del_len = rng.randint(50, max(51, len(seq) // 4))
        s = rng.randint(0, len(seq) - del_len - 1)
        return seq[:s] + seq[s + del_len:], "DEL"
    if sv_kind == "DUP":
        dup_len = rng.randint(50, max(51, len(seq) // 4))
        s = rng.randint(0, len(seq) - dup_len - 1)
        return seq[:s + dup_len] + seq[s:s + dup_len] + seq[s + dup_len:], "DUP"
    if sv_kind == "INV":
        inv_len = rng.randint(50, max(51, len(seq) // 4))
        s = rng.randint(0, len(seq) - inv_len - 1)
        return seq[:s] + revcomp(seq[s:s + inv_len]) + seq[s + inv_len:], "INV"
    return seq, "none"


# ---------------------------------------------------------------------------
# Forced-repeat (tier 9) generator
# ---------------------------------------------------------------------------

def synth_repeat_read(rng: random.Random, length: int) -> Tuple[str, str]:
    kind = rng.choice(["homopolymer", "dinuc", "trinuc", "tandem", "palindrome"])
    if kind == "homopolymer":
        return rng.choice(ALPHABET) * length, "homopolymer"
    if kind == "dinuc":
        unit = rng.choice(["CA", "GT", "AT", "CG", "TG", "AC"])
    elif kind == "trinuc":
        unit = rng.choice(["CAG", "CGG", "GAA", "CTG", "GCC", "ATT"])
    elif kind == "tandem":
        ul = rng.randint(5, 12)
        unit = "".join(rng.choices(ALPHABET, k=ul))
    else:  # palindrome
        half = random_dna(rng, length // 2)
        seq = half + revcomp(half)
        return seq[:length], "palindrome"
    reps = (length // len(unit)) + 1
    return (unit * reps)[:length], kind


# ---------------------------------------------------------------------------
# Main generation
# ---------------------------------------------------------------------------

def generate(
    contigs: List[Contig],
    regions: List[LocusRegion],
    tier_cfg: Dict,
    n_reads: int,
    seed: int,
    organism: str,
    tier: int,
    error_profile_override: Optional[Dict[str, float]] = None,
    length_dist_override: Optional[Dict[str, int]] = None,
) -> Tuple[List[Dict], Dict[str, int]]:
    rng = random.Random(seed + 1000 * tier)
    contigs_by_name = {c.name: c for c in contigs}

    bias = tier_cfg.get("region_bias")
    profile = tier_cfg["profile"]
    paired = tier_cfg.get("paired", False)
    inject = tier_cfg.get("inject_sv", False)
    force_repeat = tier_cfg.get("force_repeat", False)
    chimera_rate = tier_cfg.get("chimera_rate", 0.0)

    mean = (length_dist_override or {}).get("mean", tier_cfg["mean"])
    sd = (length_dist_override or {}).get("sd", tier_cfg["sd"])
    lo = (length_dist_override or {}).get("min", tier_cfg["min_len"])
    hi = (length_dist_override or {}).get("max", tier_cfg["max_len"])
    # paired-tier read_len cannot be overridden by registry length distribution
    if paired:
        mean = tier_cfg["mean"]
        sd = tier_cfg["sd"]
        lo = tier_cfg["min_len"]
        hi = tier_cfg["max_len"]

    biased_regions = filter_regions(regions, bias) if bias else []

    breakdown: Dict[str, int] = {}
    records: List[Dict] = []

    for i in range(n_reads):
        rid = f"{organism}_t{tier}_{i}"

        # tier 9: pure synthetic repeats
        if force_repeat:
            rl = pick_length(rng, mean, sd, lo, hi)
            seq, subclass = synth_repeat_read(rng, rl)
            sub_e, ins_e, del_e = PROFILE_ERRORS["hifi"]  # type: ignore[assignment]
            seq = errorize(rng, seq, sub_e, ins_e, del_e)  # type: ignore[arg-type]
            records.append({
                "id": rid, "seq": seq, "qual_base": "I",
                "contig": "SYNTHETIC", "start": 0, "end": len(seq),
                "strand": "+", "klass": "synth_stress", "subclass": subclass,
            })
            breakdown["synth_stress"] = breakdown.get("synth_stress", 0) + 1
            continue

        # paired short tiers
        if paired:
            read_len = tier_cfg["read_len"]
            insert = max(read_len * 2 + 20, int(rng.gauss(tier_cfg["insert"], tier_cfg["insert_sd"])))

            pos = draw_position_from_regions(rng, biased_regions, contigs_by_name, insert) if biased_regions else None
            if pos is None:
                pos = draw_position_uniform(rng, contigs, insert)
            contig, start, klass = pos
            end = min(len(contig), start + insert)
            insert_seq = contig.seq[start:end]
            if len(insert_seq) < read_len * 2:
                insert_seq = insert_seq + random_dna(rng, read_len * 2 - len(insert_seq))

            chimera_tag = ""
            if chimera_rate > 0 and rng.random() < chimera_rate:
                c2 = rng.choice(contigs)
                if len(c2) > read_len:
                    s2 = rng.randint(0, len(c2) - read_len)
                    splice_at = len(insert_seq) // 2
                    insert_seq = insert_seq[:splice_at] + c2.seq[s2:s2 + (len(insert_seq) - splice_at)]
                    chimera_tag = f";chimera={c2.name}:{s2}"

            r1 = insert_seq[:read_len]
            r2 = revcomp(insert_seq[-read_len:])
            strand1, strand2 = "+", "-"
            if rng.random() < 0.5:
                r1, r2 = r2, r1
                strand1, strand2 = strand2, strand1

            sub_e, ins_e, del_e = per_read_error_rates(profile, rng)
            if error_profile_override:
                sub_e = error_profile_override.get("sub_pct", sub_e)
                ins_e = error_profile_override.get("ins_pct", ins_e)
                del_e = error_profile_override.get("del_pct", del_e)
            r1e = errorize(rng, r1, sub_e, ins_e, del_e)
            r2e = errorize(rng, r2, sub_e, ins_e, del_e)
            klass_final = klass if not chimera_tag else "chimera"

            records.append({
                "id": rid + "/1", "seq": r1e, "qual_base": "I",
                "contig": contig.name, "start": start, "end": end,
                "strand": strand1, "klass": klass_final, "subclass": f"R1{chimera_tag}",
            })
            records.append({
                "id": rid + "/2", "seq": r2e, "qual_base": "I",
                "contig": contig.name, "start": start, "end": end,
                "strand": strand2, "klass": klass_final, "subclass": f"R2{chimera_tag}",
            })
            breakdown[klass_final] = breakdown.get(klass_final, 0) + 2
            continue

        # long-read tiers
        rl = pick_length(rng, mean, sd, lo, hi)

        pos = draw_position_from_regions(rng, biased_regions, contigs_by_name, rl) if biased_regions else None
        if pos is None:
            pos = draw_position_uniform(rng, contigs, rl)
        contig, start, klass = pos
        end = min(len(contig), start + rl)
        seq = contig.seq[start:end]
        strand = "+"
        if rng.random() < 0.5:
            seq = revcomp(seq)
            strand = "-"

        sv_tag = ""
        if inject and rng.random() < 0.4:
            seq, sv_tag = inject_sv(rng, seq)

        sub_e, ins_e, del_e = per_read_error_rates(profile, rng)
        if error_profile_override:
            sub_e = error_profile_override.get("sub_pct", sub_e)
            ins_e = error_profile_override.get("ins_pct", ins_e)
            del_e = error_profile_override.get("del_pct", del_e)
        seq_e = errorize(rng, seq, sub_e, ins_e, del_e)

        if klass == "unique" and regions:
            klass = classify_by_overlap(contig.name, start, end, regions)

        if sv_tag and sv_tag != "none":
            klass = "SV_" + sv_tag

        subclass = sv_tag if sv_tag else "long"
        records.append({
            "id": rid, "seq": seq_e, "qual_base": "I",
            "contig": contig.name, "start": start, "end": end,
            "strand": strand, "klass": klass, "subclass": subclass,
        })
        breakdown[klass] = breakdown.get(klass, 0) + 1

    return records, breakdown


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def write_fastq(path: Path, records: List[Dict], filter_suffix: Optional[str] = None) -> int:
    n = 0
    with path.open("w") as fh:
        for r in records:
            if filter_suffix and not r["id"].endswith(filter_suffix):
                continue
            seq = r["seq"]
            fh.write(f"@{r['id']}\n{seq}\n+\n{r['qual_base'] * len(seq)}\n")
            n += 1
    return n


def write_truth(path: Path, records: List[Dict]) -> None:
    with path.open("w") as fh:
        fh.write("read_id\tsource_contig\tsource_start\tsource_end\ttrue_strand\tread_class\tsubclass\n")
        for r in records:
            fh.write(
                f"{r['id']}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{r['klass']}\t{r['subclass']}\n"
            )


def write_manifest(path: Path, payload: Dict) -> None:
    path.write_text(json.dumps(payload, indent=2))


# ---------------------------------------------------------------------------
# Auto-build reference if missing
# ---------------------------------------------------------------------------

def maybe_build_reference(
    organism: str,
    cfg: Dict,
    root: Path,
    seed: int,
) -> Path:
    ref_rel = cfg["reference"]
    ref_path = (root / ref_rel) if not os.path.isabs(ref_rel) else Path(ref_rel)
    if ref_path.exists() and (ref_path.parent / "loci_index.tsv").exists():
        return ref_path
    out_dir = ref_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(Path(__file__).parent))
    import build_fake_reference as bfr  # noqa: E402

    loci_dir = root / cfg["loci_dir"]
    sub_organisms = cfg.get("sub_organisms")
    loci = bfr.collect_loci(loci_dir, cfg["assembly_key"], sub_organisms=sub_organisms, registry_root=root)
    if not loci:
        # placeholder so reads can still be drawn (uniform random sequence)
        rng = random.Random(seed)
        placeholder_seq = "".join(rng.choices(list(ALPHABET), k=200_000))
        with ref_path.open("w") as fh:
            fh.write(f">{organism}_placeholder\n")
            for i in range(0, len(placeholder_seq), 80):
                fh.write(placeholder_seq[i:i + 80] + "\n")
        (out_dir / "loci_index.tsv").write_text(
            "locus_name\tcontig\tmapped_start\tmapped_end\tsource_path\n"
        )
        return ref_path
    contigs, index = bfr.build_contigs(loci, seed=seed)
    bfr.write_fasta(ref_path, contigs)
    with (out_dir / "loci_index.tsv").open("w") as fh:
        fh.write("locus_name\tcontig\tmapped_start\tmapped_end\tsource_path\n")
        for row in index:
            fh.write("\t".join(str(c) for c in row) + "\n")
    return ref_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Generate synthetic reads for a registered organism + tier.")
    ap.add_argument("--organism", required=True)
    ap.add_argument("--tier", type=int, required=True, choices=list(TIERS.keys()))
    ap.add_argument("--reads", type=int, default=0, help="override n_reads (0 = use tier default)")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--registry", default=None, help="path to species_registry.json")
    ap.add_argument("--root", default=None, help="llmap-local repo root (default: parent of benchmarks/)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    benchmarks_dir = script_dir.parent
    root = Path(args.root).resolve() if args.root else benchmarks_dir.parent
    registry_path = Path(args.registry).resolve() if args.registry else benchmarks_dir / "species_registry.json"

    registry = json.loads(registry_path.read_text())
    if args.organism not in registry:
        avail = sorted(k for k in registry if not k.startswith("_"))
        sys.stderr.write(f"organism {args.organism!r} not in registry; have: {avail}\n")
        return 2
    cfg = registry[args.organism]

    tier_cfg = dict(TIERS[args.tier])
    n_reads = args.reads if args.reads > 0 else tier_cfg["reads"]

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_path = maybe_build_reference(args.organism, cfg, root, args.seed)
    loci_index_path = ref_path.parent / "loci_index.tsv"
    contigs = load_fasta(ref_path)
    if not contigs:
        sys.stderr.write(f"reference at {ref_path} is empty\n")
        return 3
    regions = load_regions(loci_index_path)

    records, breakdown = generate(
        contigs=contigs,
        regions=regions,
        tier_cfg=tier_cfg,
        n_reads=n_reads,
        seed=args.seed,
        organism=args.organism,
        tier=args.tier,
        error_profile_override=cfg.get("error_profile"),
        length_dist_override=cfg.get("length_distribution"),
    )

    if tier_cfg.get("paired"):
        n1 = write_fastq(out_dir / "reads.R1.fastq", records, filter_suffix="/1")
        n2 = write_fastq(out_dir / "reads.R2.fastq", records, filter_suffix="/2")
        n_written = n1 + n2
    else:
        n_written = write_fastq(out_dir / "reads.fastq", records)

    write_truth(out_dir / "truth.tsv", records)

    manifest = {
        "organism": args.organism,
        "tier": args.tier,
        "tier_profile": tier_cfg,
        "n_reads_requested": n_reads,
        "n_records_written": n_written,
        "total_bp": sum(len(r["seq"]) for r in records),
        "reference_fa": str(ref_path),
        "loci_index": str(loci_index_path) if loci_index_path.exists() else None,
        "region_breakdown": breakdown,
        "seed": args.seed,
        "paired": bool(tier_cfg.get("paired")),
        "assembly_key": cfg["assembly_key"],
        "error_profile_override": cfg.get("error_profile"),
        "length_distribution_override": cfg.get("length_distribution"),
    }
    write_manifest(out_dir / "manifest.json", manifest)

    print(json.dumps({"ok": True, "manifest": str(out_dir / "manifest.json"),
                      "n_records_written": n_written,
                      "region_breakdown": breakdown}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
