#!/usr/bin/env python3
"""Convert the 3318-row SD-pair linter TSV into per-pair JSON files matching
the LLmap nahr_flanks schema. Merges disease annotations from
genome_wide_nahr_disease.tsv (multi-row per pair_id) into a single records.

Output: one JSON per pair under linter_3318/pairs/, plus index.json.
"""
import csv, json, os, sys
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent
RAW  = ROOT / "raw"
OUT  = ROOT / "pairs"
OUT.mkdir(parents=True, exist_ok=True)

# ---- Load disease annotations (pair_id -> list of gene dicts) ---------------
disease_by_pair = defaultdict(list)
with open(RAW / "genome_wide_nahr_disease.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        pid = row["pair_id"]
        disease_by_pair[pid].append({
            "gene": row["gene"],
            "gene_type": row["gene_type"],
            "ClinGen_HI": row.get("ClinGen_HI", "") or None,
            "ClinGen_TS": row.get("ClinGen_TS", "") or None,
            "HPO_n_diseases": int(row.get("HPO_n_diseases", "0") or 0),
            "HPO_top_disease": row.get("HPO_top_disease", "") or None,
            "ClinVar_n_phenotypes": int(row.get("ClinVar_n_phenotypes", "0") or 0),
        })

# ---- Load top unique regions (rank, n_HI3, syndrome) -----------------------
top_by_pair = {}
with open(RAW / "top_unique_nahr_regions.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        top_by_pair[row["region_id"]] = {
            "rank": int(row["rank"]),
            "n_genes": int(row["n_genes"]),
            "n_disease": int(row["n_disease"]),
            "n_HI3": int(row["n_HI3"]),
            "HI_genes": row["HI_genes"],
            "known_syndrome": row.get("known_syndrome", "") or None,
            "omim_examples": row.get("omim_examples", "") or None,
        }

# ---- Map identity+propensity -> mapping_hints ------------------------------
def hints(identity: float, propensity: float, lcr_len: int):
    """Map physical SD properties to LLmap hints.
    identity in [0.90, 1.00], propensity X in ~[2..10], lcr_len in bp.
    Higher identity -> need lower anchor weight + higher lambda.
    """
    # lambda_scale: from 1.2 (low identity) to 1.8 (>= 0.995)
    if identity >= 0.995: lam = 1.8
    elif identity >= 0.99: lam = 1.6
    elif identity >= 0.98: lam = 1.5
    elif identity >= 0.96: lam = 1.4
    elif identity >= 0.93: lam = 1.3
    else:                  lam = 1.2
    # anchor_weight_scale: lower for very high identity
    if identity >= 0.995: aws = 0.25
    elif identity >= 0.99: aws = 0.30
    elif identity >= 0.98: aws = 0.35
    elif identity >= 0.96: aws = 0.40
    else:                  aws = 0.50
    # max_occ: scale with LCR length and identity
    if identity >= 0.99 and lcr_len > 200_000: mo = 120000
    elif identity >= 0.99:                       mo = 100000
    elif identity >= 0.97:                       mo = 80000
    elif identity >= 0.94:                       mo = 60000
    else:                                        mo = 40000
    # high propensity -> push lambda one notch higher (more SD-mediated risk)
    if propensity >= 5.0: lam = min(1.9, lam + 0.1)
    return {
        "lambda_scale": round(lam, 3),
        "anchor_weight_scale": round(aws, 3),
        "max_occ": mo,
        "report_multi_position": True,
        "require_psv_disambig": identity >= 0.95,
        "allow_high_mismatch": False,
    }

# ---- Build per-pair JSONs ---------------------------------------------------
manifest = []
with open(RAW / "raw_3318_sd_pairs.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for row in rd:
        pid = row["pair_id"]
        chrom = row["chrom"]
        lcr_up_s = int(row["lcr_up_start"])
        lcr_up_e = int(row["lcr_up_end"])
        lcr_dn_s = int(row["lcr_down_start"])
        lcr_dn_e = int(row["lcr_down_end"])
        ints_s   = int(row["interior_start"])
        ints_e   = int(row["interior_end"])
        ints_len = int(row["interior_len_bp"])
        avg_lcr  = int(row["avg_lcr_len_bp"])
        identity = float(row["identity"])
        propX    = float(row["propensity_X"])
        # Coordinates: the full extent from LCR up start to LCR down end is the
        # full NAHR-vulnerable window.
        full_start = min(lcr_up_s, lcr_dn_s)
        full_end   = max(lcr_up_e, lcr_dn_e)
        diseases   = disease_by_pair.get(pid, [])
        hi3_genes  = sorted({d["gene"] for d in diseases
                              if d.get("ClinGen_HI") == "3"})
        ts3_genes  = sorted({d["gene"] for d in diseases
                              if d.get("ClinGen_TS") == "3"})
        top_info   = top_by_pair.get(pid, {})
        is_top     = bool(top_info)

        rec = {
            "name": pid,
            "taxonomy_id": "paralog_family",
            "description": (
                f"Linter NAHR-prone interval (genome-wide SD-pair scan v2). "
                f"{chrom}:{ints_s}-{ints_e} ({ints_len/1e6:.2f} Mb interior) "
                f"flanked by paralogous {avg_lcr/1000:.1f} kb LCRs at "
                f"{identity*100:.2f}% identity, propensity X={propX:.2f}. "
                + (f"Known syndrome: {top_info.get('known_syndrome')}. " if is_top and top_info.get('known_syndrome') else "")
                + (f"HI=3 genes: {','.join(hi3_genes[:6])}. " if hi3_genes else "")
            ).strip(),
            "coordinates": {
                "grch38": {"chr": chrom, "start": full_start, "end": full_end}
            },
            "size_bp": full_end - full_start,
            "nahr_class": "recurrent_del_dup",
            "architecture": {
                "type": "direct_orientation_block",
                "BP_mediating_SDs": [
                    f"{chrom}:{lcr_up_s}-{lcr_up_e}",
                    f"{chrom}:{lcr_dn_s}-{lcr_dn_e}",
                ],
                "homology_percent": round(identity * 100, 3),
                "sd_length_kb": round(avg_lcr / 1000.0, 2),
                "interior_bp": ints_len,
                "propensity_X": propX,
                "internal_structure": ", ".join(
                    sorted({d["gene"] for d in diseases
                            if d.get("gene_type") == "protein_coding"})[:20]
                ) or "no protein-coding genes annotated",
                "copy_number_in_population": "variable_CNV",
            },
            "mapping_hints": hints(identity, propX, avg_lcr),
            "related_diseases": sorted({
                d["HPO_top_disease"] for d in diseases
                if d.get("HPO_top_disease")
            })[:20],
            "recurrence_evidence": {
                "gnomad_sv_count": None,
                "dbvar_count": None,
                "clingen_HI3_genes": hi3_genes,
                "clingen_TS3_genes": ts3_genes,
                "n_protein_coding_genes": sum(
                    1 for d in diseases if d.get("gene_type") == "protein_coding"
                ),
                "n_total_features": len(diseases),
                "top_region_rank": top_info.get("rank"),
                "known_syndrome": top_info.get("known_syndrome"),
                "tumor_recurrence": None,
            },
            "confidence": 0.95 if is_top else (0.85 if hi3_genes else 0.7),
            "sources": [
                "schlein-lab eye_evolution NAHR linter v2 (genome_wide_nahr_v2)",
                "raw_3318_sd_pairs.tsv 2026-05",
            ],
        }
        # Write per-pair JSON
        out_path = OUT / f"{pid}.json"
        with open(out_path, "w") as g:
            json.dump(rec, g, indent=2)
        manifest.append({
            "pair_id": pid, "chrom": chrom, "size_bp": rec["size_bp"],
            "identity": identity, "propensity_X": propX,
            "is_top_region": is_top, "n_HI3": len(hi3_genes),
        })

with open(ROOT / "index.json", "w") as g:
    json.dump({
        "schema_version": "1",
        "source": "schlein-lab eye_evolution NAHR linter v2",
        "n_pairs": len(manifest),
        "pairs": manifest,
    }, g, indent=2)

print(f"Wrote {len(manifest)} per-pair JSONs to {OUT}")
print(f"Index: {ROOT / 'index.json'}")
