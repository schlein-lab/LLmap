# LLmap SegDup Catalog — Specification

**Status:** Draft v0.1 — 2026-06-02
**Decision context:** [intrinsic detection rejected](../../memory/llmap_intrinsic_segdup_detector.md); curated catalog with database-backed provenance is the agreed approach. See companion docs:
- [segdup_databases_inventory.md](segdup_databases_inventory.md) — source database analysis
- [segdup_class_taxonomy.md](segdup_class_taxonomy.md) — 3-axis classification scheme

## 1. Three-Tier Architecture

Total SD content (~5–6.6 % of mammalian genomes, ~5 000–10 000 distinct loci) is too large for full manual curation. The catalog is split into three tiers:

| Tier | Count | Granularity | Population | Mapping behaviour |
|------|------:|---|---|---|
| **T1 Curated** | ~50–300 | per-locus JSON, full PSV catalogue, mapping strategy, DOI provenance | manual / PR-based | locus-specific resolver, PSV disambiguation, class-aware fallback chain |
| **T2 Bulk** | ~5 000–10 000 | one JSONL line per locus, coords + 3-axis tag + DB-provenance | auto-import from UCSC / T2T / GIAB / gnomAD-SV | "be careful": keep top-K hits, lower MAPQ floor, no class-specific PSV step |
| **T3 Detect** | rest of genome | none | on-the-fly k-mer multiplicity anomaly | warning + reduced MAPQ, no special routing |

T1 grows organically only where evidence justifies it. T2 grows automatically each catalog release. T3 is mapper-internal, not stored in the catalog.

### 1.1 Promotion criterion (T2 → T1)

A locus is a T1-promotion candidate when the priority score meets a threshold:

- OMIM/ClinGen disease gene overlap: +5
- Antibody / immune-receptor locus: **auto-T1** (no score threshold)
- gnomAD-SV allele frequency > 1 % in any population: +2
- HPRC pangenome haplotype-variability flag: +2
- GIAB low-mappability flag: +1

Suggested threshold: score ≥ 5. Promotion requires a curator-reviewed PR.

## 1.2 Classification axes — two-axis structural × haplotype split (since v0.2)

A locus is classified along **two independent axes**, not one. This separation is essential because the same structural architecture (e.g. NoDup single-copy) can carry different haplotype-level PSV patterns (e.g. canonical-PSV vs chimdup-PSV in IGHG4) — and the same haplotype-class can appear on different structural architectures.

| Axis | What it captures | Schema field | Vocabulary |
|---|---|---|---|
| **Structural** | copy number + spatial arrangement | `structural_architecture` (one value, enum) | `single_copy`, `tandem_dup_small/large`, `intrachromosomal_interspersed`, `interchromosomal`, `subtelomeric`, `pericentromeric`, `palindromic_inverted`, `inverted_dispersed`, `complex_mosaic`, `unresolved` |
| **Haplotype** | which PSV pattern sits on the copies | `haplotype_class` (free-form string) | locus-specific, follow `<locus_short>_<descriptor>` convention. Meaning defined in `diagnostic_features.discriminating_snps` of the same entry. |

`haplotype_class` is **intentionally not a global enum** — each locus has its own PSV combinatorics, and forcing them into a global vocabulary would either bloat the enum or lose information. Validator only checks `^[A-Za-z][A-Za-z0-9_]+$` (non-empty identifier). See [feedback_avoid_artificial_scoring](../../memory/feedback_avoid_artificial_scoring.md) — same principle: no global classification where locus-specificity is real.

T2 bulk records carry `structural_architecture` (inferred from coordinates) but **not** `haplotype_class` — the latter requires PSV / sequence-level resolution which T2 doesn't have.

### 1.2.1 IGHG4 worked example (three T1 entries from the same 295-bp CH1)

| T1 entry | structural_architecture | haplotype_class | n haps in 104-sample sweep |
|---|---|---|---:|
| `IGHG4_chimdup_tandem` | `tandem_dup_small` | `IGHG4_chimdup_homozygous` | 22 |
| `IGHG_canondup_nahr_block` | `tandem_dup_large` | `IGHG_canonical_homozygous` | 15 (14/15 canonical, 1 mixed) |
| `IGHG4_chimdup_canonical_arch` | `single_copy` | `IGHG4_chimdup_single_copy` | 39 |

Three structurally and haplotype-distinct classes from one ~120 kb region — would have been collapsed into one `IGHG4_chimdup` entry under the single-axis scheme.

## 2. T1 (Curated) JSON Schema

One file per locus at `catalog/curated/<locus_id>.json`. Required fields below; full JSON-Schema lives at `catalog/schema/curated.schema.json`.

```jsonc
{
  "locus_id": "IGHG4_chimdup_tandem",          // snake_case, stable
  "human_name": "IGHG4 chimeric tandem duplication",
  "version": "v2026.Q2",                       // catalog release this entry belongs to
  "schema_version": "0.2",

  // === Two-axis structural × haplotype classification (§1.2) ===
  "structural_architecture": "tandem_dup_small",       // exactly one (enum, axis A)
  "haplotype_class": "IGHG4_chimdup_homozygous",       // free-form locus-specific tag
  "mechanism": ["gene_conversion_prone"],              // zero or more (axis B)
  "clinical_function": ["immune_receptor_locus", "pangenome_only"], // zero or more (axis C)
  "nahr_status": "ambiguous_low_identity",             // free-form refinement of mechanism

  // === Coordinates per reference assembly ===
  "coords": {
    "GRCh38": { "chrom": "chr14", "start": 105625772, "end": 105645272, "strand": "-" },
    "CHM13":  { "chrom": "chr14", "start":  99891444, "end":  99910944, "strand": "-" },
    "HPRC_pangenome_v1.1": null                // null if not yet resolved
  },

  // === Structural parameters ===
  "unit_size_bp": 19500,
  "copy_count_typical": "2-3",
  "identity_pct_within_unit": 99.99,           // see ighg4_sgamma4_identical_in_tandem_dup.md
  "boundary_sharpness": "sharp",               // sharp | gradual | unknown

  // === Diagnostic features (what discriminates copies) ===
  "diagnostic_features": {
    "discriminating_snps": [
      { "id": "IGHG4_CH1_pos69",
        "position": { "transcript": "IGHG4-201", "cds_offset": 69, "grch38_chr14": 105625998 },
        "ref": "C", "alt": "G", "consequence": "p.Ala23=",
        "class": "DUP_fixed", "frequency_in_carriers": 1.0,
        "source": "this work + per_read_genotyping" },
      { "id": "IGHG4_CH1_pos252",
        "position": { "transcript": "IGHG4-201", "cds_offset": 252, "grch38_chr14": 105625815 },
        "ref": "C", "alt": "T", "consequence": "p.Asn84=",
        "class": "DUP_fixed", "frequency_in_carriers": 1.0,
        "source": "this work + per_read_genotyping" }
      // pos 21, 57, 189 omitted here — they are ORIG_polymorph, see entry
    ],
    "promoter_signature": {
      "window_relative_to_anchor": "-100..0",
      "anchor": "IGHG4_CH1_start",
      "canonical_motif": "CGGTTCTT...GTCTATCTGCGATGG",
      "duplicate_motif":  "GGGTTCTT...AACTGTCCGCGAGG",
      "note": "100 bp promoter / 5'UTR; NOT the Sγ4 switch region"
    },
    "switch_region": {
      "name": "Sgamma4",
      "discriminates_copies": false,
      "identity_between_copies_pct": 99.99,
      "note": "Cannot be used as discriminator; canonical and tandem-dup copies are essentially identical."
    },
    "ch2_ch3_drift": { "present": true, "per_sample_variable": true,
                       "evidence": "HG00329 hap1 G030290 vs G030293" }
  },

  // === Mapper strategy hints ===
  "mapping_strategy": {
    "primary": {
      "kmer_size": 25,
      "max_mismatch": 2,
      "include_flanking_bp": 200,                // ±200 bp around CH1 to capture promoter
      "require_unique_chain": false              // tandem dup → expect tied chains by design
    },
    "fallback_chain": [
      { "stage": 1, "name": "relaxed_mismatch", "kmer_size": 25, "max_mismatch": 4 },
      { "stage": 2, "name": "chain_only",       "use_extension": false },
      { "stage": 3, "name": "multi_position",   "report_all_top_k": true, "k": 4 },
      { "stage": 4, "name": "llm_checkpoint",   "opt_in_flag": "--llm-fallback" },
      { "stage": 5, "name": "novel_haplotype_flag", "emit_warning": true }
    ]
  },

  // === Provenance ===
  "provenance": {
    "primary_sources": [
      { "type": "preprint",    "id": "10.1101/2025.10.09.25337669v1", "year": 2025,
        "note": "medRxiv Belios et al. — IGHG4 ChimDup characterisation" }
    ],
    "supporting_databases": [
      "UCSC genomicSuperDups (hg38, frozen 2014)",
      "HPRC R2 pangenome (Liao 2023 + 2026 update)",
      "T2T-CHM13 v2.0 SD-BED (Vollger 2022)"
    ],
    "curator": "schlein-lab",
    "first_added": "2026-06-02",
    "last_reviewed": "2026-06-02"
  },

  // === Cross-references to other catalog entries ===
  "related_loci": ["IGHG_canondup_nahr_block", "IGH_constant_region"]
}
```

## 3. T2 (Bulk) JSONL Schema

One line per locus, one file per assembly per release:
`catalog/bulk/v2026.Q2.grch38.jsonl`, `.../v2026.Q2.chm13.jsonl`, etc.

```jsonc
{"locus_id":"ucsc_sd_grch38_0001","architecture":"tandem_dup_small",
 "coords":{"chrom":"chr1","start":12345,"end":34567},
 "identity_pct":94.5,"source":"UCSC_genomicSuperDups",
 "version":"v2026.Q2","promotion_score":2}
```

Required: `locus_id`, `architecture`, `coords`, `identity_pct`, `source`, `version`.
Optional: `mechanism[]`, `clinical_function[]`, `promotion_score`, `notes`.

## 4. Layout

```
catalog/
  schema/
    curated.schema.json        # JSON-Schema for T1
    bulk.schema.json           # JSON-Schema for T2
  curated/
    IGHG4_chimdup_tandem.json
    IGHG_canondup_nahr_block.json
    HLA_class_I_block.json
    ...
  bulk/
    v2026.Q2.grch38.jsonl
    v2026.Q2.chm13.jsonl
    v2026.Q2.mm39.jsonl
    ...
  index.json                   # top-level: lists releases, T1 entries, source-DB versions
  CHANGELOG.md
```

## 5. Release Process

1. **Q-1 month before release:** poll UCSC, T2T, gnomAD-SV, GIAB, HPRC. Diff against previous version.
2. **Q-2 weeks:** regenerate T2 JSONL files. Compute promotion-scores. Open PRs for score-≥5 candidates.
3. **Release day:** tag `catalog_v2026.Q<N>`. CI validates against JSON-Schema, runs LLmap smoke-tests for all T1 loci against a bundled mini-truth-set.
4. **Post-release:** publish CHANGELOG with new/promoted/deprecated locus list.

## 6. What this does NOT do

- Does **not** classify novel SegDup-character regions on-the-fly (T3 only emits warnings, no class).
- Does **not** include IPD-IMGT/HLA-derived sequence data in the catalog bundle (CC-BY-ND licence — external reference only).
- Does **not** override LLmap's general mapping logic outside catalog-covered loci.

## 7. Open questions

- T1-PSVs vs IGHG4↔IGHGP pseudocaller PSVs: keep entirely separate catalogs or cross-reference?
- Cross-species T2: BISER-self-computed for organisms without dedicated SD-track (mouse / rat / zebrafish / drosophila) — separate release cadence?
- LLM-fallback (stage 4): per-locus opt-in default, or global user toggle? Currently spec'd as user-toggle.
