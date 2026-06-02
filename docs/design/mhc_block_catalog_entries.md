# MHC-Block T1 Catalog Entries — Selection Rationale

**Status:** v2026.Q2 release notes — 2026-06-02
**Curator:** schlein-lab
**Schema:** `catalog/schema/curated.schema.json` (v0.2, two-axis)

This document justifies the six MHC-region T1 catalog entries added in v2026.Q2 and notes which additional MHC loci were *deliberately deferred* to a future backlog.

## 1. Selected T1 Entries (6)

| Entry | structural_architecture | haplotype_class | clinical_function | Reason for T1 promotion |
|---|---|---|---|---|
| `HLA_A_block` | `complex_mosaic` | `HLA_A_high_diversity` | `mhc_locus`, `mendelian_disease_locus` | Classical class-I anchor; >9 000 alleles; primary KIR ligand for HLA-A*03/11 supertype; canonical example for `gene_conversion_prone` |
| `HLA_B_block` | `complex_mosaic` | `HLA_B_high_diversity` | `mhc_locus`, `mendelian_disease_locus` | Most polymorphic human gene (>11 000 alleles); HLA-B*57:01 abacavir pharmacogenomic anchor; B*15:02 carbamazepine PGx |
| `HLA_C_block` | `complex_mosaic` | `HLA_C_high_diversity` | `mhc_locus`, `mendelian_disease_locus` | KIR ligand C1/C2 dimorphism; HLA-C*06:02 psoriasis; tight LD with HLA-B; ~88 % nt identity to HLA-B → severe cross-mapping |
| `HLA_DRB1_block` | `complex_mosaic` | `HLA_DRB1_high_diversity` | `mhc_locus`, `mendelian_disease_locus` | Class-II anchor; >4 000 alleles; covers DR-region CNV (DRB1/3/4/5) + pseudogene paralogs (DRB2/6/7/8/9); RA / T1D / MS / narcolepsy associations |
| `MICA_MICB_locus` | `complex_mosaic` | `MICA_polymorphic` | `mhc_locus`, `mendelian_disease_locus` | NKG2D-ligand paralog pair; ~84 % nt identity (genuine paralog cross-map); exon-5 microsatellite needs STR-aware handling; MICA-Del/null alleles |
| `TAP1_TAP2_locus` | `complex_mosaic` | `TAP1_TAP2_polymorphic` | `mhc_locus` (no `mendelian_disease_locus` — TAP-deficiency BLS-1 OMIM #604571 IS Mendelian, see §3) | Antigen-processing transporter pair; head-to-head paralogs; TAP2*Bky2 LoF allele; ABC-cassette cross-maps to ABCB-family elsewhere |

### Structural × Haplotype assignment summary

All six entries use `structural_architecture = complex_mosaic` per Norman 2017 / segdup_class_taxonomy.md §HLA-A/B/C entry. Rationale:

1. Pure HLA-A/B/C paralog pairs typically fall *below* the Bailey & Eichler 2002 ≥ 90 % / ≥ 1 kb SD threshold across full gene bodies — so they are not "SegDups" in the strict identity sense.
2. The classical MHC sub-region as a whole (chr6:28,510,120-33,480,577, ~5 Mb, per Horton 2008) IS a complex_mosaic of paralogous gene blocks formed by ancient duplication + extensive gene conversion + selection.
3. Inclusion in the LLmap SegDup catalog is *pragmatic* — extreme allelic diversity and paralog cross-mapping break single-reference mappers, regardless of strict SD-definition compliance.

`haplotype_class` is **deliberately a free-form per-locus tag** (schema field is `^[A-Za-z][A-Za-z0-9_]+$`, no global enum). This matches the IGHG-entry convention from v2026.Q1 — each MHC locus has its own combinatorics of alleles, so any global enum would either bloat or lose information. See spec §1.2 (axis split rationale).

All six entries share `mechanism = ["gene_conversion_prone"]` (Parham 1989; Chen 2007; Erlich 2012 reviews) and `nahr_status = "non_nahr"` — MHC diversification is gene-conversion + balancing-selection driven, NOT NAHR-driven. This is structurally distinct from the IGHG NAHR-block entry (`nahr_status = "ambiguous_low_identity"`).

## 2. Mapping Strategy

All six entries use a common MHC-specific strategy that diverges from IGHG4 in three points:

| Parameter | IGHG4 entries | MHC entries | Rationale |
|---|---|---|---|
| `primary.max_mismatch` | 2 | **4** | Per-position diversity in HLA exons 2/3 is ~10× the genome average; mm=2 misses too many true hits |
| `primary.require_unique_chain` | true (single-copy) / false (tandem) | **false** (all) | HLA/MIC/TAP have legitimate paralog cross-maps; uniqueness requirement would mask them |
| `fallback_chain[0]` relaxed mismatch | 4 | **6** | Higher headroom for extreme allelic diversity beyond standard |

DRB1 additionally uses `multi_position.k = 8` (vs 4-6 elsewhere) because the DRB-paralog group (DRB1 + DRB3/4/5 functional + DRB2/6/7/8/9 pseudogene) routinely produces ≥ 8 plausible hits.

`include_flanking_bp = 200` and `include_flanking_anchor = "<gene>_exon2_start"` are uniform — exon 2 is the hyperpolymorphic anchor for all class-I/II genes.

## 3. Sequence/PSV License — IPD-IMGT/HLA

**IPD-IMGT/HLA is CC-BY-ND.** No allele-specific PSVs, no reference sequences, and no PSV-frequency tables from IPD-IMGT have been embedded in the catalog. Each entry:

- Lists `diagnostic_features.diversity_hotspots` as *structural* descriptions (exon 2 = α1 domain peptide groove, etc.) with no inline sequence.
- Carries `diagnostic_features.discriminating_snps = []` (empty array) — schema-valid, signals "see IPD-IMGT for PSVs" rather than fabricated content.
- Provides `diagnostic_features.ipd_imgt_reference` block with the URL and license note, and recommends runtime API resolution.
- Adds `provenance.license_note` reiterating that IPD-IMGT sequences/PSVs are NOT embedded.

This satisfies CC-BY-ND (linking and citing are permitted; deriving / redistributing modified copies is not).

### Note on `clinical_function`

`mhc_locus` is applied to all six entries. `mendelian_disease_locus` is added to HLA-A/B/C/DRB1 + MICA/MICB because they have direct OMIM Mendelian associations (e.g. ankylosing spondylitis HLA-B*27, narcolepsy DRB1*15:01-DQB1*06:02, MICA in psoriatic arthritis). TAP1/TAP2 was a borderline call — BLS-1 (OMIM #604571) is a Mendelian recessive transporter deficiency, but most TAP polymorphism work is LD-driven via class-II haplotypes. Decision: do **not** add `mendelian_disease_locus` to TAP1/TAP2 in v2026.Q2 because the polymorphic TAP alleles themselves are not the primary disease driver — the BLS-1 phenotype involves rare LoF rather than common polymorphism. Revisit if a TAP-LoF-specific entry is needed.

## 4. Backlog — Deliberately Deferred from v2026.Q2

The following MHC loci were considered but **not** promoted to T1 in this release. They should be candidates for v2026.Q3 or later, with explicit evidence triggers.

### 4.1 Class-II beyond DRB1

| Locus | Reason for deferral | Promotion trigger |
|---|---|---|
| `HLA_DPB1` | Less clinically anchored than DRB1; nonetheless ~2 000+ alleles and HSCT-mismatch relevance. Deferred to keep v2026.Q2 reviewable. | Add when v2026.Q3 includes HSCT-typing use case OR when DRB1 entry exposes shared infrastructure. |
| `HLA_DPA1` | Pairs with DPB1; less polymorphic (~80 alleles). | Add jointly with DPB1 as `HLA_DP_paralog_pair` (paired entry similar to MICA/MICB pattern). |
| `HLA_DQB1` | ~2 200 alleles; T1D / celiac / narcolepsy core. Genuinely belongs in T1. | **Strong v2026.Q3 candidate** — listed first in backlog. |
| `HLA_DQA1` | Pairs with DQB1; ~100 alleles. | Joint entry with DQB1, similar to DPA1-DPB1. |
| `HLA_DRB3` / `DRB4` / `DRB5` | Functional paralogs of DRB1 only present on some haplotypes. Currently covered as "co_duplicated_genes" in `HLA_DRB1_block`. | Separate entries only if the joint DRB1-block treatment proves insufficient. |

### 4.2 Class-I beyond A/B/C

| Locus | Reason for deferral | Promotion trigger |
|---|---|---|
| `HLA_E` | Non-classical class-I, NKG2A ligand; only ~300 alleles. Less mapping-confusion. | Add if NKG2A-pharmacology use case emerges. |
| `HLA_F` / `HLA_G` | Non-classical, oligomorphic. Limited polymorphism makes mapping tractable. | Likely T2-bulk-only forever. |

### 4.3 Class-III and other MHC genes

| Locus | Reason for deferral | Promotion trigger |
|---|---|---|
| `C4A` / `C4B` | Complement C4 has tandem-dup CNV (1-4 copies) + HERV-K insertion polymorphism — this IS a real SegDup case and a strong T1 candidate. Deferred only because it needs its own pangenome characterization which is not in v2026.Q2. | **v2026.Q3 priority** — separate entry with `tandem_dup_small` + `nahr_prone` likely. |
| `CYP21A2` / `CYP21A1P` | Classical pseudogene-paralog pair, CAH disease locus; pure `pseudogene_paralog_pair` case rather than MHC-mosaic. | Add as separate `CYP21A2_CYP21A1P_pair` entry. |
| `PSMB8` / `PSMB9` | Co-located with TAP1/TAP2 (`co_duplicated_genes` field captures them). Polymorphism is modest. | Add only if proteasome-allele use case justifies. |
| `BTNL2` | Sarcoidosis, IBD associations. Modest polymorphism. | Promote only on score ≥ 5 (OMIM + AF). |

## 5. Cross-references

All six entries reference each other (and IGHG entries are not cross-referenced because they are mechanistically distinct — IGHG is NAHR-block / tandem-dup architecture, MHC is gene-conversion / complex-mosaic). Conceptual relationships are captured in `related_loci`:

```
HLA_A_block        ← → HLA_B_block / HLA_C_block / HLA_DRB1_block / MICA_MICB_locus / TAP1_TAP2_locus
HLA_B_block        ← → (same)
HLA_C_block        ← → (same)
HLA_DRB1_block     ← → (same)
MICA_MICB_locus    ← → (same)
TAP1_TAP2_locus    ← → (same)
```

This is *intra-MHC* clique cross-referencing. An MHC↔IGH cross-reference would be misleading because the mechanisms are fundamentally different (no NAHR in MHC, no Bailey-SD in classical HLA-A/B/C).

## 6. Validation

```
$ python3 catalog/validate.py
✓ catalog/curated/HLA_A_block.json: OK
✓ catalog/curated/HLA_B_block.json: OK
✓ catalog/curated/HLA_C_block.json: OK
✓ catalog/curated/HLA_DRB1_block.json: OK
✓ catalog/curated/IGHG4_chimdup_canonical_arch.json: OK
✓ catalog/curated/IGHG4_chimdup_tandem.json: OK
✓ catalog/curated/IGHG_canondup_nahr_block.json: OK
✓ catalog/curated/MICA_MICB_locus.json: OK
✓ catalog/curated/NBPF_olduvai_cluster.json: OK
✓ catalog/curated/TAP1_TAP2_locus.json: OK

10/10 valid; 0 failed.
```

10/10 valid (6 new MHC + 3 IGHG + 1 NBPF in-flight from parallel work). Brief stated 9/9 (3 IGHG + 6 new); the NBPF entry appeared during this work cycle as an independent v2026.Q2 addition and validates cleanly.

## 7. References

- **Norman PJ et al. (2017)** *Sequences of 95 human MHC haplotypes reveal extreme coding variation in genes other than highly polymorphic HLA class I and II.* Genome Res. 27(5):813-823. doi:10.1101/gr.213538.116
- **Horton R et al. (2008)** *Variation analysis and gene annotation of eight MHC haplotypes: The MHC Haplotype Project.* Immunogenetics. doi:10.1007/s00251-007-0262-2
- **Robinson J et al.** IPD-IMGT/HLA Database. https://www.ebi.ac.uk/ipd/imgt/hla/ (CC-BY-ND, external reference only)
- **Parham P, Lawlor DA (1989)** Diversity and diversification of HLA-A,B,C alleles. *J Immunol.* 142(11):3937-3950.
- **Chen JM, Cooper DN, Chuzhanova N, Férec C, Patrinos GP (2007)** Gene conversion: mechanisms, evolution and human disease. *Nat Rev Genet.* doi:10.1038/nrg2193
- **Bahram S et al. (1996)** MIC, a class I gene family of major histocompatibility complex molecules. *Immunol Rev.*
- **Erlich HA (2012)** HLA DNA typing: past, present, and future. *Tissue Antigens.*
- **OMIM #604571** Bare lymphocyte syndrome, type I (TAP1/TAP2 deficiency).
