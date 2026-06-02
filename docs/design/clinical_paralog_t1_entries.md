# Clinical Paralog T1 Catalog Entries

**Status:** initial release, v2026.Q2
**Author:** schlein-lab
**Date:** 2026-06-02
**Companion to:** [segdup_catalog_spec.md](segdup_catalog_spec.md), [segdup_class_taxonomy.md](segdup_class_taxonomy.md)

This document records the rationale and the empirical caveats for the five
clinical-paralog T1 catalog entries added in v2026.Q2:

1. `NBPF_olduvai_cluster` — chr1 Olduvai/DUF1220 family (~30 paralogs)
2. `NPHP1_locus` — 2q13 tandem locus, recurrent ~290 kb deletion
3. `FCGR_cluster` — 1q23.3 low-affinity Fcγ-receptor cluster
4. `NCF1_cluster` — 7q11.23 NCF1 / NCF1B / NCF1C (CGD)
5. `PMS2_PMS2P3_pair` — 7p22 PMS2 / PMS2CL (Lynch syndrome)

They complement the three IGHG entries (`IGHG4_chimdup_tandem`,
`IGHG4_chimdup_canonical_arch`, `IGHG_canondup_nahr_block`) and the MHC /
Y-palindrome blocks added in the same release.

---

## 1. Why these five

These are the highest-priority *non-IGH* clinical paralog loci called out
in the SegDup taxonomy doc (§6 — the "eight LLmap anchor loci"). They cover
four distinct mechanistic regimes:

| Locus | Mechanism regime | Why it earns a T1 slot |
|---|---|---|
| NBPF | Core-duplicon expansion | Canonical high-copy paralog family; primate-evolution literature anchor (Jiang 2007); largest human-lineage CNV |
| NPHP1 | Classical Lupski-style DP-LCR NAHR | Recurrent ~290 kb deletion in ~80% of NPHP1 disease cases (Saunier 2000); textbook NAHR substrate |
| FCGR | Tandem NAHR + IGC + chimerism | Clinically actionable CNV (FCGR3B → autoimmunity); good PSV literature (Niederer 2010); chimeric NAHR products characterised (Mueller 2015) |
| NCF1 | Gene-conversion → disease | Prototype gene-conversion-disease locus; >90% of NCF1-CGD alleles are pseudogene-converted (Vázquez 2001) |
| PMS2 / PMS2CL | Ongoing IGC + clinical reporting trap | Hayward 2007 showed PSVs are "unsafe"; Lynch-syndrome reporting hotspot; the *hardest* locus for short-read calling |

---

## 2. Architecture / mechanism / clinical tags chosen

| Locus | structural_architecture | mechanism | clinical_function | nahr_status |
|---|---|---|---|---|
| NBPF_olduvai_cluster | `intrachromosomal_interspersed` | `nahr_prone`, `core_duplicon` | `dispersed_paralog_gene_family`, `pangenome_only` | `high_identity_nahr_prone` |
| NPHP1_locus | `tandem_dup_large` | `nahr_prone` | `mendelian_disease_locus`, `microdel_microdup_syndrome` | `high_identity_nahr_prone` |
| FCGR_cluster | `tandem_dup_large` | `nahr_prone`, `gene_conversion_prone` | `dispersed_paralog_gene_family`, `immune_receptor_locus` | `high_identity_nahr_prone` |
| NCF1_cluster | `tandem_dup_small` | `gene_conversion_prone` | `mendelian_disease_locus`, `pseudogene_paralog_pair` | `high_identity_nahr_prone` |
| PMS2_PMS2P3_pair | `intrachromosomal_interspersed` | `gene_conversion_prone` | `mendelian_disease_locus`, `pseudogene_paralog_pair`, `gene_conversion_disease_locus` | `ambiguous_low_identity` |

Notes:
- NPHP1 architecture is set to `tandem_dup_large` referring to the 45 kb
  flanking LCRs (the NAHR substrate). The taxonomy doc §6 also discusses a
  `tandem_dup_small` reading for the same 45 kb LCRs — we have chosen the
  `large` bin because the 50 kb threshold puts 45 kb right at the boundary
  and the recurrent-deletion product (~290 kb) is firmly in the
  large-segdup-block regime.
- FCGR is tagged `tandem_dup_large` (82.5 kb unit) per Niederer 2010 /
  Breunis 2009.
- NCF1 is tagged `tandem_dup_small` (~15 kb unit) per Heyworth 2002 — but
  embedded in a much larger WBS LCR mosaic.
- PMS2/PMS2CL nahr_status is `ambiguous_low_identity` because the two
  paralogs are in *inverted* orientation and the dominant mechanism is IGC
  not NAHR DEL/DUP.

---

## 3. Coordinate sources (GRCh38)

All coordinates retrieved 2026-06-02 from Ensembl REST API
(`rest.ensembl.org/lookup/symbol/homo_sapiens/<GENE>`). Assembly = GRCh38.

| Locus | Region | Coordinate source |
|---|---|---|
| NBPF | chr1:16.56-149.8 Mb (bounding-box, not contiguous) | NBPF1 anchor from Ensembl |
| NPHP1 | chr2:110,122,311-110,205,066 (-strand) | Ensembl |
| FCGR cluster | chr1:161,505,415-161,710,000 | FCGR2A, FCGR2B, FCGR3B from Ensembl |
| NCF1 cluster | chr7:73,220,646-75,172,044 | NCF1, NCF1B, NCF1C from Ensembl |
| PMS2 / PMS2CL | chr7:5,970,925-6,751,392 | PMS2 + PMS2CL from Ensembl |

⚠ **PMS2 naming caveat (catalogued):** Ensembl uses `PMS2P3` for a *different*
pseudogene at chr7:75.5 Mb. The clinically relevant PMS2 pseudogene
(Hayward 2007's "PMS2CL") is at chr7:6.73 Mb. LLmap therefore hard-codes
coordinates rather than relying on gene-symbol lookup. The file is named
`PMS2_PMS2P3_pair.json` to match the user task wording; the JSON body
documents the disambiguation.

---

## 4. Validation status

```
$ cd /home/<user>/llmap-local && python3 catalog/validate.py
✓ catalog/curated/FCGR_cluster.json: OK
✓ catalog/curated/NBPF_olduvai_cluster.json: OK
✓ catalog/curated/NCF1_cluster.json: OK
✓ catalog/curated/NPHP1_locus.json: OK
✓ catalog/curated/PMS2_PMS2P3_pair.json: OK
... (plus the 14 pre-existing entries)
19/19 valid; 0 failed.
```

---

## 5. Cross-references between entries

- `NBPF_olduvai_cluster.related_loci` → `NPHP1_locus` (both are
  archetypal SD-clusters of opposite mechanisms: core-duplicon expansion
  vs Lupski-NAHR recurrent-DEL).
- `NPHP1_locus.related_loci` → `NBPF_olduvai_cluster`
- `FCGR_cluster.related_loci` → `NCF1_cluster` (both 7q11/1q23 region
  pseudogene-flanked immune-effector clusters)
- `NCF1_cluster.related_loci` → `FCGR_cluster`, `PMS2_PMS2P3_pair`
  (gene-conversion-disease cluster)
- `PMS2_PMS2P3_pair.related_loci` → `NCF1_cluster` (the two
  gene-conversion-disease locus entries cross-reference each other)
- The three IGHG entries already cross-reference each other; we did not
  add a cross-link from IGHG to FCGR or NCF1 even though all three are
  immune-relevant SD clusters — those would be functional, not
  mechanistic, links and we follow the catalog convention of
  mechanism-based cross-references.

---

## 6. Discriminating-SNP / PSV evidence per locus

### NBPF
**PSV table:** not provided. Per-paralog point-SNP discrimination is not
feasible at the catalog level because intra-cluster identity exceeds 99%
over the Olduvai HOR and paralog correspondence is not stable across
haplotypes (250-350 Olduvai domains per diploid, vs ~30 named NBPF genes
in GRCh38). The entry instead specifies a *HOR-pattern* discriminator
strategy (Olduvai triplet composition + flanking unique anchor).
**Honesty flag:** This is one of the loci where PSVs are genuinely hard to
find — the literature does not catalog per-paralog SNPs in usable tables.
Vollger 2022 + HPRC R2 + T2T-CHM13 are the only assemblies that resolve
the cluster.

### NPHP1
**PSV table:** not provided. The locus is diagnosed structurally
(NPHP1 copy-number 0/1/2 + flanking 45 kb LCR junction signature), not by
SNP. Saunier 2000 reports breakpoint mapping but no SNP-level PSV table
useful for read-anchored mapping. Entry uses depth + junction-search
strategy.

### FCGR
**PSV table:** ✓ available. Niederer 2010 (doi:10.1093/hmg/ddq216)
supplementary contains the canonical PSV table. We list three
representative coding-region positions
(rs1801274 FCGR2A H131R, rs396991 FCGR3A V158F, FCGR3B NA1/NA2 block)
as illustrative entries; the full Niederer supplementary should be
ingested by the mapper as a downstream resource.

### NCF1
**PSV table:** primary diagnostic = the GTGT-vs-GT signature at exon 2
start (Vázquez 2001, Heyworth 2002). Additional positions across the
15 kb unit are catalogued in Hoffmeister 2025 Frontiers
supplementary (doi:10.3389/fimmu.2025.1640496). We list the GTGT-vs-GT
as the canonical PSV and cite the supplementary for the rest.

### PMS2 / PMS2CL
**PSV table:** ✓ Hayward 2007 supplementary BUT explicitly flagged
*unsafe* by the same paper. We list exon 10 (absent from pseudogene) as
the only structurally-fixed anchor, and the 3'-exon PSV block as
`class: uncertain` with a strong note that LLmap should NOT auto-call
variants there. This is the most clinically delicate entry — the
catalog actively documents what NOT to call confidently.

---

## 7. Mapping-strategy choices (k-mer, mismatch)

| Locus | kmer_size | max_mismatch | Rationale |
|---|---|---|---|
| NBPF | 25 | 4 | >99% intra-cluster identity; ~3 SNPs/window expected; multi-mapping is correct |
| NPHP1 | 25 | 2 | 99.5% LCR identity; mapping is depth-driven not SNP-driven |
| FCGR | 25 | 3 | 98.5% identity; 82.5 kb unit; needs more tolerance than IGHG4 |
| NCF1 | 25 | 2 | 98% identity; flanking unique anchors close by (GTF2IRD2 etc.) |
| PMS2 | 25 | 2 | 98% identity over 100 kb; flag rather than force in 3' exons |

All five use the standard 5-stage `fallback_chain`
(relaxed_mismatch → chain_only / locus-specific anchor → multi_position →
llm_checkpoint → novel_haplotype_flag), customized where the locus
biology demands it (NPHP1 adds `depth_check` and `junction_search`; PMS2
adds `exon10_anchor_walk`; NCF1 adds `anchor_walk` from GTF2IRD2 flanks;
FCGR adds `junction_search_chimera`).

---

## 8. Loci where PSV lists were hard to find

Three of the five entries lacked a clean, reusable PSV table:

1. **NBPF** — high copy number plus haplotype-unstable paralog
   correspondence means no static PSV table exists. The right primitive
   is a HOR-pattern matcher, not a SNP set.
2. **NPHP1** — diagnosis is structural (deletion), not SNP-level. The
   literature does not publish a flanking-LCR PSV table because it is
   not needed for clinical detection.
3. **PMS2 / PMS2CL** — Hayward 2007 *does* publish a PSV table, but the
   same paper *explicitly* warns it is unsafe (IGC drift). This is the
   one locus where the PSV list exists but cannot be used naively.

FCGR (Niederer 2010 supplementary Table S2) and NCF1 (Vázquez 2001 +
Hoffmeister 2025) have usable PSV tables.

---

## 9. References (DOIs added in this batch)

- Sikela JM et al. *Changing the name of the NBPF/DUF1220 domain to the Olduvai domain.* F1000Research. doi:[10.12688/f1000research.13586.2](https://doi.org/10.12688/f1000research.13586.2)
- Jiang Z et al. *Ancestral reconstruction of segmental duplications reveals punctuated cores of human genome evolution.* Nat Genet. 2007. doi:[10.1038/ng.2007.9](https://doi.org/10.1038/ng.2007.9)
- Vollger MR et al. *Segmental duplications and their variation in a complete human genome.* Science. 2022. doi:[10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)
- Saunier S et al. *Characterization of the NPHP1 locus: mutational mechanism involved in deletions in familial juvenile nephronophthisis.* Am J Hum Genet. 2000. doi:[10.1086/302819](https://doi.org/10.1086/302819)
- Yu H et al. *Comparative genomic analyses of the human NPHP1 locus reveal complex genomic architecture and its regional evolution in primates.* PLOS Genet. 2015. doi:[10.1371/journal.pgen.1005686](https://doi.org/10.1371/journal.pgen.1005686)
- OMIM 256100 (NPHP1 disease) and OMIM 607100 (NPHP1 gene)
- Niederer HA et al. *Copy number, linkage disequilibrium and disease association in the FCGR locus.* Hum Mol Genet. 2010. doi:[10.1093/hmg/ddq216](https://doi.org/10.1093/hmg/ddq216)
- Breunis WB et al. *Copy number variation at the FCGR locus includes FCGR3A, FCGR2C and FCGR3B but not FCGR2A and FCGR2B.* Hum Mutat. 2009. doi:[10.1002/humu.20997](https://doi.org/10.1002/humu.20997)
- Mueller M et al. *Nonallelic homologous recombination of the FCGR2/3 locus results in copy number variation and novel chimeric FCGR2 genes with aberrant functional expression.* Genes Immun. 2015. doi:[10.1038/gene.2015.25](https://doi.org/10.1038/gene.2015.25)
- Vázquez N et al. *Mutational analysis of patients with p47-phox-deficient chronic granulomatous disease: significance of recombination events between the p47-phox gene (NCF1) and its highly homologous pseudogenes.* Exp Hematol. 2001. doi:[10.1016/s0301-472x(00)00646-9](https://doi.org/10.1016/s0301-472x(00)00646-9)
- Heyworth PG, Noack D, Cross AR. *Identification of a novel NCF-1 (p47-phox) pseudogene not containing the signature GT deletion: significance for A47° chronic granulomatous disease carrier detection.* Blood. 2002. doi:[10.1182/blood-2002-03-0861](https://doi.org/10.1182/blood-2002-03-0861)
- Hoffmeister F et al. *Reliable genetic diagnosis of NCF1 (p47phox)-deficient chronic granulomatous disease using high-throughput sequencing.* Front Immunol. 2025. doi:[10.3389/fimmu.2025.1640496](https://doi.org/10.3389/fimmu.2025.1640496)
- Hayward BE et al. *Extensive gene conversion at the PMS2 DNA mismatch repair locus.* Hum Mutat. 2007. doi:[10.1002/humu.20457](https://doi.org/10.1002/humu.20457)
- Bouras A et al. *PMS2 or PMS2CL? Characterization of variants detected in the 3′ of the PMS2 gene.* Genes Chromosomes Cancer. 2024. doi:[10.1002/gcc.23193](https://doi.org/10.1002/gcc.23193)

---

## 10. Caveats and known gaps

1. **GRCh38 vs CHM13 vs HPRC R2 coordinate reconciliation** for NBPF is
   incomplete in this release. The Olduvai cluster is one of the loci where
   GRCh38 collapses several HOR copies; T2T-CHM13 + HPRC R2 should be
   added as `coords` entries in a future release once Vollger 2022 /
   HPRC R2 NBPF coordinate tables have been parsed.
2. **PMS2 Ensembl naming bug** (PMS2P3 at Ensembl chr7:75.5 Mb is a
   different pseudogene) is documented in the entry's `coords.note` but
   should be cross-flagged in the catalog README so users don't blindly
   trust the Ensembl record.
3. **NPHP1 test cases** are placeholders (no real patient samples in the
   schlein-lab dataset for this locus); update when real CGD/NPHP/Lynch
   patient samples enter the catalog test set.
4. **FCGR Mueller-2013 citation** in the original task was actually
   Mueller et al. 2015 *Genes Immun* (the matching FCGR-NAHR-chimera
   paper); fixed in the JSON `provenance.primary_sources` note.
5. **Vázquez 2001 PNAS** in the original task was actually Vázquez et al.
   2001 *Exp Hematol*; fixed.
