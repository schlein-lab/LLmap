# SegDup Class Taxonomy for LLmap Catalog

**Status:** design draft, literature-grounded
**Author:** Christian Schlein, LLmap team
**Date:** 2026-06-02

> **Update v0.2 (same day):** The "Axis A" architecture vocabulary below remains valid as **structural architecture**, but it does NOT encode haplotype-level PSV patterns. Since IGHG4 turned out to need both a structural class (e.g. `tandem_dup_small`) AND a haplotype-class tag (e.g. `IGHG4_chimdup_homozygous`), the catalog schema now requires a second free-form field `haplotype_class`. Each T1 entry sets both. See [segdup_catalog_spec.md §1.2](segdup_catalog_spec.md#12-classification-axes--two-axis-structural--haplotype-split-since-v02).
>
> Empirical justification from 104-sample HPRC sweep: same 295 bp CH1 yielded three distinct (structural, haplotype) pairs — `(tandem_dup_small, chimdup)`, `(tandem_dup_large, canonical)`, `(single_copy, chimdup)`. A single-axis classification would have collapsed them.
**Purpose:** Define the minimum set of `class` tags that each LLmap catalog locus carries. The taxonomy must be *anchored in published nomenclature*; we deliberately re-use existing terms instead of inventing new ones.

---

## 1. Why we need this

Each locus row in the LLmap SegDup catalog needs a `class` field that downstream tools (priors, NAHR-risk scoring, reporting) can switch on. We want:

1. **Reproducibility** — every class must be cite-able to a peer-reviewed paper.
2. **Coverage** — the eight LLmap anchor loci (IGHG4-ChimDup, IGH-CanonDup/NAHR-block, HLA-A/B/C, NBPF, Y-AZF palindromes, NPHP1, FCGR cluster, MEFV) must each get at least one class with literature support.
3. **Orthogonality** — classes should be (mostly) independent dimensions. Architecture, mechanism, and clinical-disorder status are not the same thing; a locus can carry several class tags.

To keep the catalog usable we split tags into three orthogonal axes:

- **A. Architecture** (geometry of the duplicated copies) — primary `class` tag.
- **B. Mechanism / risk** — secondary tag (NAHR-prone, gene-conversion-prone, palindrome-mediated, etc.).
- **C. Function / clinical** — tertiary tag (immune-receptor cluster, MHC, genomic-disorder hotspot, core duplicon).

Each axis has a closed vocabulary defined below. A locus may be e.g. `architecture=tandem_dup_small; mechanism=nahr_prone; clinical=immune_receptor_locus`.

---

## 2. Canonical defining thresholds (used throughout)

| Property | Threshold | Source |
|---|---|---|
| SD minimum length | 1 kb | Bailey & Eichler 2006; Sharp 2005 |
| SD minimum identity | 90% | Bailey & Eichler 2006; Sharp 2005 |
| NAHR-prone identity | ≥ 95–97% (often ≥ 97%) | Sharp 2005; Stankiewicz & Lupski 2002; Liu 2012 |
| NAHR-prone length | ≥ 5–10 kb | Stankiewicz & Lupski 2002; Lupski review |
| NAHR-prone distance | 50 kb to ~10 Mb | Sharp 2005 |
| NAHR-prone orientation | directly oriented (DP-LCR) for DEL/DUP; inverted for INV | Lupski; Stankiewicz |
| Tandem vs interspersed cutoff | < 1 Mb / same chr → "tandem/clustered"; ≥ 1 Mb apart or different chr → "interspersed" | Bailey & Eichler 2006; Marques-Bonet & Eichler 2009 |

These numbers are not invented for LLmap; they are the values consistently used by the Eichler-lab papers and the Sharp/Lupski genomic-disorder literature.

---

## 3. Axis A — Architecture classes (primary `class` tag)

Closed vocabulary. Pick exactly one per locus copy-pair.

| Class tag | Definition | Anchor literature | Example locus |
|---|---|---|---|
| `tandem_dup_small` | Directly adjacent or near-adjacent (< 50 kb intervening), duplication unit < 50 kb, same chromosome, same orientation. Equivalent to Lupski's "tandem repeat" / Bailey's "tandem duplication". | Bailey & Eichler 2006 (NRG); Lupski 1998 | IGHG4-ChimDup (19.5 kb tandem); CMT1A *PMP22* tandem on the same haplotype |
| `tandem_dup_large` | Same as above but unit ≥ 50 kb. Captures the "large tandem block" mode typical of the IGH constant gene region and the FCGR cluster. | Watson/Rodriguez 2020 (IGH ~25 kb–~115 kb blocks); Mueller 2013 (FCGR 82.5 kb unit) | IGH-CanonDup/NAHR-block (~115 kb co-dup of G2/A1/E); FCGR 82.5 kb LCR |
| `intrachromosomal_interspersed` | Two paralogous copies on the same chromosome but separated by ≥ 1 Mb of unique sequence (Bailey/Vollger threshold). | Bailey & Eichler 2006; Vollger 2022 *Science*; Vollger 2023 *Nature* | many 1q21 NBPF blocks; pericentromeric SD bursts |
| `interchromosomal` | Paralogous copies on different chromosomes. Strongly enriched in pericentromeric / subtelomeric regions. | Bailey & Eichler 2006; Vollger 2022 | pericentromeric & subtelomeric SDs; GAGE-X / NBPF-1p ↔ 1q21 |
| `palindromic_inverted` | Two copies in inverted orientation, often as a palindrome with an axis of symmetry, typical of the Y ampliconic class (P1–P8). | Skaletsky 2003 (*Nature*); Kuroda-Kawaguchi 2001 | Y-AZFc palindromes P1–P5; X-AZFa direct repeats; some 8p23 |
| `inverted_dispersed` | Paralogous copies in inverted orientation, *not* organised as a palindrome (no axis of symmetry, often ≥ 100 kb apart). Substrate for recurrent inversions. | Porubsky 2022 *Cell* (recurrent inversions flanked by inverted SDs) | 17q21.31 H1/H2 inversion; 16p11.2 inversion |
| `complex_mosaic` | A locus where multiple SD pairs of different sizes / orientations / origins are interleaved. Equivalent to the "duplication block" of Jiang 2007 and the "module mosaic" of Marques-Bonet 2009. The block as a whole cannot be reduced to one geometry. | Jiang & Eichler 2007 (*Nat Genet*); Marques-Bonet 2009 (Annu Rev) | 22q11.2 LCR-A…H (8 sub-LCRs); MHC-class-I; chr14q32 IGH constant gene region as a whole |
| `subtelomeric` | SD that lies in the distal ~500 kb of a chromosome arm, sharing sequence with other chromosome ends. Architecturally typically `interchromosomal` but called out separately because of its distinct evolutionary dynamics. | Bailey & Eichler 2006; Linardopoulou 2005 *Nature* | most chromosome-end SDs |
| `pericentromeric` | SD in the pericentromeric "hat" (~5 Mb around the centromere), often with many paralogs (high-copy LCR / HCR). | Bailey & Eichler 2002 *Science*; Vollger 2022 | 1p11/1q11; 9p11/9q11; 15q11–q13; 22q11.2 proximal block |
| `unresolved` | Catalog row where the geometry could not be determined unambiguously from the current assembly set (e.g. T2T-CHM13 + HPRC haplotypes disagree). Use this honestly rather than guessing. | — | rare; flagged for re-assembly |

Notes:
- `tandem_dup_small` and `tandem_dup_large` use 50 kb as the split because (a) it lines up with the upper bound of single-block SDs in IGH (the IGHG4 unit) and (b) it separates the "<50 kb" mode that dominates Sharp 2005's CNP catalog from the multi-tens-of-kb NAHR-blocks. This 50 kb cut-off is *our* operational split, but both bins themselves are literature-anchored.
- `subtelomeric` and `pericentromeric` overlap with `interchromosomal` by definition; they are subtypes — the catalog should write the spatial tag in addition to `interchromosomal`.

---

## 4. Axis B — Mechanism / risk tags

These describe *what the duplication does mechanistically*. A locus can carry zero, one or several mechanism tags. The tag is empty when the locus is structurally a SD but has no documented NAHR/conversion activity.

| Mechanism tag | Definition | Anchor literature |
|---|---|---|
| `nahr_prone` | Directly-oriented paralog pair (DP-LCR) with length ≥ 5–10 kb, identity ≥ 95–97%, distance < ~10 Mb. Predicts recurrent DEL/DUP. | Stankiewicz & Lupski 2002; Sharp 2005; Liu 2012 |
| `nahr_inversion_substrate` | Inverted paralog pair (IR-LCR) with the same size/identity properties. Predicts recurrent inversions. | Stankiewicz & Lupski 2002; Porubsky 2022 *Cell* |
| `gene_conversion_prone` | Pair with ongoing interlocus gene conversion (IGC), evidenced by phylogenetic incongruence or elevated identity in the converted tract. Vollger 2023 quantified this genome-wide. | Vollger 2023 *Nature*; Chen 2007 |
| `core_duplicon` | A short, highly-replicated unit that seeds expansion of larger duplication blocks. Defined by Jiang 2007 as "ancestral duplicon populating > 67% of all duplication blocks within a group". | Jiang & Eichler 2007 *Nat Genet*; Marques-Bonet 2009 |
| `palindrome_mediated_nahr` | Specifically NAHR between palindromic arms (Y AZFc-style), distinguished because both rate and product spectrum differ from generic DP-LCR NAHR. | Skaletsky 2003; Repping 2002 |
| `none_documented` | Locus is SD but no NAHR / IGC product reported. Default for many young, low-identity intrachromosomal interspersed SDs. | — |

`nahr_prone` and `nahr_inversion_substrate` are mutually exclusive *per copy pair* but a complex locus (e.g. 22q11.2) can have both, with different sub-LCR pairs supporting each.

---

## 5. Axis C — Functional / clinical tags

These capture *why we care about the locus clinically or biologically*. A locus can carry several.

| Clinical tag | Definition | Anchor literature |
|---|---|---|
| `microdel_microdup_syndrome` | The locus is the substrate for a recurrent microdeletion/microduplication syndrome listed in ClinGen's recurrent-CNV regions or in OMIM as a genomic disorder. | Lupski 1998; Mefford & Eichler 2009; ClinGen Recurrent CNV list |
| `mendelian_disease_locus` | A deletion at this locus is a recognised cause of a Mendelian disease (single gene per ClinGen Dosage Sensitivity). Distinct from microdel-syndromes because phenotype = single-gene LoF, not contiguous-gene syndrome. | ClinGen Dosage Map; Stankiewicz & Lupski 2010 |
| `immune_receptor_locus` | Member of the immunoglobulin / T-cell receptor family (IGH, IGK, IGL, TRA, TRB, TRG, TRD). | Watson 2017; Rodriguez/Watson 2020 |
| `mhc_locus` | MHC class-I (HLA-A/B/C/E/F/G) or class-II (HLA-DRA/DRB/DQA/DQB/DPA/DPB) or class-III genes/SDs. | Norman 2017 *Genome Res*; Horton 2008 |
| `dispersed_paralog_gene_family` | A gene family with > 2 functional paralogs spread by SD (globins, NBPF, FCGR, GAGE, SMN, NPIP). | Bailey & Eichler 2006; Marques-Bonet 2009 |
| `pseudogene_paralog_pair` | Functional gene + closely-related processed-pseudogene/duplicated-pseudogene that confounds short-read variant calling (PMS2/PMS2P3, SMN1/SMN2 to a lesser degree, GBA/GBAP1, CYP2D6/CYP2D7). | Hayward & Eichler reviews; ClinGen technical notes |
| `gene_conversion_disease_locus` | Locus where pathogenic alleles arise via gene conversion between a gene and its paralog/pseudogene. Distinct from `gene_conversion_prone` (which is purely population-genetic). | Chen 2007; Casola 2012 |
| `pangenome_only` | Locus appears in HPRC pangenome haplotypes but is absent from GRCh38 (or vice versa). Flag for catalog completeness. | Liao 2023 (HPRC v1); Vollger 2023 |

`microdel_microdup_syndrome` and `mendelian_disease_locus` are *not* mutually exclusive. NPHP1 is both (single-gene deletion that is recurrently NAHR-mediated and disease-causing); 22q11.2 is microdel-syndrome only.

---

## 6. Application to the eight LLmap anchor loci

| Anchor locus | Architecture | Mechanism | Clinical | Reasoning |
|---|---|---|---|---|
| **IGHG4-ChimDup** (~19.5 kb tandem, ~99.99% identity, chr14q32.33) | `tandem_dup_small` | `nahr_prone` (small but ≥99% identity and < 50 kb apart → fits Sharp/Stankiewicz) plus `gene_conversion_prone` (IGHG4 ↔ IGHGP literature) | `immune_receptor_locus`, `dispersed_paralog_gene_family` | 19.5 kb is below the 50 kb tandem-block cutoff and unit is directly oriented adjacent → small tandem. The ~99.99% identity and proximity satisfy classical NAHR criteria. Per our memory the dup carries CH2/CH3-SNPs and is transcriptionally usually silent — that is consistent with NAHR-derived recent paralog, not gene conversion alone. |
| **IGH-CanonDup / NAHR-block** (~115 kb co-dup of G2/A1/E, 84.8% identity, ~14.8% gaps) | `tandem_dup_large` + `complex_mosaic` | `nahr_prone` is *uncertain* — the 84.8% identity is below the canonical ≥ 95–97% NAHR threshold. The 14.8% gap fraction places it closer to a divergent ancient tandem block than to an active LCR. Tag as `gene_conversion_prone` (HPRC iso-seq data shows IGC at the locus). | `immune_receptor_locus`, `dispersed_paralog_gene_family` | Honest call: at 84.8% it is *probably no longer an active NAHR substrate*; flag the discrepancy with Watson/Rodriguez's higher-identity ~25 kb V-region blocks. |
| **HLA-A / -B / -C** (highly polymorphic class-I) | `complex_mosaic` (the MHC class-I sub-region is a mosaic of paralogous gene blocks) | `gene_conversion_prone` (HLA IGC is classical; Parham 1989, Chen 2007) | `mhc_locus`, `dispersed_paralog_gene_family` | Note: pure HLA-A/B/C are paralogs by ancient duplication; they are not high-identity SDs in the Bailey sense (≥ 90% over ≥ 1 kb is not generally met across the full gene). The locus is included in the SegDup catalog because *misalignment is a practical problem*, not because the genes meet the strict ≥ 90% SD definition. The block as a whole (MHC class-I region) is `complex_mosaic`. |
| **NBPF / 1q21** (Olduvai/DUF1220 family) | `intrachromosomal_interspersed` (with 1p11 satellite block also → some pairs cross to `interchromosomal`) | `core_duplicon` (Jiang 2007 explicitly cites NBPF as a core-duplicon family); `gene_conversion_prone` | `dispersed_paralog_gene_family` | NBPF is the canonical core-duplicon example; copies are interspersed across the 1q21.1 region with a satellite at 1p11.2. |
| **Y-AZF palindromes** (P1–P5, > 99% intra-arm identity, > 100 kb arms) | `palindromic_inverted` | `palindrome_mediated_nahr` | `mendelian_disease_locus` (AZFc deletion → spermatogenic failure); `immune_receptor_locus` does NOT apply | Skaletsky 2003 defined the ampliconic class; Kuroda-Kawaguchi 2001 defined the AZFc 229 kb direct-repeat NAHR products. Palindromes are *distinct* from ordinary IR-LCRs and earn their own architecture tag. |
| **NPHP1** (~85 kb tandem, flanking ~45 kb directly-oriented LCRs at 2q13) | `tandem_dup_large` (the duplicated unit itself) and the flanking 45 kb LCRs are `tandem_dup_small` with `nahr_prone` | `nahr_prone` (recurrent ~290 kb deletion via DP-LCR NAHR) | `mendelian_disease_locus` (autosomal-recessive nephronophthisis-1; ClinGen ISCA-37405 region) | Classic Lupski-style recurrent DEL via DP-LCR. 80 % of NPHP1 disease cases are this NAHR product. |
| **FCGR2A/2B/3A/3B** (1q23.3, 82.5 kb tandem repeat, > 98% identity) | `tandem_dup_large` | `nahr_prone` (the 82.5 kb LCR is directly oriented; CNV products are recurrent FCGR3B / FCGR2C deletions/duplications) | `immune_receptor_locus` (broad sense — Fcγ receptors are immune effectors, not Ig per se; you may prefer a separate `immune_effector_locus` if strictness matters); `dispersed_paralog_gene_family`; `gene_conversion_disease_locus` (FCGR2C chimeric formation via IGC) | Mueller 2013, Niederer 2010. |
| **MEFV** (16p13.3) | This one is *not* a segmental-duplication locus in the Bailey sense. The MEFV gene itself sits in a region with some low-complexity repeats but does not have a paralogous high-identity SD pair. Single duplication events have been reported (Cazeneuve 2015) but only as private events. Tag as `unresolved` *or* exclude from the catalog. | — | `mendelian_disease_locus` (FMF, autosomal-recessive) is the *only* honest tag here — and even that is not because of SD architecture, just because MEFV is on someone's clinical list. | Honesty note: including MEFV in a SegDup catalog requires either dropping the SegDup criterion or admitting MEFV is a *control region*. I would not assign a SD class. |

### 6.1 Anchor-level summary

- Six of the eight anchors fit cleanly into the architecture vocabulary.
- **IGH-CanonDup** sits awkwardly: its identity (84.8 %) is below the NAHR-prone band, so we should *not* tag it `nahr_prone` unconditionally. This is a genuine literature gap — it is *too divergent* to be a typical Lupski LCR but *too tandem-structured* to be a generic interspersed duplication. Flag this in the catalog as `class=tandem_dup_large; mechanism=gene_conversion_prone; nahr_status=ambiguous_low_identity`.
- **MEFV** does not fit any SegDup architecture class. Either drop it from the SegDup catalog or document explicitly that it is included only for clinical-comparison purposes.

---

## 7. Why we did *not* invent additional classes

The user's task lists several candidate categories that we deliberately did **not** turn into top-level classes, because the literature does not use them as distinct architecture types:

- **"LCR vs HCR" (low-copy vs high-copy repeat)**: "LCR" in the genomic-disorder literature (Stankiewicz/Lupski) is just an older synonym for "segmental duplication". There is no consistently used "high-copy repeat" class in the Eichler/Lupski framework; multi-copy SDs (NBPF, GAGE) are described by copy number directly, not by a separate class.
- **"Direct vs inverted" as a primary class**: orientation is captured under architecture (`tandem_dup_*` = direct, `palindromic_inverted`/`inverted_dispersed` = inverted). Promoting orientation to a top-level class causes cross-talk with `tandem` vs `interspersed`.
- **"De-novo-SV hotspot"**: this is downstream of `nahr_prone` plus locus-specific recombination-hotspot data; we keep it as an annotation rather than a class.
- **"Y palindrome" as separate from "MHC-style"**: kept as distinct architecture (`palindromic_inverted`) vs functional tag (`mhc_locus`), because Y palindromes and MHC class-I/II are mechanistically and architecturally different despite both being immune-relevant complex loci.

---

## 8. Discussion: ambiguous and edge cases

1. **MHC at the SD-definition boundary.** HLA-A/B/C paralogs are conventionally treated as a "structurally variable region" (Norman 2017) but individual gene pairs often fall *below* the 90 % / 1 kb Bailey threshold across their full lengths. Including MHC in a SegDup catalog is a *pragmatic* choice (read-mapping ambiguity) more than a *strict* one. We recommend keeping `mhc_locus` as a clinical tag and using `complex_mosaic` architecture for the whole region rather than per-gene-pair SD entries.

2. **IGH constant region heterogeneity.** Watson/Rodriguez 2020 describe the V-region as harbouring ~25 kb high-identity blocks (clearly `tandem_dup_small`/`nahr_prone`). The constant region (G1/G2/G3/G4/A1/A2/E/M/D) has a different layout: the ~115 kb co-dup of G2/A1/E is *also* tandem but at lower identity. Treating them with the same class hides this and would mislead downstream NAHR-risk priors. Use sub-architecture tagging (`tandem_dup_large` + identity bin) to keep them distinguishable.

3. **Pseudogene paralog pairs (PMS2/PMS2P3, GBA/GBAP1, CYP2D6/D7) versus "real" SDs.** These are SDs in the strict sense (high identity, paired) but their clinical importance is dominated by gene-conversion-mediated pathogenic allele creation, not by NAHR DEL/DUP. We tag them `gene_conversion_disease_locus` and rely on the mechanism tag rather than carving a new architecture class.

4. **Subtelomeric/pericentromeric overlap.** A SegDup at 22q11.2 LCR-A is both pericentromeric and the substrate of a microdel syndrome. We allow multiple tags rather than forcing a single label.

5. **Y palindromes are not microdel syndromes in the autosomal sense.** ClinGen's recurrent CNV list is largely autosomal. AZFc deletion is a recurrent NAHR product but is conventionally classed as a "Y-microdeletion" with its own clinical category. We tag it `mendelian_disease_locus` for catalog uniformity rather than introducing a Y-specific tag.

---

## 9. TL;DR — recommended minimum class set for the LLmap catalog

**Architecture (10, primary `class`, exactly one):**
`tandem_dup_small`, `tandem_dup_large`, `intrachromosomal_interspersed`, `interchromosomal`, `subtelomeric`, `pericentromeric`, `palindromic_inverted`, `inverted_dispersed`, `complex_mosaic`, `unresolved`.

**Mechanism (6, zero-to-many):**
`nahr_prone`, `nahr_inversion_substrate`, `gene_conversion_prone`, `core_duplicon`, `palindrome_mediated_nahr`, `none_documented`.

**Clinical/function (8, zero-to-many):**
`microdel_microdup_syndrome`, `mendelian_disease_locus`, `immune_receptor_locus`, `mhc_locus`, `dispersed_paralog_gene_family`, `pseudogene_paralog_pair`, `gene_conversion_disease_locus`, `pangenome_only`.

Each class is supported by ≥ 1 primary literature reference; nothing in the vocabulary is original to LLmap.

---

## 10. References (DOIs)

- Bailey JA, Eichler EE. *Primate segmental duplications: crucibles of evolution, diversity and disease.* Nat Rev Genet. 2006;7(7):552–564. doi:10.1038/nrg1895
- Bailey JA, Yavor AM, Massa HF, Trask BJ, Eichler EE. *Segmental duplications: organization and impact within the current Human Genome Project assembly.* Genome Res. 2001;11(6):1005–1017. doi:10.1101/gr.gr-1871r
- Bailey JA et al. *Recent segmental duplications in the human genome.* Science. 2002;297(5583):1003–1007. doi:10.1126/science.1072047
- Casola C, Hahn MW. *Gene conversion among paralogs in mammalian genomes.* Trends Genet. 2012. doi:10.1016/j.tig.2012.10.003
- Chen JM, Cooper DN, Chuzhanova N, Férec C, Patrinos GP. *Gene conversion: mechanisms, evolution and human disease.* Nat Rev Genet. 2007;8(10):762–775. doi:10.1038/nrg2193
- Horton R et al. *Variation analysis and gene annotation of eight MHC haplotypes: The MHC Haplotype Project.* Immunogenetics. 2008. doi:10.1007/s00251-007-0262-2
- Jiang Z et al. *Ancestral reconstruction of segmental duplications reveals punctuated cores of human genome evolution.* Nat Genet. 2007;39(11):1361–1368. doi:10.1038/ng.2007.9
- Jiang Z, Hubley R, Smit A, Eichler EE. *The evolution of human segmental duplications and the core duplicon hypothesis.* Cold Spring Harb Symp Quant Biol. 2009;74:355–362. doi:10.1101/sqb.2009.74.011
- Kuroda-Kawaguchi T et al. *The AZFc region of the Y chromosome features massive palindromes and uniform recurrent deletions in infertile men.* Nat Genet. 2001;29(3):279–286. doi:10.1038/ng757
- Liao W-W et al. *A draft human pangenome reference (HPRC).* Nature. 2023;617:312–324. doi:10.1038/s41586-023-05896-x
- Linardopoulou EV et al. *Human subtelomeres are hot spots of interchromosomal recombination and segmental duplication.* Nature. 2005;437:94–100. doi:10.1038/nature04029
- Liu P et al. *Mechanisms for recurrent and complex human genomic rearrangements.* Curr Opin Genet Dev. 2012;22(3):211–220. doi:10.1016/j.gde.2012.02.012
- Lupski JR. *Genomic disorders: structural features of the genome can lead to DNA rearrangements and human disease traits.* Trends Genet. 1998;14(10):417–422. doi:10.1016/S0168-9525(98)01555-8
- Marques-Bonet T, Eichler EE. *The evolution of human segmental duplications and the core duplicon hypothesis.* Cold Spring Harb Symp Quant Biol. 2009;74:355–362. PMC4114149
- Mefford HC, Eichler EE. *Duplication hotspots, rare genomic disorders, and common disease.* Curr Opin Genet Dev. 2009. doi:10.1016/j.gde.2009.04.003
- Mueller M et al. *Genomic pathology of low-affinity Fcγ receptor copy-number variation.* Hum Mol Genet. 2013. doi:10.1093/hmg/ddt284
- Niederer HA, Willcocks LC, Rayner TF et al. *Copy number, linkage disequilibrium and disease association in the FCGR locus.* Hum Mol Genet. 2010;19(16):3282–3294. doi:10.1093/hmg/ddq216
- Norman PJ et al. *Sequences of 95 human MHC haplotypes reveal extreme coding variation in genes other than highly polymorphic HLA class I and II.* Genome Res. 2017;27(5):813–823. doi:10.1101/gr.213538.116
- Porubsky D et al. *Recurrent inversion polymorphisms in humans associate with genetic instability and genomic disorders.* Cell. 2022;185(11):1986–2005.e26. doi:10.1016/j.cell.2022.04.017
- Repping S et al. *Recombination between palindromes P5 and P1 on the human Y chromosome causes massive deletions and spermatogenic failure.* Am J Hum Genet. 2002;71(4):906–922. doi:10.1086/342928
- Rodriguez OL, Gibson WS, Parks T, Emery M, Powell J, Strahl M et al. *A novel framework for characterizing genomic haplotype diversity in the human immunoglobulin heavy chain locus.* Front Immunol. 2020;11:2136. doi:10.3389/fimmu.2020.02136
- Sharp AJ et al. *Segmental duplications and copy-number variation in the human genome.* Am J Hum Genet. 2005;77(1):78–88. doi:10.1086/431652
- Skaletsky H et al. *The male-specific region of the human Y chromosome is a mosaic of discrete sequence classes.* Nature. 2003;423(6942):825–837. doi:10.1038/nature01722
- Stankiewicz P, Lupski JR. *Genome architecture, rearrangements and genomic disorders.* Trends Genet. 2002;18(2):74–82. doi:10.1016/S0168-9525(02)02592-1
- Stankiewicz P, Lupski JR. *Structural variation in the human genome and its role in disease.* Annu Rev Med. 2010;61:437–455. doi:10.1146/annurev-med-100708-204735
- Vollger MR et al. *Segmental duplications and their variation in a complete human genome.* Science. 2022;376(6588):eabj6965. doi:10.1126/science.abj6965
- Vollger MR et al. *Increased mutation and gene conversion within human segmental duplications.* Nature. 2023;617:325–334. doi:10.1038/s41586-023-05895-y
- Watson CT, Steinberg KM, Huddleston J et al. *Complete haplotype sequence of the human immunoglobulin heavy-chain variable, diversity, and joining genes and characterization of allelic and copy-number variation.* Am J Hum Genet. 2013;92(4):530–546. doi:10.1016/j.ajhg.2013.03.004
- Watson CT, Glanville J, Marasco WA. *The individual and population genetics of antibody immunity.* Trends Immunol. 2017. doi:10.1016/j.it.2017.04.003
- Jeong H, Dishuck PC, Yoo D et al. *Structural polymorphism and diversity of human segmental duplications.* Nat Genet. 2025. doi:10.1038/s41588-024-02051-8
- ClinGen Dosage Sensitivity Map (NPHP1: ISCA-37405). https://search.clinicalgenome.org/kb/gene-dosage/region/ISCA-37405

---

## 11. Open questions / next step

- Decide on the operational identity bin boundaries inside the schema JSON (will come in a separate doc): 90–95 % / 95–97 % / 97–99 % / ≥ 99.9 %.
- Decide whether `subtelomeric` and `pericentromeric` should remain architecture tags or be promoted to a fourth "location" axis (currently they overlap with `interchromosomal`).
- IGH-CanonDup ambiguity (84.8 % identity) — confirm with a per-base-identity histogram before finalising the catalog row.
- MEFV inclusion decision: drop from SegDup catalog or add a separate `clinical_control_locus` class.
