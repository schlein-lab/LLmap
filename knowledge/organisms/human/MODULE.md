# Human (Homo sapiens) — LLmap region knowledge module

- **Organism:** Homo sapiens (taxid 9606)
- **Reference assemblies:**
  - **Primary:** GRCh38.p14 / hg38 (GCA_000001405.29). All region definitions in this module are tuned against GRCh38 coordinates and k-mer multiplicity statistics.
  - **Alternate:** T2T-CHM13 v2.0 (GCA_009914755.4). Where CHM13 materially differs from GRCh38, this module's region descriptions note the divergence — most notably for: centromeric satellite arrays (only fully resolved in CHM13), the IGH locus (CHM13 carries a duplicated IGHG4 paralog absent from the GRCh38 primary haplotype; see `regions/paralog_family_immunoglobulin.json`), the acrocentric short arms (β-satellite + rDNA + γ-satellite), and Yq12 heterochromatin (CHM13 resolves the AT-rich gap).
- **Curator:** Christian Schlein (christian.schlein@googlemail.com)
- **Last updated:** 2026-05-15
- **Schema version:** 0.1.0

## Coverage map

This module classifies windows of the human reference into the following region types. Total: 17 region entries grouped into 7 universal taxonomy categories.

| File | Universal taxonomy | What it covers | Approx. genome fraction |
|------|--------------------|----------------|-------------------------|
| `regions/unique_single_copy.json` | unique_single_copy | Most of the euchromatic genome; intergenic + intronic single-copy DNA | ~55% |
| `regions/coding.json` | coding | Protein-coding exons of ~20000 genes | ~1.5% |
| `regions/pseudogene.json` | pseudogene | ~14000 processed + duplicated pseudogenes | ~1% |
| `regions/low_complexity.json` | low_complexity | Homopolymers, simple repeats, AT-rich heterochromatin (e.g. Yq12) | ~3% |
| `regions/tandem_repeat_alpha_satellite.json` | tandem_repeat | 171 bp centromeric monomer, hierarchical HORs | ~3% |
| `regions/tandem_repeat_beta_satellite.json` | tandem_repeat | 68 bp acrocentric short-arm satellite | <1% |
| `regions/tandem_repeat_gamma_satellite.json` | tandem_repeat | 220 bp pericentromeric satellite | <1% |
| `regions/tandem_repeat_telomeric.json` | terminal_repeat | TTAGGG telomere hexamer | <0.1% |
| `regions/centromere_like.json` | centromere_like | Composite centromere/pericentromere region (broader; contains the satellite entries) | ~5-6% |
| `regions/paralog_family_immunoglobulin.json` | paralog_family | IGH (chr14), IGK (chr2), IGL (chr22) | <0.2% |
| `regions/paralog_family_mhc.json` | paralog_family | HLA class I + II on chr6p21 | <0.2% |
| `regions/paralog_family_kir.json` | paralog_family | KIR cluster on chr19q13.4 (LRC) | <0.01% |
| `regions/paralog_family_olfactory.json` | paralog_family | ~400 OR genes + ~600 OR pseudogenes across 50+ clusters | <1% |
| `regions/dispersed_repeat_alu.json` | dispersed_repeat | Alu SINE, ~1.1M copies | ~11% |
| `regions/dispersed_repeat_line1.json` | dispersed_repeat | L1 LINE, ~500k copies | ~17% |
| `regions/dispersed_repeat_ltr.json` | dispersed_repeat | HERVs + solo LTRs | ~8% |
| `regions/dispersed_repeat_sine.json` | dispersed_repeat | Non-Alu SINEs (MIR, MIR3) | ~2.5% |

Together repetitive + paralogous + low-complexity regions account for ~50% of the human genome (Lander 2001), which is why this module emphasises multi-position output and PSV-aware mapping for those classes.

## Other module files

- `classifier_rules.json` — priority-ordered decision rules mapping feature vectors to region names.
- `codon_table.json` — standard human codon usage frequencies; consumed by the universal `orf_density` extractor.
- `specific_loci/` — **Layer 2** per-locus database with concrete GRCh38 (and CHM13 where divergent) coordinates for 130 named loci: 24 centromeres, 43 telomeres (24 q-arms + 19 p-arms; acrocentric p-arms omitted), 3 immunoglobulin loci (IGH/IGK/IGL), 3 MHC sub-regions (class I/II/III), KIR cluster, 3 olfactory clusters, 5 acrocentric NORs (45S rDNA), Y-chromosome PAR1/PAR2/AZF, chrM, 24 segmental-duplication / recurrent-CNV regions, 12 low-complexity pericentric and acrocentric satellite blocks, and 8 additional difficult-region entries (NUMTs, palindromes, etc.). See `specific_loci/README.md` for the file format and source list.

## Notable omissions / candidates for future curation

The following are not curated as standalone entries in this v0.1 module and fall back to the universal taxonomy or to one of the existing entries:

- **Segmental duplications (SDs) > 1 kb / >90% identity** — currently absorbed by `paralog_family_*` where the SD overlaps a named locus and by `unique_single_copy` elsewhere. A dedicated `segmental_duplication.json` would be a natural Phase-2 addition (Bailey 2002, Eichler 2019).
- **rDNA arrays (45S/5S)** — large tandem repeats on the acrocentric short arms and chr1q42; deserve a dedicated entry once T2T-resolved coordinates are integrated.
- **Y-chromosome ampliconic palindromes (P1-P8)** — chrY-specific paralog family with very high intra-palindrome identity; needs its own entry.
- **Subtelomeric TAR1 / degenerate telomere variants** — currently fall under `low_complexity`.
- **Mitochondrial chromosome (chrM)** — separate genome; not handled by this nuclear-genome module.

## Primary sources

- Lander et al. 2001, doi:10.1038/35057062 — initial human genome sequencing and repeat composition.
- Nurk et al. 2022 (T2T consortium), doi:10.1126/science.abj6987 — complete sequence of a human genome.
- Bailey et al. 2002, doi:10.1126/science.1072047 — segmental duplications.
- Sudmant/Eichler et al. 2019 (and related), doi:10.1038/s41576-019-0180-9 — SD updates.
- Willard 1985, PMID:2581134 — alpha satellite chromosome-specific organisation.
- Moyzis et al. 1988, PMID:3413074 — TTAGGG telomere repeat identification.
- Batzer & Deininger 2002, doi:10.1038/nrg798 — Alu elements.
- Waterston (MGSC) 2002, doi:10.1038/nature01601 — for IGH genomic organisation context.
- The MHC sequencing consortium 1999, doi:10.1038/44853 — MHC genomic sequence.
- Trowsdale 2001, doi:10.1038/35038512 — KIR/LRC organisation.
- Niimura & Nei 2003, doi:10.1101/gr.1326003 — olfactory receptor family.
- Brown et al. 1979, PMID:225837 — early Sau3A/beta-satellite characterisation.
