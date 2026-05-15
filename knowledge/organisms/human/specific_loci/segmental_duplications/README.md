# Segmental duplications — recurrent CNV / NAHR hotspots

This directory catalogues human segmental-duplication (SD) clusters that act as
NAHR (non-allelic homologous recombination) substrates and produce the
canonical recurrent microdeletion / microduplication syndromes, plus inversion
polymorphisms mediated by inverted SD pairs. Each file describes ONE named SD
locus on GRCh38 with its architecture (direct vs inverted, homology, copy
number) and locus-specific mapping hints for LLmap.

These entries are Layer-2 instance loci (per the
[parent README](../README.md)). The runtime classifier triggers
`require_psv_disambig` and `report_multi_position` on any read that lands in
one of these intervals, and tightens its expectation of multi-position chains.

## Files

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `1q21_1_distal_NBPF.json` | 1q21.1 BP3-BP4 | direct SDs | NBPF/NOTCH2NL |
| `1q21_1_proximal_TAR.json` | 1q21.1 BP2-BP3 | direct SDs | RBM8A (TAR syndrome) |
| `1q21_2_BCL9_region.json` | 1q21.2 | complex SD | continuation of 1q21 |
| `2q13_recurrent_CNV.json` | 2q13 | direct SDs | NPHP1 region |
| `3p25_3_VHL.json` | 3p25.3 | direct SDs | VHL, ITPR1 |
| `5p15_33_PROK2_TERT.json` | 5p15.33 | complex SD | TERT, SDHA paralogs |
| `7q11_23_Williams_Beuren.json` | 7q11.23 | direct SDs | ELN, GTF2I |
| `7q22_q31.json` | 7q22-q31 | complex | CFTR region |
| `8p23_1_defensin_inverted.json` | 8p23.1 | inverted SDs | DEFB cluster |
| `10q11_22_recurrent.json` | 10q11.22 | direct SDs | GPRIN2, NPY4R |
| `15q11_q13_PWS_AS.json` | 15q11-q13 BP1-BP3 | direct SDs | SNRPN, UBE3A |
| `15q13_3_CHRNA7.json` | 15q13.3 BP4-BP5 | direct SDs | CHRNA7 |
| `16p11_2_distal.json` | 16p11.2 distal | direct SDs | SH2B1 |
| `16p11_2_proximal.json` | 16p11.2 proximal | direct SDs | TBX6, KIF22 |
| `16p12_1_recurrent.json` | 16p12.1 | direct SDs | (two-hit ID) |
| `16p13_11_recurrent.json` | 16p13.11 | direct SDs | MYH11, NDE1 |
| `17p11_2_Smith_Magenis.json` | 17p11.2 | direct SMS-REPs | RAI1 |
| `17q11_2_NF1.json` | 17q11.2 | direct NF1-REPs | NF1 |
| `17q12_RCAD.json` | 17q12 | direct SDs | HNF1B |
| `17q21_31_inversion.json` | 17q21.31 | inverted SDs | KANSL1, MAPT |
| `22q11_21_DiGeorge.json` | 22q11.2 LCR22 A-D | direct SDs | TBX1, CRKL |
| `22q13_Phelan_McDermid.json` | 22q13.33 | subtelomeric | SHANK3 |
| `Xp22_3_STS.json` | Xp22.31 | direct VCX SDs | STS |
| `Xq28_recurrent.json` | Xq28 int22h / F8 | inverted SDs | F8, MECP2 |

### chr19 (most SD-dense autosome, largest ZNF paralog family)

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `chr19_19p13_3_subtelomeric_SD_cluster.json` | 19p13.3 tip | subtelomeric SDs | WASH, FAM138, DUX4-like |
| `chr19_19p13_3_GAMT_SD.json` | 19p13.3 | complex SDs | GAMT, MUM1 |
| `chr19_19p13_3_LILR_LAIR_xref.json` | 19q13.42 (LRC) | tandem paralog | LILR/LAIR (LRC) |
| `chr19_19p13_2_KISS1R_SD.json` | 19p13.2 | complex SDs | KISS1R |
| `chr19_19p13_2_NOTCH3_SD.json` | 19p13.12 | complex SDs | NOTCH3 |
| `chr19_19p13_11_INSR_SD.json` | 19p13.2 | complex SDs | INSR |
| `chr19_19p13_11_SLC25A23_SD.json` | 19p13.11 | direct SDs | SLC25A23, KLF1 |
| `chr19_19p12_ZNF_block.json` | 19p12 | tandem KRAB-ZNF | proximal ZNF cluster |
| `chr19_centromere.json` | 19cen | alpha-sat + SDs | D19Z3 |
| `chr19_19q11_pericentric_inversion.json` | 19q11 | inverted SDs | pericentric inv(19) |
| `chr19_19q13_11_CYP4F_SD.json` | 19q13.11 | tandem paralog | CYP4F2/3/8/11/12/22 |
| `chr19_19q13_12_ZNF45_cluster.json` | 19q13.12 | tandem KRAB-ZNF | ZNF45 cluster |
| `chr19_19q13_2_EML2_SD.json` | 19q13.2 | complex SDs | EML2, BCKDHA |
| `chr19_19q13_2_PSG_xref.json` | 19q13.2 | tandem paralog | PSG1-PSG11 |
| `chr19_19q13_31_ZNF_largest_cluster.json` | 19q13.31 | tandem KRAB-ZNF | **largest ZNF cluster** (~40-50) |
| `chr19_19q13_32_APOE_CETP_region.json` | 19q13.32 | tandem paralog | APOE/APOC, CETP, LIPC |
| `chr19_19q13_32_KIR_xref.json` | 19q13.42 (KIR) | xref → /kir/ | KIR cluster |
| `chr19_19q13_33_ZNF_tail_cluster.json` | 19q13.33-q13.41 | tandem KRAB-ZNF | ~30 KRAB-ZNF tail paralogs |
| `chr19_19q13_33_TGFB1_SD.json` | 19q13.2 | complex SDs | TGFB1, BCL3 |
| `chr19_19q13_41_CEACAM_SIGLEC_cluster.json` | 19q13.41 | tandem paralog | CEACAM / PSG / SIGLEC |
| `chr19_19q13_41_FCGBP_SD.json` | 19q13.41 | tandem paralog | FCGBP |
| `chr19_19q13_42_NLRP_cluster.json` | 19q13.42 | tandem paralog | NLRP2/4/5/7/8/9/11/12/13 |
| `chr19_19q13_42_NKG2_LRC_xref.json` | 19q13.42 (LRC) | tandem paralog | LRC overview |
| `chr19_19q13_42_43_PRSS_family.json` | 19q13.42-43 | complex SDs | PRSS / KLK paralogs |
| `chr19_19q13_43_TRIM_KRTAP_xref.json` | 19q13.43 | tandem paralog | KRTAP4/5 cluster |
| `chr19_19q13_43_subtelomeric_SD.json` | 19q13.43 tip | subtelomeric SDs | ZIM2 / PEG3 |

### chr20

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `chr20_20p13_subtelomeric_SD.json` | 20p13 tip | subtelomeric SDs | DEFB125-132 cluster |
| `chr20_20p13_OXT_AVP_SD.json` | 20p13 | inverted paralog | OXT/AVP |
| `chr20_20p12_3_recurrent_CNV.json` | 20p12.3 | direct SDs | JAG1 (Alagille) |
| `chr20_20p12_BMP2_SD.json` | 20p12 | direct SDs | BMP2 |
| `chr20_20p11_23_PAX1_SD.json` | 20p11.22 | complex SDs | PAX1, FOXA2 |
| `chr20_20p11_21_CDC25B_SD.json` | 20p11.21 | direct SDs | CDC25B |
| `chr20_20p11_21_FOXA2_FERMT1.json` | 20p11.21 | complex SDs | FERMT1, BMP2K |
| `chr20_centromere_pericentric.json` | 20cen | alpha-sat + SDs | D20Z2 |
| `chr20_20q11_21_pseudocentromeric_SD.json` | 20q11.1-q11.21 | complex SDs | pericentric belt |
| `chr20_20q11_21_pericentric_inversion.json` | 20q11.21 | inverted SDs | inv(20)(p12q11) + microdel/dup |
| `chr20_20q11_22_BCL2L1_SD.json` | 20q11.21-22 | direct SDs | BCL2L1, TPX2, ID1 |
| `chr20_20q11_22_pericentric_SD_belt.json` | 20q11.22 | complex SDs | distal pericentric belt |
| `chr20_20q13_12_HNF4A_SD.json` | 20q13.12 | complex SDs | HNF4A (MODY1) |
| `chr20_20q13_13_NCOA3_SD.json` | 20q13.12 | direct SDs | NCOA3 (AIB1) |
| `chr20_20q13_2_PTPN1_SD.json` | 20q13.13 | direct SDs | PTPN1 |
| `chr20_20q13_2_ZNF217_SD.json` | 20q13.2 | direct SDs | ZNF217, CYP24A1, BCAS1 |
| `chr20_20q13_32_GHRH_SD.json` | 20q13.32 | complex SDs | GHRH |
| `chr20_20q13_33_KCNQ2_SD.json` | 20q13.33 | complex SDs | KCNQ2, CHRNA4 |
| `chr20_20q13_33_subtelomeric_SD.json` | 20q13.33 tip | subtelomeric SDs | OPRL1, COL9A3 region |

Total: 24 (original) + 26 (chr19) + 19 (chr20) = 69 SD-locus files.

## Primary sources

- Bailey et al. 2002 — Recent segmental duplications in the human genome.
  Science 297:1003. doi:10.1126/science.1072047
- Sharp et al. 2005 / 2006 — discovery papers for many recurrent CNVs.
- Eichler 2019 / Sudmant et al. — SD update.
  doi:10.1038/s41576-019-0180-9 ; doi:10.1038/nature19075
- Nurk et al. 2022 — T2T-CHM13. doi:10.1126/science.abj6987
- Individual disorder references are listed in each per-locus file.

## Conventions

- All coordinates are GRCh38. Where T2T-CHM13 differs materially, see the
  parent module's `MODULE.md`.
- `homology_percent` is the **inter-SD** identity between the flanking duplicon
  copies that mediate NAHR; not the intra-SD self-similarity.
- `mapping_hints.report_multi_position: true` and `require_psv_disambig: true`
  are set on every entry — these are the defining behaviours for SD regions.
- `confidence` reflects how well-established the locus is in the literature
  (0.95 = canonical recurrent CNV, well-mapped breakpoints; 0.6-0.7 = SD
  hotspot recognised but architecture less precisely defined in GRCh38).
