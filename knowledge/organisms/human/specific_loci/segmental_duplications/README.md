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

### chr16 (densest SD chromosome in the genome; LCR16a duplicon family)

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `16p11_2_ATXN2L_SD.json` | 16p11.2 distal | direct SDs | ATXN2L, SH2B1 |
| `16p11_2_BOLA2_BP4_BP5.json` | 16p11.2 BP4-BP5 duplicons | direct SDs | BOLA2A/B, SLX1, SULT1A3/A4 |
| `16p11_2_BP2_BP3_SD.json` | 16p11.2-p12.1 BP2-BP3 | direct SDs | LCR16a core |
| `16p11_2_NPIPB_SD.json` | 16p NPIPB array | direct SDs | NPIPB1-B15 |
| `16p11_2_SULT1A_paralog.json` | 16p11.2 SULT1A cluster | direct SDs | SULT1A1/A2/A3/A4 |
| `16p11_2_TBC1D10B_SD.json` | 16p11.2 BP4-BP5 internal | direct SDs | TBC1D10B, CORO1A |
| `16p11_2_TP53TG3_SD.json` | 16p11.2 TP53TG3 array | direct SDs | TP53TG3A-F |
| `16p12_1_MYH11_BP1_BP2.json` | 16p12.1 BP1-BP2 | direct SDs | MYH11 paralogs |
| `16p12_2_SD_blocks.json` | 16p12.2 (Bailey 2002) | direct SDs | OTOA, LCR16a |
| `16p13_1_LCR16a_NDE1.json` | 16p13.1 LCR16a | mosaic SDs | NDE1, MYH11 |
| `16p13_2_A2BP1_SD.json` | 16p13.2 RBFOX1 | intragenic SDs | RBFOX1 |
| `16p13_3_alpha_globin.json` | 16p13.3 HBA cluster | direct SDs (Z/X) | HBA1, HBA2 |
| `16p13_3_subtelomeric_SD.json` | 16p13.3 PKD1 paralogs | direct SDs | PKD1, TSC2 |
| `16p13_3_subtelomere_cap.json` | 16pter cap | subtelomeric SDs | - |
| `16p13_3_Rubinstein_Taybi.json` | 16p13.3 CREBBP | direct SDs | CREBBP |
| `16p11_2_distal_220kb_SH2B1.json` | 16p11.2 distal recurrent | direct SDs | SH2B1 (recurrence) |
| `16q11_2_pericentric_heterochromatin.json` | 16q11.2 pericentric | satellite + SD mosaic | (heterochromatin) |
| `16q12_SD_blocks.json` | 16q12 | mosaic SDs | CYLD, NKD1, SALL1 |
| `16q22_CTCF_SD.json` | 16q22 | direct SDs | CTCF |
| `16q23_3_WWOX_SD.json` | 16q23.3 FRA16D | intragenic SDs | WWOX |
| `16q24_ANKRD11_SD.json` | 16q24.3 | direct SDs | ANKRD11 |

### chr17 (SMS-REPs, NF1-REPs, CMT1A-REPs, KANSL1 H1/H2)

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `17p11_2_Potocki_Lupski.json` | 17p11.2 dup (reciprocal) | direct SMS-REPs | RAI1 |
| `17p11_2_distal_LCR.json` | 17p11.2 distal LCR | direct SDs | (atypical SMS BP) |
| `17p11_2_FLCN_SD.json` | 17p11.2 FLCN | direct SDs | FLCN |
| `17p11_2_LLGL1_SD.json` | 17p11.2 middle SMS-REP | direct SDs | LLGL1, FLII |
| `17p11_2_TRIM16_SD.json` | 17p11.2 TRIM16/16L | direct SDs | TRIM16, TRIM16L |
| `17p12_PMP22_CMT1A.json` | 17p12 CMT1A-REPs | direct SDs | PMP22 |
| `17p13_1_BCL6B_SD.json` | 17p13.1 | direct SDs | BCL6B, TP53 |
| `17p13_2_PRPF8_SD.json` | 17p13.2 | direct SDs | PRPF8 |
| `17p13_3_LIS1_SD.json` | 17p13.3 | direct SDs | PAFAH1B1 (LIS1), YWHAE |
| `17p13_3_Miller_Dieker.json` | 17p13.3 Miller-Dieker | direct SDs | PAFAH1B1, YWHAE |
| `17p13_subtelomeric_SD.json` | 17pter cap | subtelomeric SDs | - |
| `17q11_2_NF1_REPa.json` | 17q11.2 NF1-REPa | direct SDs | SUZ12P1 paralog |
| `17q11_2_NF1_REPc.json` | 17q11.2 NF1-REPc | direct SDs | SUZ12 paralog |
| `17q11_2_pericentric_SD.json` | 17q11.1-q11.2 pericentric | mosaic SDs | - |
| `17q21_31_KANSL1_SD.json` | 17q21.31 KANSL1 flank | inverted/direct (H1/H2) | KANSL1 |
| `17q21_31_Koolen_de_Vries.json` | 17q21.31 microdeletion | direct SDs (on H2) | KANSL1 |
| `17q21_32_LRRC37_SD.json` | 17q21.31-q21.32 LRRC37 | direct SDs | LRRC37A/A2/A3/B |
| `17q22_NOG_SD.json` | 17q22 | direct SDs | NOG, TRIM37 |
| `17q23_2_TBX2_TBX4_SD.json` | 17q23.2 SD flanks | direct SDs | TBX2, TBX4 |
| `17q23_TBX2_TBX4.json` | 17q23 microdeletion (rec.) | direct SDs | TBX2, TBX4 |
| `17q24_q25_NPLOC4_SD.json` | 17q24-q25 | mosaic SDs | NPLOC4 |
| `17q25_3_subtelomeric_SD.json` | 17qter cap | subtelomeric SDs | KRT cluster, SECTM1 |

### chr18

| File | Region | Mechanism | Key dosage genes |
|------|--------|-----------|-------------------|
| `18p11_32_subtelomeric_SD.json` | 18pter cap | subtelomeric SDs | - |
| `18p11_31_STK11_paralog.json` | 18p11.31 | direct SDs | STK11 paralog |
| `18p11_21_pericentric_SD.json` | 18p11.21 pericentric | mosaic SDs | - |
| `18q11_2_SD_module.json` | 18q11.2 | direct SDs | GATA6 |
| `18q12_DTNA_SD.json` | 18q12 | intragenic SDs | DTNA |
| `18q21_SMAD7_SD.json` | 18q21 | direct SDs | SMAD7 |
| `18q21_2_TCF4_SD.json` | 18q21.2 | intragenic SDs | TCF4 |
| `18q22_TNFRSF11A_SD.json` | 18q22 | direct SDs | TNFRSF11A (RANK) |
| `18q23_ZADH2_SD.json` | 18q23 subtelomeric | mosaic SDs | ZADH2, TSHZ1, CTDP1 |

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

Total: 24 (original) + 21 (chr16 new) + 22 (chr17 new) + 9 (chr18) + 26 (chr19) + 19 (chr20) = 121 SD-locus files (chr16-18 add ~52 new entries; chr16 alone is the densest SD chromosome in the human genome).

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
