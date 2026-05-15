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

Total: 24 SD-locus files.

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
