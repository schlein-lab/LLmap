# Difficult-to-map regions — miscellaneous hard cases

This directory catalogues miscellaneous loci that are difficult for short-
and long-read mappers but don't fit neatly into the `segmental_duplications/`,
`low_complexity/`, or the canonical loci handled by sister directories
(`immunoglobulin/`, `mhc/`, `kir/`, `olfactory/`, `centromeres/`,
`telomeres/`, `ribosomal/`, `mitochondrial/`, `y_chromosome/`).

The themes here are:

1. **NUMTs** — nuclear copies of mitochondrial DNA. These create the
   chrM ↔ nuclear mismapping artefact responsible for a large fraction
   of false mtDNA heteroplasmy calls in clinical sequencing.
2. **V(D)J-recombined immune loci other than IGH/IGK/IGL** — the T-cell
   receptor loci TRA/TRD/TRB/TRG.
3. **Tandem paralog clusters with extreme intra-cluster identity** —
   histone HIST1, MUC mucin VNTR, ZNF KRAB-zinc-finger.
4. **Globin clusters** — HBA on 16p13.3 (alpha-globin) and HBB on
   11p15.5 (beta-globin).
5. **Y palindromes P1-P8** — AZF microdeletion drivers.
6. **HLA class I high divergence** — peptide-binding-groove exon
   identity-drop problem.
7. **Pericentric inversion polymorphisms** — chr1, chr8, chr9, chr16
   pericentric inversions.
8. **Heterochromatin size polymorphisms (qh+ / qh-)** — 1qh, 9qh, 16qh,
   Yqh classical C-band polymorphisms with megabase-scale length variation.

## Files

| File | Region | Why it's hard |
|------|--------|---------------|
| `NUMT_chr5_major.json` | chr5q33.1 (~6 kb) | largest single NUMT |
| `NUMT_genomewide_pattern.json` | pattern marker | landscape-level chrM paralogs |
| `MUC_cluster_11p15.json` | 11p15.5 | mucin VNTR exons |
| `ZNF_cluster_19q13_4.json` | 19q13.31-q13.43 | >200 KRAB-ZNF paralogs |
| `TRA_TRD_locus_chr14.json` | 14q11.2 | V(D)J + paralog density |
| `TRB_locus_chr7.json` | 7q34 | V(D)J + paralog density |
| `TRG_locus_chr7.json` | 7p14 | V(D)J + paralog density |
| `HIST1_cluster_chr6.json` | 6p22 | identical H4 paralogs |
| `HBB_globin_chr11.json` | 11p15.5 | HBG1/HBG2 near-identity |
| `HBA_globin_chr16.json` | 16p13.3 | -alpha3.7 / -alpha4.2 SDs |
| `chrY_AZF_palindromes.json` | Yq11.2 | P1-P8 palindromes |
| `HLA_class_I_divergence.json` | 6p22.1 + 6p21.33 | high SNP density |
| `chr1_pericentric_inversion.json` | chr1 | pericentric inv polymorphism |
| `chr8_pericentric_inversion.json` | chr8 | inv(8) polymorphism |
| `chr9_pericentric_inversion.json` | chr9 | inv(9)(p11q13) common |
| `chr16_pericentric_inversion.json` | chr16 | inv(16) polymorphism |
| `chr1qh_polymorphism.json` | 1q12 | qh+/qh- HSat2/3 size polymorphism |
| `chr9qh_polymorphism.json` | 9q12 | qh+/qh- HSat3 size polymorphism |
| `chr16qh_polymorphism.json` | 16q11.2 | qh size polymorphism |
| `chrYqh_polymorphism.json` | Yq12 | DYZ1/DYZ2 10-40 Mb polymorphism |

Total: 20 difficult-region files.

## Overlap with other directories

- **HLA class I**: full MHC handled in `../mhc/`.
- **chrY palindromes**: chrY MSY architecture in `../y_chromosome/`.
- **NUMTs**: chrM itself in `../mitochondrial/`.
- **TR loci**: parallel to the IG loci in `../immunoglobulin/`.
- **Heterochromatin qh polymorphisms**: the static cytogenetic extent
  of each block lives in `../low_complexity/`; this directory captures
  the **population-level size variation** that complicates SV calling.

## Primary sources

- Hazkani-Covo 2010 — NUMT catalogue. doi:10.1371/journal.pgen.1000834
- Skaletsky et al. 2003 — Y MSY palindromes. doi:10.1038/nature01722
- Rhie et al. 2023 — T2T-Y (HG002). doi:10.1038/s41586-023-06457-y
- Nurk et al. 2022 — T2T-CHM13. doi:10.1126/science.abj6987
- Altemose et al. 2022 — Complete centromere catalog. doi:10.1126/science.abl4178
- Horvath et al. 2003 — pericentromeric LCR architecture. doi:10.1101/gr.1531003
- She et al. 2004 — pericentromeric SD shores. doi:10.1038/nature02564
- Sudmant et al. 2013 — SD copy-number variation. doi:10.1126/science.1234148
- Lefranc et al. — IMGT TR locus references. doi:10.1093/nar/gkh752
- Higgs et al. — globin cluster organisation. doi:10.1146/annurev.genom.7.080505.115700
- Marques-Bonet & Eichler — KRAB-ZNF clusters. doi:10.1101/gr.6395807
