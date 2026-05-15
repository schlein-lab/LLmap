# Low-complexity / heterochromatin blocks

This directory catalogues large blocks of low-complexity / classical-satellite /
acrocentric heterochromatin that are either gaps in GRCh38 or are present but
unmappable with standard seed-chain-extend due to monotonous tandem
satellite content. Each file describes ONE named block on GRCh38 with its
satellite composition and locus-specific mapping hints for LLmap.

Reads landing in these regions almost always need
`report_multi_position: true`, very high `max_occ`, and
`allow_high_mismatch: true`, because the underlying sequence is internally
near-identical across megabases of array.

## Files

| File | Region | Composition |
|------|--------|-------------|
| `chr1q12_pericentric_HSat2.json` | 1q12 | HSat2 + HSat3 |
| `chr9q12_pericentric_HSat3.json` | 9q12 | HSat3 (largest autosomal C-band) |
| `chr16q11_2_pericentric.json` | 16q11.2 | mixed alpha + HSat2/3 |
| `chrYq12_heterochromatin.json` | Yq12 | DYZ1 + DYZ2 |
| `acrocentric_13p.json` | 13p | alpha + beta-sat + 45S NOR |
| `acrocentric_14p.json` | 14p | alpha + beta-sat + 45S NOR |
| `acrocentric_15p.json` | 15p | alpha + beta-sat + 45S NOR |
| `acrocentric_21p.json` | 21p | alpha + beta-sat + 45S NOR |
| `acrocentric_22p.json` | 22p | alpha + beta-sat + 45S NOR |
| `chr3p_subtelomeric.json` | 3p subtelomere | TAR1 + degenerate TTAGGG |
| `chr4p_subtelomeric.json` | 4p subtelomere | TAR1 + degenerate TTAGGG |
| `chr15q26_subtelomeric.json` | 15q26 subtelomere | shared subterminal mosaic |

Total: 12 low-complexity-block files.

## Primary sources

- Nurk et al. 2022 — T2T-CHM13 (resolves all these blocks). doi:10.1126/science.abj6987
- Altemose et al. 2022 — Complete centromere catalog. doi:10.1126/science.abl4178
- Brown 1979 — early Sau3A/HSat characterisation. PMID:225837
- Riethman et al. 2005 — subtelomere mapping. doi:10.1101/gr.3284105

## Conventions

- All coordinates are GRCh38. Many of these blocks are partly or fully
  represented as `gap` features in GRCh38; the start/end ranges given here
  cover the published cytogenetic extent of the block including assembled
  flanking sequence. The actual satellite array is shorter and is fully
  resolved only in T2T-CHM13.
- `mapping_hints.allow_high_mismatch: true` and `lambda_scale: 3.0` are
  consistent across the satellite arrays — chain extension simply cannot
  recover a unique alignment here.
- Acrocentric short arms (13p, 14p, 15p, 21p, 22p) are nearly indistinguishable
  by sequence between chromosomes; a read mapping to any acrocentric p-arm
  should be flagged ambiguous unless a chromosome-specific HOR or pseudogene
  pins it.
