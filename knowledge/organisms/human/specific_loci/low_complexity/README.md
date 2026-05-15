# Low-complexity / heterochromatin blocks

This directory catalogues large blocks of low-complexity / classical-satellite /
acrocentric heterochromatin and the pericentromeric sub-centromeric LCR
shorelines that flank every human centromere. Many are gaps in GRCh38 or
are present but unmappable with standard seed-chain-extend due to monotonous
tandem satellite content. Each file describes ONE named block on GRCh38 with
its satellite composition and locus-specific mapping hints for LLmap.

Reads landing in these regions almost always need
`report_multi_position: true`, very high `max_occ`, and
`allow_high_mismatch: true`, because the underlying sequence is internally
near-identical across megabases of array.

## Files

### Classical heterochromatin C-bands (qh blocks)

| File | Region | Composition |
|------|--------|-------------|
| `chr1q12_pericentric_HSat2.json` | 1q12 | HSat2 + HSat3 |
| `chr9q12_pericentric_HSat3.json` | 9q12 | HSat3 (largest autosomal C-band) |
| `chr16q11_2_pericentric.json` | 16q11.2 | mixed alpha + HSat2/3 |
| `chrYq12_heterochromatin.json` | Yq12 | DYZ1 + DYZ2 |

### Pericentromeric sub-centromeric LCR shorelines (Horvath/She 2003-2004)

| File | Region | Composition |
|------|--------|-------------|
| `chr1p11_pericentric.json` | 1p11 | alpha-sat D1Z5 + SD |
| `chr2p11_pericentric.json` | 2p11 | alpha-sat D2Z + SD |
| `chr2q11_pericentric.json` | 2q11 | alpha-sat D2Z + SD |
| `chr3p11_pericentric.json` | 3p11 | alpha-sat D3Z1 + SD |
| `chr3q11_pericentric.json` | 3q11 | alpha-sat D3Z1 + SD |
| `chr4p11_pericentric.json` | 4p11 | alpha-sat D4Z1 + SD |
| `chr4q11_pericentric.json` | 4q11 | alpha-sat D4Z1 + SD |
| `chr5p11_pericentric.json` | 5p11 | alpha-sat D5Z2 + SD |
| `chr5q11_pericentric.json` | 5q11 | alpha-sat D5Z2 + SD |
| `chr6p11_pericentric.json` | 6p11 | alpha-sat D6Z1 + SD |
| `chr6q11_pericentric.json` | 6q11 | alpha-sat D6Z1 + SD |
| `chr7p11_pericentric.json` | 7p11 | alpha-sat D7Z1/D7Z2 + SD |
| `chr7q11_pericentric.json` | 7q11 | alpha-sat D7Z1/D7Z2 + SD |
| `chr8p11_pericentric.json` | 8p11 | alpha-sat D8Z2 + SD |
| `chr8q11_pericentric.json` | 8q11 | alpha-sat D8Z2 + SD |
| `chr9p11_pericentric.json` | 9p11 | alpha-sat D9Z1 + SD |
| `chr10p11_pericentric.json` | 10p11 | alpha-sat D10Z1 + SD |
| `chr10q11_pericentric.json` | 10q11 | alpha-sat D10Z1 + SD |
| `chr11p11_pericentric.json` | 11p11 | alpha-sat D11Z1 + SD |
| `chr11q11_pericentric.json` | 11q11 | alpha-sat D11Z1 + SD |
| `chr12p11_pericentric.json` | 12p11 | alpha-sat D12Z3 + SD |
| `chr12q11_pericentric.json` | 12q11 | alpha-sat D12Z3 + SD |
| `chr16p11_pericentric.json` | 16p11 | alpha-sat D16Z2 + SD |
| `chr17p11_pericentric.json` | 17p11 | alpha-sat D17Z1/D17Z1-B + SD |
| `chr17q11_pericentric.json` | 17q11 | alpha-sat D17Z1 + SD |
| `chr18p11_pericentric.json` | 18p11 | alpha-sat D18Z1 + SD |
| `chr18q11_pericentric.json` | 18q11 | alpha-sat D18Z1 + SD |
| `chr19p11_pericentric.json` | 19p11 | alpha-sat D19Z1/D19Z2/D19Z3 + SD |
| `chr19q11_pericentric.json` | 19q11 | alpha-sat D19Z2 + SD |
| `chr20p11_pericentric.json` | 20p11 | alpha-sat D20Z1/D20Z2 + SD |
| `chr20q11_pericentric.json` | 20q11 | alpha-sat D20Z1 + SD |
| `chrXp11_pericentric.json` | Xp11 | alpha-sat DXZ1 + SD |
| `chrXq11_pericentric.json` | Xq11 | alpha-sat DXZ1 + SD |
| `chrYp11_pericentric.json` | Yp11 | alpha-sat DYZ3 + SD |

### Acrocentric p-arms (NOR-bearing)

| File | Region | Composition |
|------|--------|-------------|
| `acrocentric_13p.json` | 13p | alpha + beta-sat + 45S NOR |
| `acrocentric_14p.json` | 14p | alpha + beta-sat + 45S NOR |
| `acrocentric_15p.json` | 15p | alpha + beta-sat + 45S NOR |
| `acrocentric_21p.json` | 21p | alpha + beta-sat + 45S NOR |
| `acrocentric_22p.json` | 22p | alpha + beta-sat + 45S NOR |

### Subtelomeric mosaics

| File | Region | Composition |
|------|--------|-------------|
| `chr3p_subtelomeric.json` | 3p subtelomere | TAR1 + degenerate TTAGGG |
| `chr4p_subtelomeric.json` | 4p subtelomere | TAR1 + degenerate TTAGGG |
| `chr15q26_subtelomeric.json` | 15q26 subtelomere | shared subterminal mosaic |

Total: 46 low-complexity-block files.

## Primary sources

- Nurk et al. 2022 — T2T-CHM13. doi:10.1126/science.abj6987
- Altemose et al. 2022 — Complete centromere catalog. doi:10.1126/science.abl4178
- Horvath et al. 2003 — pericentromeric LCR architecture. doi:10.1101/gr.1531003
- She et al. 2004 — pericentromeric SD shores. doi:10.1038/nature02564
- Sudmant et al. 2013 — SD evolution. doi:10.1126/science.1234148
- Rhie et al. 2023 — T2T-Y (HG002). doi:10.1038/s41586-023-06457-y
- Brown 1979 — Sau3A/HSat characterisation. PMID:225837
- Riethman et al. 2005 — subtelomere mapping. doi:10.1101/gr.3284105

## Conventions

- All coordinates are GRCh38; CHM13 v2.0 boundaries are given where the
  GRCh38 placeholder is mostly Ns.
- `mapping_hints.allow_high_mismatch: true` and `lambda_scale: 3.0` are
  consistent across the satellite arrays — chain extension cannot recover
  a unique alignment here.
- Acrocentric p-arms and the sub-centromeric SD shorelines are nearly
  indistinguishable between chromosomes; reads should be flagged ambiguous
  unless a chromosome-specific HOR variant or pseudogene pins them.
