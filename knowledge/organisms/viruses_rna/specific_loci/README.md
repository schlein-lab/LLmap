# RNA viruses — Layer 2 region database

Instance-level locus annotations for RNA viruses of medical importance.
Each JSON file describes ONE named region on a specific RefSeq genome
accession; coordinates are given relative to that accession (no
multi-assembly mapping — viral RefSeq is the reference).

This mirrors the Layer 2 model used for the human/great_apes/mouse
databases (see `../../human/specific_loci/README.md`) with two
viral-specific additions to every record:

- `virus_type` — Baltimore class string, e.g. `ssRNA(+)`, `ssRNA(-)`,
  `ssRNA(-)_ambisense`, `ssRNA_RT`.
- `segmented: true` on segmented genomes (Influenza A/B/C, Lassa).
  Multi-segment "whole virus" records carry a synthetic
  `multi_segment` coordinate stub.

Lineage / subtype / genotype entries (Pango lineages, HA/NA subtypes,
HCV genotypes, Lassa lineages I-VII, HIV-1 groups M/N/O/P) use
`taxonomy_id: paralog_family` to express the "many related variants
of one reference genome" relationship; lineage-defining mutations are
listed in `lineage_defining_mutations` where applicable.

## Directory layout

```
specific_loci/
├── README.md                 # this file
├── coronaviridae/            # 54 — SARS-CoV-2 (genome+regions+lineages), SARS-CoV-1, MERS
├── retroviridae/             # 31 — HIV-1 (HXB2 ref + regions + groups/subtypes), HIV-2
├── orthomyxoviridae/         # 46 — Flu A 8 segments + H1-H18 + N1-N11 + combos, Flu B, Flu C
├── flaviviridae/             # 38 — HCV (genome+regions+gts 1-7), Dengue 1-4, Zika, WNV
├── filoviridae/              # 16 — Ebola Zaire (genome+regions+editing site) + 5 species
└── arenaviridae/             # 15 — Lassa L+S segments + features + lineages I-VII
```

Total: **200** locus files (excluding README).

## Mapping-hint conventions

Three preset patterns dominate this database; pick by region biology:

- **Conserved structural / catalytic** — `allow_high_mismatch: false`,
  `require_psv_disambig: true`, `max_occ: 3000-8000`. Used for: RdRp,
  protease, IRES, frameshift elements, LTR core, packaging signals,
  nucleocapsid, M/E. Roughly 90 of 200 records.
- **Lineage / subtype / genotype** — `report_multi_position: true`,
  `require_psv_disambig: true`, `allow_high_mismatch: false`,
  `max_occ: 10000-20000`. Used for: SARS-CoV-2 Pango lineages,
  influenza HA/NA subtypes, HCV genotypes, HIV-1 groups, Lassa
  lineages. Roughly 102 records.
- **Hypervariable** — `allow_high_mismatch: true`,
  `report_multi_position: true`, `lambda_scale: 1.5-2.0`. Used for:
  HIV-1 env V1-V5, HCV HVR1/HVR2. Roughly 8 records — RNA viruses
  store most of their variability in lineage diversity rather than
  intra-locus hypervariability, so this bucket is small by design.

Tiny genomes (≤30 kb) keep `max_occ` modest (≤20k); the budget per
hit is set against the expected per-host viral copy number rather
than mapping-target multiplicity.

## Scope coverage

- **SARS-CoV-2**: full Wuhan reference + 21 functional regions
  (5'UTR/TRS-L, ORF1ab + -1 PRF + RdRp + helicase, full spike with
  NTD/RBD/RBM/furin/S2/HR1-HR2, ORF3a-10, 3'UTR/s2m) + 20 Pango
  lineages from B.1 to XEC + XBB recombination breakpoint atlas.
- **SARS-CoV-1, MERS**: reference + spike (RBD), frameshift, N gene.
- **HIV-1**: HXB2 reference + 21 regions including LTRs, TAR, Ψ,
  RRE, gag-pol PRF, protease/RT/integrase, env V1-V5 + CD4bs + MPER,
  tat/rev/nef + 7 group/subtype entries; HIV-2 reference.
- **Influenza A**: 8 segment references + H1-H18 + N1-N11 + 7
  pandemic/zoonotic reassortment strain records (1918, 1957, 1968,
  pdm09, H5N1 2.3.4.4b, H7N9, H3N2 seasonal). Influenza B (Victoria
  + Yamagata) and C references.
- **HCV**: H77 reference + 15 regions (IRES, all NS proteins, HVR1/2,
  3'UTR X-tail) + 7 genotype entries.
- **Flaviviruses**: Dengue 1-4 references + 4 UTR/sfRNA structural
  features (DENV-2 coords); Zika reference + envelope/NS1/3'UTR;
  West Nile reference + envelope/3'UTR-sfRNA.
- **Ebola**: Zaire reference + 10 regions including the GP poly-A
  editing site + 4 additional species (Sudan, Bundibugyo, Taï Forest,
  Reston) + Marburg.
- **Lassa**: L + S segment references, 6 functional sub-regions
  (ambisense intergenic hairpins, NP, GPC, L, Z), and 7 phylogeographic
  lineages.

## Sources

- NCBI Virus / RefSeq accessions cited per-file
- ICTV Master Species List
- Pango lineage designations + NextStrain Nextclade
- Los Alamos HIV Sequence Database
- WHO/OIE/FAO influenza nomenclature
- Smith 2013 HCV genotypes (doi:10.1002/hep.26744)
- Carroll 2015 Ebola Makona (doi:10.1038/nature14594)
- Lauring 2020 SARS-CoV-2 evolution
- doi:10.1126/science.aaf6116 (Zika)
- doi:10.1016/j.cell.2015.07.020 (Lassa lineages)
