# Architectural archetypes — type-level, not coordinate-bound

This directory contains **pure descriptive entries** for repeat-architecture
archetypes. Unlike the rest of `specific_loci/`, these entries are NOT
coordinate-bound. They describe the **architecture** of a repeat family
that occurs at many places in the genome, and serve as a reference for the
classifier when it has detected the architecture but the exact coordinate
is irrelevant or the relevant copies are scattered.

Each file uses `"coordinates": {"grch38": {"chr": "architectural_archetype", "start": 0, "end": 0}}`
as a sentinel — the runtime should ignore these entries during coordinate
intersection lookups and instead consult them only when the universal
classifier produces a partial match and needs to name the architectural
family.

## Files

| File | Archetype |
|------|-----------|
| `alpha_satellite_higher_order_repeat.json` | 171 bp monomer → chromosome-specific HOR → tandem array |
| `alpha_satellite_suprachromosomal.json` | suprachromosomal families SF1-SF5 (monomer-type strata) |
| `beta_satellite_archetype.json` | 68 bp Sau3A pentamer arrays (acrocentric short arms) |
| `gamma_satellite_archetype.json` | ~220 bp pericentric GC-rich tandem |
| `hsat2_hsat3_archetype.json` | CATTC / ATTCC pentamer classical satellites (C-band) |
| `mer22_subterminal.json` | MER22 ~600 bp subterminal dispersed family |
| `tar1_telomeric_repeat.json` | TAR1 subterminal repeat family |
| `5s_rdna_cluster.json` | 5S rRNA cluster, chr1q42.13 (coordinate-bound) |

Total: 8 architecture files (one of which — `5s_rdna_cluster.json` — IS
coordinate-bound because the 5S array has a single defined location on chr1).

## Why a separate directory

The other `specific_loci/` subdirectories all contain coordinate-bound
instance loci consumed by the Layer-2 lookup at mapping time. These
architecture entries are consumed only by:

- the `consensus_match` feature extractor (Layer 1) when it needs to attach
  a human-readable family name to a feature-only match;
- the Layer-3 active agent, when it wants to describe what it has found
  to the user in `XR:Z:` tags.

They do NOT trigger Layer-2 coordinate-intersection mapping behaviour.

## Primary sources

- Alexandrov et al. 1988 — alpha-satellite suprachromosomal families.
  PMID:3209091
- Willard 1985 — chromosome-specific alpha-sat HORs. PMID:2581134
- Brown 1979 — Sau3A / classical satellite. PMID:225837
- Moyzis et al. 1988 — TTAGGG telomere. PMID:3413074
- Altemose et al. 2022 — complete centromere catalogue. doi:10.1126/science.abl4178
- Nurk et al. 2022 — T2T-CHM13. doi:10.1126/science.abj6987
- Riethman et al. 2005 — subtelomeres / TAR1. doi:10.1101/gr.3284105
- Stults et al. 2008 — rDNA copy number. doi:10.1101/gr.072348.107
