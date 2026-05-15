# Human specific_loci — Layer 2 region database

This directory contains **instance-level** locus annotations for the human
genome. Each JSON file describes ONE named, coordinate-bound locus on
GRCh38 (with CHM13 coordinates where they materially differ).

This is Layer 2 in the LLmap three-layer region-knowledge model
(see `../../../EXTENDING.md`):

- **Layer 1** — Universal taxonomy (`../regions/*.json`): describes
  region *types* (`paralog_family`, `centromere_like`, etc.) via
  feature signatures. Type-level, not instance-level.
- **Layer 2 — this directory:** named real-world loci with concrete
  GRCh38 coordinates, expected paralog structure, locus-specific
  mapping hints, sources. Consumed at mapping time when a read lands
  in a flagged locus.
- **Layer 3** — Active agent (`--llm=auto`): live LLM consultation at
  mapping checkpoints.

## File format

Each per-locus JSON file extends the universal `region_entry_schema`
with two locus-specific fields: `coordinates` and `size_bp`, and
typically a `expected_paralogs` block where relevant.

```json
{
  "name": "IGH_locus_chr14",
  "taxonomy_id": "paralog_family",
  "description": "Immunoglobulin heavy chain locus; ~4 IGHG paralogs ...",
  "coordinates": {
    "grch38": {"chr": "chr14", "start": 105500000, "end": 107300000},
    "chm13":  {"chr": "chr14", "start":  99670000, "end": 101470000}
  },
  "size_bp": 1800000,
  "expected_paralogs": [
    {"name": "IGHM", "approx_offset_in_locus": 0},
    ...
  ],
  "expected_copy_count": "1 per haplotype but heterogeneous ...",
  "mapping_hints": {
    "lambda_scale": 1.2,
    "anchor_weight_scale": 0.5,
    "max_occ": 50000,
    "report_multi_position": true,
    "require_psv_disambig": true,
    "VDJ_special_handling": true
  },
  "sources": ["doi:..."],
  "confidence": 0.95
}
```

`taxonomy_id` must reference one of the universal entries in
`../../../schema.json` (`paralog_family`, `centromere_like`,
`terminal_repeat`, `tandem_repeat`, `custom`, etc.).

`feature_signature` is **omitted** from these per-locus files —
classification of a read into one of these loci happens by coordinate
intersection (Layer 2 lookup), not by feature predicate (which is
Layer 1's job).

## Coordinate provenance

GRCh38 ranges are drawn from UCSC's `cytoBand` and `gap` tracks
(centromere gap intervals, p- and q-telomere gap intervals) and from
published literature for the named loci (IGH, MHC, KIR, NORs).
CHM13 coordinates, where given, come from Nurk et al. 2022
(doi:10.1126/science.abj6987) and the T2T-CHM13 v2.0 assembly.

Acrocentric chromosomes (13, 14, 15, 21, 22) carry only a q-arm
telomere as a useful mapping target — their p-arms are
satellite/rDNA repeat arrays without a defined p-telomere boundary
in GRCh38. p-telomere files are therefore omitted for these five
chromosomes.

## Directory layout

```
specific_loci/
├── README.md                         # this file
├── centromeres/        # 24 files (chr1..22, chrX, chrY)
├── telomeres/          # 43 files (22 chrs × p,q + chrX p,q + chrY p,q
│                       #            minus 5 acrocentric p-arms)
├── immunoglobulin/     # IGH, IGK, IGL — 3 files
├── mhc/                # MHC class I, II, III — 3 files
├── kir/                # KIR cluster — 1 file
├── olfactory/          # large OR clusters — 3 files
├── ribosomal/          # 45S rDNA NORs on the 5 acrocentrics — 5 files
├── y_chromosome/       # PAR1, PAR2, AZF — 3 files
└── mitochondrial/      # chrM — 1 file
```

Total: 86 locus files (excluding README).

## Mapping-hint conventions used here

- `centromere_like` loci — `max_occ: 200000`, `lambda_scale: 3.0`,
  `report_multi_position: true`, `allow_high_mismatch: true`.
  Tandem-satellite heavy; chains often cannot uniquely place a read.
- `terminal_repeat` loci (telomeres) — `max_occ: 100000`,
  `lambda_scale: 2.0`, `report_multi_position: true`,
  `expected_unique: false`. TTAGGG hexamer is identical at every
  chromosome end.
- `paralog_family` loci (IGH/IGK/IGL, MHC, KIR, OR, Y palindromes) —
  `report_multi_position: true`, `require_psv_disambig: true`,
  `lambda_scale: 1.0-1.4` depending on intra-cluster identity.
- `tandem_repeat` loci (rDNA NORs) — very high `max_occ` (the 45S
  unit is present in 300-400 copies per haploid genome).
- `custom: mitochondrial` (chrM) — `circular: true`,
  `expected_copy_count: 100-10000` (mtDNA per cell varies by tissue),
  `report_multi_position: false` (chrM is unique within its 16.6 kb
  contig); flag NUMTs separately.

## Extending

To add a new locus:

1. Pick a directory (or create one if a new category is justified).
2. Copy an existing file with similar taxonomy_id as template.
3. Fill in `name`, `description`, `coordinates`, `size_bp`,
   `expected_paralogs` (if relevant), and cite sources with DOIs.
4. Validate with:
   ```bash
   python3 -c "import json; json.load(open('your_file.json'))"
   ```
5. Update this README's count and the entry in `../MODULE.md`.

## Sources at a glance

- Nurk et al. 2022 — T2T-CHM13 v2.0 — doi:10.1126/science.abj6987
- Miga et al. 2020 — centromeres of chrX, chr8 — doi:10.1038/s41586-020-2547-7
- Altemose et al. 2022 — complete centromere catalog — doi:10.1126/science.abl4178
- Moyzis et al. 1988 — TTAGGG telomere identification — PMID:3413074
- Watson et al. 2013 — IGH locus organisation — doi:10.1093/nar/gks1313
- Horton et al. 2004 — MHC sequence — doi:10.1038/nrg1489
- The MHC Sequencing Consortium 1999 — doi:10.1038/44853
- Trowsdale et al. 2001 — KIR/LRC — doi:10.1038/35038512
- Pyatkin et al. / Niimura & Nei 2003 — OR family — doi:10.1101/gr.1326003
- Stults et al. 2008 — rDNA copy number variation — doi:10.1101/gr.072348.107
- Ross et al. 2005 — chrY sequence — doi:10.1038/nature03440
- Skaletsky et al. 2003 — chrY MSY — doi:10.1038/nature01722
- Andrews et al. 1999 — chrM revised Cambridge reference — doi:10.1038/13779
