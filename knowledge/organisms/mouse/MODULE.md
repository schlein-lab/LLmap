# Mus musculus knowledge module

- **organism**: *Mus musculus* (laboratory mouse, C57BL/6J reference background)
- **reference build**: GRCm39 (Genome Reference Consortium Mouse Build 39, GCA_000001635.9), released 2020-06-24 by the Genome Reference Consortium
- **last_updated**: 2026-05-15
- **maintainer**: LLmap project (knowledge layer)

## Reference notes

GRCm39 differs from older mouse builds in ways that matter for mapping:

- **vs mm10 / GRCm38** (2011): GRCm39 fixes hundreds of mis-assemblies, adds the
  mitochondrial genome to the main FASTA, closes many gaps in centromeric and
  pericentromeric regions, and corrects the chromosome 13/14 minor-satellite
  arrays. Coordinates are NOT one-to-one with mm10; a liftover is required
  for any externally-curated annotation track.
- **vs mm9 / NCBI Build 37** (2007): much larger differences in repeat content,
  Y chromosome, and IGH/IGK/IGL loci. Coordinates are completely incompatible.
  Any pre-2012 paper coordinates need liftover via UCSC chain files.

When pinning region coordinates, this module assumes GRCm39 chromosome
naming (`chr1`, ..., `chr19`, `chrX`, `chrY`, `chrM`). Source citations are
given as DOIs of foundational papers; coordinate-level data was curated
against UCSC mm39 / Ensembl 110 tracks.

## Region inventory

### Universal-taxonomy regions

- `regions/unique_single_copy.json` — sequences with a single best position in GRCm39
- `regions/low_complexity.json` — homopolymer / short-tandem stretches
- `regions/coding.json` — protein-coding intervals using mouse codon usage
- `regions/pseudogene.json` — disrupted coding-like sequence
- `regions/centromere_like.json` — pericentromeric satellite-block assemblies

### Tandem repeats specific to mouse

- `regions/tandem_repeat_minor_satellite.json` — 120 bp centromere-core repeat
- `regions/tandem_repeat_major_satellite.json` — 234 bp pericentromeric repeat
- `regions/tandem_repeat_telomeric.json` — TTAGGG telomere (vertebrate-conserved)

### Paralog families

- `regions/paralog_family_immunoglobulin.json` — IgH (chr12), IgK (chr6), IgL (chr16)
- `regions/paralog_family_mhc_h2.json` — H-2 complex on chr17
- `regions/paralog_family_olfactory.json` — olfactory receptor clusters (genome-wide)

### Dispersed (interspersed) repeats — mouse-specific

- `regions/dispersed_repeat_b1.json` — B1 SINE (mouse Alu-equivalent)
- `regions/dispersed_repeat_b2.json` — B2 SINE
- `regions/dispersed_repeat_line1.json` — L1Md (mouse-specific L1 lineage)
- `regions/dispersed_repeat_iap.json` — Intracisternal A-Particle ERV

## Files

- `MODULE.md` (this file)
- `regions/*.json` (15 region entries — see inventory above)
- `classifier_rules.json` — deterministic feature -> region rules
- `codon_table.json` — *Mus musculus* codon usage (GenBank standard table 1, mouse-weighted frequencies)

## Sources

Core foundational references cited across the region files:

- Waterston RH et al. (2002) *Initial sequencing and comparative analysis
  of the mouse genome*. **Nature** 420:520-562. doi:10.1038/nature01262
- Mouse Genome Sequencing Consortium / Church DM et al. (2009) *Lineage-specific
  biology revealed by a finished genome assembly of the mouse*. **PLoS Biology**
  7:e1000112. doi:10.1371/journal.pbio.1000112
- Frankish A et al. (2023) *GENCODE 2023*. **Nucleic Acids Research**
  51:D942-D949. doi:10.1093/nar/gkac1071 (GENCODE M33 = mouse GENCODE 33 for GRCm39)
- Pertea M, Salzberg SL (2010) *Between a chicken and a grape: estimating
  the number of human genes*. **Genome Biology** 11:206.
  doi:10.1186/gb-2010-11-5-206 (gene-count context applied to mouse here)
