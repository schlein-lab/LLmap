# Provenance Mode — Genome-Test Results (2026-06-13)

End-to-end validation of the provenance/contamination mode (Agent 1 detectors +
Agent 2 core/resolver/CLI). Goal: every read → exactly one provenance flag,
lossless Σ-invariant, quantifiable spectrum — tested on real genomes.

## Validation summary

| Test | Data | Result |
|---|---|---|
| **Synth ground-truth** | 100 reads (70 host / 15 dup / 15 chim) | **host=70, dup=15, chim=15, Σ=100** — derive-path + Σ deterministically correct |
| **Real whole-genome** | HiFi vs GRCh38, 143.318 reads | host 92.8% / rdna 6.5% / chim 0.7% / pseudo ~0 — **Σ==N lossless** |
| **Matched HG002 unique** | chr1:20-25Mb, HiFi vs Illumina-300x (streamed) | HiFi: host 98.3/chim 1.7/rdna 0.005 · Illumina: host ~100/rdna 0.004 — **Σ OK both** |
| **Matched HG002 repeat** | chr1:145-148Mb (1q21 segdup) | HiFi rdna 2.08% / chim 0.54 · Illumina rdna 0.025 — **Σ OK both** |

## Robust finding (across all platforms + regions)
**The Σ-invariant holds everywhere** — the lossless framework is proven e2e on
real matched data: every read lands in exactly one Layer-1 class, no read
dropped, on both HiFi and Illumina, unique and repeat regions.

## Honest finding: the raw BAM-derive heuristics are platform/aligner-dependent
The matched test (same sample, same region, two platforms) isolates this:
- **`chim` (SA-tag) is long-read-biased.** HiFi 1.7% vs Illumina ~0 — but this is
  REAL SV-spanning biology (long reads cross breakpoints → SA tags), NOT artefact
  chimera. Same artefact-vs-biology doctrine as VDJ/numt-mthet → needs the
  Layer-1 (`chim` artefact: cross-chrom/PCR) vs Layer-3 (`bio:sv` real SV) split.
  [Agent 1 extending `chimera_provenance`.]
- **`rdna` (MAPQ<5) is aligner-MAPQ-convention-dependent.** HiFi 2.08% vs Illumina
  0.025% in the SAME 1q21 segdup — ~80× — but driven by aligner MAPQ scales, not a
  clean repeat-ambiguity measure. The true repeat/TE surface is genome-intrinsic
  → needs a Dfam/RepeatMasker/segdup catalog + WaveCollapse spread-mass, NOT raw
  MAPQ. [Block 2.]
- **`dup` (0x400) needs duplicate-marking.** Both GIAB BAMs are unmarked → dup=0 on
  both (a data property, not a platform difference). Needs a position-heuristic
  (same start+strand+CIGAR) for unmarked BAMs. [Agent 2, small.]

## Data-driven priorities (confirmed by the test)
1. **Block 2 — TE/repeat (the addressable surface):** ~6.5% whole-genome /
   2% in segdup are repeat-ambiguous → genome-intrinsic TE-family catalog +
   WaveCollapse spread-mass (honest positional uncertainty vs faked unique hit).
   This is where LLmap structurally beats minimap2/BWA.
2. **Heuristic refinement:** `bio:sv` Layer-3 split (Agent 1) + position-dup
   (Agent 2).
3. **Block 3 — catalog classes:** numt (chrM+NUMT-BED), para (PSV run), exo
   (contaminant panel), mei (TE catalog) — currently correctly inert (`None`).

## Reproduce
`samtools view <GIAB-HG002-BAM-url> <region> | llmap provenance-spectrum --out-prefix P`
GIAB HG002 HiFi: AshkenazimTrio/.../PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.GRCh38.duplomap.bam
GIAB HG002 Illumina: .../NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam
