# Genome-mode long-read mapping: locus selection by alignment quality

Status: **bug diagnosed, fix specified, NOT yet implemented.** main (`4e76de8`) is
clean — B2 transcript recall is 0.90, nothing broken committed. This note preserves
the diagnosis for a fresh focused session.

## Symptom

`llmap align --mode reads` (genome) maps **0%** of real HiFi reads >~2 kb (identity
0.000), while minimap2 maps the same reads 100%. transcript-mode (B2) is unaffected.

## Diagnosis chain (exemplary — each step refuted a hypothesis)

1. **Length-dependent:** a read works at 2 kb, fails at ≥5 kb (identity → 0.000).
2. **Refuted `max_alignment_score`** (raised to 1e6 → still 0.000).
3. **Refuted seed-explosion / `max_occ`** (lowered 5000→200: hits 128k→1892 but the
   chain + result were IDENTICAL → not a frequency problem; the spurious anchors
   survive any cutoff → a LOCAL low-complexity repeat, global occ < 200).
4. **Refuted soft-masking / case** (uppercased ref+read → still 0.000).
5. **Diagonal dump:** the chain diagonal (q−r) is ~constant (drift ~7 over 5 kb) →
   the chain *looks* colinear. But the alignment accounting gave M=3001/I=1999 for a
   q-span≈r-span≈4955 → the per-anchor `min(q,r)` interpolation FABRICATED matches.
6. **Honest WFA (per-window) revealed the truth: M=1** per 1000 bp window → the
   anchors' surrounding sequence does NOT align. The earlier "45 kb → 1.000" was a
   DEGENERATE 1-match alignment (1/1), not a real fix. (The aligner is fine — B2
   uses it at 0.90.)
7. **mm2 truth comparison (the decisive step):** read `m64043_200711_235708/36242789`
   (29352 bp) maps **reverse (FLAG 16) @ chr1:19973045** (span 19973045..20002397).
   llmap's chain is **forward @ 19997333** — *inside* mm2's span (NOT a 24 kb-distant
   wrong locus; the right footprint, wrong strand). A ±90 bp search of read[4:4+k]
   against the ref around 19997333 finds **no match in either orientation** → the
   forward anchors are spurious (chance k-mer hits in a local repeat), and the read's
   true REVERSE chain is out-ranked.

## Root cause

genome-mode ranks candidate chains by **anchor count** and extends the top one. A
read crossing a local repeat gets a **dense spurious chain** (many chance k-mer hits,
here a forward chain where the read truly maps reverse) that out-scores the read's
**true (sparse-but-real) chain**. The spurious chain extends to a **degenerate M=1**
alignment → rejected on identity → and the true chain is either out-ranked beyond
`max_chains_to_extend` or never tried → 0% mapped. genome-mode has **no
locus-quality selection** (`SelectTranscriptLocus` is transcript-only). This is why
B2/transcript works (per-exon, strand-aware extension) but genome long-read doesn't.

`max_occ` / frequency cutoffs do NOT help — the repeat is local (global occ < 200).

## Fix (specified, minimap2-style)

Rank by **post-extension alignment quality, not anchor count**:
- keep top-N candidate chains (not just #1),
- extend each, compute identity,
- choose the best; a chain extending to ~0% identity (the degenerate M=1) is a
  **reject signal** → escalate to the next candidate.
- unify with the existing `sufficiency` work: extend top chain → sufficiency check
  (identity good enough?) → if not, escalate. The "solve-if-sufficient" loop, one
  stage earlier (locus choice).

This is **foundational for genomic assemble-then-map**: the `ProbeMapper` injected
into `MapClusteredReads` must do this quality locus-selection, else the probe also
lands on the repeat locus.

## What landed this session (for context)

- **B2 splice-aware chaining (Prio-1): transcript spliced-recall 0.20 → 0.90** —
  committed, the session's headline win.
- Provenance/M(pos)/pangenome baseline, Block-3 catalogs, the assemble-then-map
  module chain (cluster + MinHash sketch, lossless layout, propagation, WaveCollapse
  guard), sufficiency/consensus_resolution, Idee-2 planted-truth QC — committed/ready.
