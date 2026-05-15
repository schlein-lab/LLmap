# Layer 4 — Population variant priors

This directory adds a fourth knowledge layer to LLmap, sitting alongside the
three layers documented in `knowledge/EXTENDING.md`:

```
Layer 1: Universal taxonomy        (feature -> region type)
Layer 2: Per-locus database        (named loci, paralog structure, PSVs)
Layer 3: Live agent (--llm=auto)   (on-demand LLM consultation)
Layer 4: Population variant priors (THIS layer)
```

Layer 4 fits between Layer 2 and Layer 3 in priority: it is consulted at
mapping time after the deterministic taxonomy/locus rules but before the
live agent is invoked. Unlike Layer 3, Layer 4 is **fully deterministic**
once the variant catalogs have been ingested into a local prior table.

## Why population variant frequencies matter for mapping

When a read disagrees with the reference at a position, a classical mapper
treats every mismatch (or soft-clip, or split) as an alignment penalty. This
is correct in homozygous-reference loci, but wrong at polymorphic sites:

- A SNV at MAF 0.45 is roughly as likely to be alt as ref. Penalising it
  costs sensitivity for no information gain.
- A known 5 kb deletion at gnomAD AF 0.003 explains a soft-clip block in
  a long read without invoking a novel-SV hypothesis.
- A polymorphic STR (e.g. CAG repeat in HTT) has a population length
  distribution. Reads whose length falls inside that distribution should
  not be penalised against the single reference allele.

If the mapper knows the expected variant landscape at a locus, it can
discount expected disagreement and conserve its penalty budget for
disagreement that is genuinely novel (i.e. evidence of a real variant in
the sample, or a mis-mapping).

## What "consume at runtime" means here

Variant ingestion is a **one-time activity at module-build time**, very much
like Layer 1 region-knowledge curation. Pipelines look like:

1. Operator runs `llmap variant-ingest --source gnomad_sv_v4.1 --vcf <...>
   --output gnomad_sv.priors`.
2. The ingest pass scans the VCF/BCF/BED, bucketises per chromosome at
   100 bp resolution (configurable), aggregates AF and SV-type counts, and
   writes a flat indexed table (see `PRIOR_FORMAT.md`).
3. At mapping time, `llmap align --variant-priors gnomad_sv.priors,
   dbsnp.priors,...` memory-maps the table(s) and uses them as a read-only
   lookup. No live API call is ever made for Layer 4.

The runtime is fully offline and deterministic. The same input + the same
prior table gives the same output forever.

## Relationship to the other layers

| Layer | Determinism | Network? | Trigger |
|-------|-------------|----------|---------|
| 1 — universal taxonomy | deterministic | offline | every read |
| 2 — per-locus database | deterministic | offline | reads in flagged loci |
| 3 — live agent | non-deterministic | online | ambiguity / SD / unknown |
| 4 — variant priors | deterministic | offline | every aligned position (cheap lookup) |

Priority order when Layer-4 priors are loaded:

```
classical scoring
    + Layer 1 mapping_hints
    + Layer 2 locus overrides
    + Layer 4 per-position prior discount
    -> if still ambiguous and --llm != off -> Layer 3
```

Layer 4 acts as a **scoring modifier**, not a classifier — it changes
penalty weights, it doesn't change which region a position belongs to.

## Activation flags

```
llmap align ... --variant-priors PATH[,PATH,...]
             [--variant-snv-discount 0.7]   # fraction of normal mismatch
             [--variant-sv-discount  0.5]   # fraction of normal SV penalty
             [--variant-min-af 0.001]       # ignore ultra-rare variants
```

If no `--variant-priors` is passed, Layer 4 is entirely inactive and the
mapper behaves as before.

## Neural encoder training note

The per-position prior vector (n_snv_alt_alleles, snv_max_af,
n_sv_overlapping, sv_max_af, sv_summary one-hot) is itself a feature.
Sequence-to-region neural encoders trained for LLmap can take this as an
additional input channel; variant density patterns are themselves
informative about region biology (centromeric satellites carry far higher
SV density than housekeeping CDS, dispersed repeats have characteristic
SNV densities, segmental duplications have characteristic AF spectra).

See `NEURAL_TRAINING.md` for details.
