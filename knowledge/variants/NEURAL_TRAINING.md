# Using variant priors as features for the LLmap neural encoder

LLmap's V1.0 design includes a sequence-to-region neural encoder (the
"AI mandatory" pillar). The variant-prior table from Layer 4 is a
natural auxiliary input for that encoder.

## Why variant density is a useful learned signal

Population variant density is **not** redundant with sequence content.
Two windows with identical k-mer composition can have very different
variant landscapes:

- A centromeric satellite array carries high SV density (DEL/INS/INV
  churn) but low SNV density per bp.
- A housekeeping CDS carries low SV density and low SNV density.
- An MHC class-II exon carries moderate SV density but extreme SNV
  density (driven by balancing selection).
- A segmental duplication carries elevated SV density due to NAHR.
- An olfactory-receptor cluster carries high pseudogenisation noise
  visible as a peculiar SNV/SV co-occurrence pattern.

The model can learn these signatures and use them to disambiguate
regions whose raw sequence signatures alone are confusable (e.g.
distinguish an SD from an old dispersed-repeat copy).

## Suggested feature vector per position bucket

```
[ shannon_5mer,            # from Layer-1 feature extractor
  gc_content,
  kmer_multiplicity_p95,
  log1p(n_snv_alt_alleles), # from Layer 4
  snv_max_af,
  log1p(n_sv_overlapping),
  sv_max_af,
  sv_type_onehot_DEL,        # 6 channels for VCF SV vocab
  sv_type_onehot_INS,
  sv_type_onehot_DUP,
  sv_type_onehot_INV,
  sv_type_onehot_CNV,
  sv_type_onehot_BND ]
```

The `sv_type_onehot_*` channels carry the per-type AF from
`sv_summary`, zero if the type is absent from this bucket.

## Training data

The training corpus is bucket-level, not read-level. For each position
bucket in the reference, the input is the feature vector above plus a
fixed window of one-hot encoded sequence; the target is the region
label produced by the deterministic classifier (Layer 1) augmented with
specific_loci labels (Layer 2) where applicable.

This is **distant supervision** — labels come from the rule-based
classifier, the model learns to reproduce them. Once accuracy is high
enough, the model also serves as a soft-label propagator for buckets
where the rules are ambiguous.

## Why this matters for the user's intuition

The user's original observation:

> If we know a 5 kb deletion exists at 0.3 % allele frequency in gnomAD
> at chr14:105.6 M, a read with a corresponding 5 kb soft-clip should
> not be penalised.

The deterministic Layer 4 lookup handles this case directly. The
neural encoder generalises the same logic: even at a position not
present in the catalog, the model can learn "this region looks
SV-prone given its sequence and the surrounding variant density,
relax the SV penalty here too". That is the bridge from a static
prior to a learned prior — closing the gap on novel populations and
on non-human genomes where catalog coverage is sparse.

## Out of scope (for this stub)

- Concrete model architecture (CNN vs Transformer vs Mamba) — decided
  in `LLmap/training/` not here.
- Pre-training corpus selection — needs separate proposal.
- Encoder distillation back to a rule-set for the deterministic Layer
  1 — interesting but later.
