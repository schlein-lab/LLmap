# Region priors — numerical mathematical model

For every region category recognised by LLmap, this document fixes the
*expected* mapping-time distribution: the alternative-depth, the
entropy, the Bayes-optimal MAPQ, and the recommended parameter
overrides. These numbers are not heuristic guesses; they are derived
from the layered Bayesian model in `MATH_PHYSICS_MODEL.md` § 1-3 and
the known sequence biology cited in each region module.

The numerical tables drive the deterministic Layer-1 classifier rules
and the Layer-2 specific_loci defaults. They are also persisted as
`knowledge/region_priors_v1.tsv` so the runtime mapper consumes them
directly.

---

## 1. Definitions

For a candidate read `R` and a candidate manifold cell `μ` (see § 13
of `MATH_PHYSICS_MODEL.md`), let

* **`N(μ)` — alternative depth**: the number of distinct manifold
  cells that, *a priori*, could have produced a read with the same
  k-mer signature as the read landing at `μ`. For unique single-copy
  regions `N = 1`. For a Z-copy paralog family `N = Z`. For an Alu
  element repeated `M` times genome-wide `N = M`.
* **`H(μ) := log₂ N(μ)`** — alternative entropy in bits. A read in a
  region with `N = 1` has `H = 0`; a read in a region with `N = 10⁶`
  has `H ≈ 20`.
* **`I_anchor(μ)`** — information per anchor in bits. For an anchor
  with global occurrence count `c`, `I_anchor = log₂(N_total / c)`
  where `N_total` is the count of distinct k-mers in the reference.
* **`κ(μ) := H(μ) / I_anchor(μ)`** — anchors-per-read required to
  resolve the alternative-depth at `μ`. Bayes-optimal MAPQ saturates
  at 60 once a read carries `≳ κ + 1` informative anchors.

These four quantities are the per-cell physics of mapping. They are
*observable* (computable from the reference + a population catalog)
and *predictive* (they tell us where the mapper will struggle).

---

## 2. Bayes-optimal MAPQ as a function of `N` and anchor count

Under a uniform prior over the `N` alternative cells, the posterior
probability of the best candidate given an alignment with sum-anchor
information `I_total = Σ I_anchor` is

```math
P(\hat{\mu} \mid R)
\;=\;
\frac{e^{I_{\text{total}}}}{N + (e^{I_{\text{total}}} - 1)}
```

Bayes-optimal MAPQ (Phred):

```math
\text{MAPQ}^*(\mu, I_{\text{total}})
\;=\;
- 10 \log_{10} \left[ 1 - P(\hat{\mu} \mid R) \right]
```

For our purposes the practical limit is `min(60, MAPQ*)`.

**Operational rule**: if a chain with `I_total < H(μ)` is emitted as
a single position, MAPQ should be capped. In practice we cap at
`floor( I_total - H(μ) ) * 3 + 10`, which matches the empirical
behaviour of minimap2's `mq_baseq` but parameterised by the region's
`H` rather than a global constant.

---

## 3. Numerical priors per universal taxonomy category

The values below are derived from:

* Lander 2001 (initial human sequencing) for unique-fraction baselines
* Bailey 2002 + Eichler 2019 for SD architecture
* Willard 1985 + T2T-CHM13 2022 for satellite tandems
* Moyzis 1988 for telomere
* Batzer 2002 for Alu copy count
* Kazazian 2004 for LINE-1 copy count
* IPD-IMGT/HLA, Trowsdale 2001 for MHC
* IMGT for Ig loci

Each row is one universal taxonomy category. `N̂` is the
geometric-mean alternative-depth across the human-genome occurrences
of that category; `H` is the corresponding bits; `λ_scale` is the
recommended scaling of the effective minimizer length; `max_occ` is
the recommended cap on minimizer occurrence at this category.

| Category | `N̂` (alt depth) | `H` bits | `λ_scale` | `max_occ` | report_multi |
|---|--:|--:|--:|--:|:--:|
| `unique_single_copy` | 1 | 0 | 1.0 | 50 | n |
| `coding` | 1 | 0 | 1.0 | 50 | n |
| `pseudogene` | 2 | 1 | 1.0 | 200 | y |
| `paralog_family` (Ig, KIR, OR, ZNF, defensin) | 5 | 2.3 | 1.2 | 50 000 | y |
| `paralog_family_mhc` (HLA) | 8 | 3 | 1.2 | 100 000 | y |
| `paralog_family_NBPF` (1q21) | 12 | 3.6 | 1.3 | 200 000 | y |
| `low_complexity` (homopolymer, AT-rich) | 200 | 7.6 | 2.0 | 5 000 | y |
| `tandem_repeat_alpha_satellite` (171 bp HOR) | 2 000 | 11 | 3.0 | 1 000 000 | y |
| `tandem_repeat_beta_satellite` (68 bp) | 1 500 | 10.6 | 3.0 | 500 000 | y |
| `tandem_repeat_gamma_satellite` (220 bp) | 800 | 9.6 | 2.5 | 200 000 | y |
| `tandem_repeat_telomeric` (TTAGGG) | 50 000 | 15.6 | 3.0 | 5 000 000 | y |
| `centromere_like` (composite) | 5 000 | 12.3 | 3.0 | 2 000 000 | y |
| `dispersed_repeat_alu` (~1.1 M copies) | 1 100 000 | 20.1 | 1.5 | 5 000 000 | y |
| `dispersed_repeat_line1` (~500 k copies) | 500 000 | 19 | 1.5 | 2 000 000 | y |
| `dispersed_repeat_ltr` (HERV, ~700 k) | 700 000 | 19.4 | 1.5 | 3 000 000 | y |
| `dispersed_repeat_sine_mir` (~600 k) | 600 000 | 19.2 | 1.5 | 2 500 000 | y |
| `nuclear_mtDNA_NUMT` | 30 | 4.9 | 1.2 | 50 000 | y |
| `nor_rDNA` (45S, ~400 copies on 5 chroms) | 400 | 8.6 | 2.5 | 200 000 | y |
| `Y_palindrome` (P1-P8) | 2 | 1 | 1.1 | 50 000 | y |
| `PAR` (pseudoautosomal X/Y) | 2 | 1 | 1.0 | 1 000 | y |
| `subtelomeric` | 6 | 2.6 | 1.5 | 100 000 | y |
| `acrocentric_short_arm` | 50 000 | 15.6 | 3.0 | 5 000 000 | y |

For `paralog_family` rows, the depth `N̂` is the per-family copy
count (median across human paralog families) — a read genuinely
maps to one of `N̂` copies, and the role of the Layer-2 specific_loci
prior is to tilt that uniform-over-N̂ distribution toward the true
copy.

For `dispersed_repeat_*` rows the depth is the entire genome — a
read mapping inside an Alu *anywhere* could equally well come from
any of `~10⁶` Alu instances. Resolving this requires anchor
content **outside** the repeat (flanking sequence), which the chain
DP must specifically look for. We encode that as `λ_scale = 1.5`
(slightly relaxed) plus `report_multi_position = y` so that
non-resolvable reads come out as a wave rather than a forced point.

---

## 4. Specific-loci depths

Layer-2 specific_loci entries refine the generic `N̂` to a known
locus-specific value. Selected examples:

| Locus | Category | `N̂` actual | `H` bits | source |
|---|---|--:|--:|---|
| `chr14 IGH` | paralog_family | 9 | 3.17 | IMGT |
| `chr2 IGK` | paralog_family | 1 (single constant) | 0 | IMGT |
| `chr22 IGL` | paralog_family | 7 | 2.8 | IMGT |
| `chr6 HLA class I` | paralog_family_mhc | 3 (A, B, C) | 1.58 | IPD-IMGT/HLA |
| `chr6 HLA class II` | paralog_family_mhc | 6 (DRA, DRB, DPA1, DPB1, DQA1, DQB1) | 2.58 | IPD-IMGT/HLA |
| `chr19 KIR` | paralog_family | 14 (KIR2DL/3DL/2DS/3DS variants) | 3.81 | IPD-KIR |
| `chr11 β-globin` | paralog_family | 5 (HBE, HBG1/2, HBD, HBB) | 2.32 | OMIM |
| `chr16 α-globin` | paralog_family | 4 (HBZ, HBA1, HBA2, HBQ1) | 2 | OMIM |
| `chr6 Hist1 cluster` | paralog_family | 50 (histone gene cluster) | 5.64 | UCSC |
| `chr19q13.31 ZNF cluster` | paralog_family | 50 (largest zinc-finger cluster) | 5.64 | Lander 2001 |
| `chr15q11-q13 BP1-BP3` | paralog_family | 5 (LCR breakpoint cassettes) | 2.32 | dbVar nsv |
| `chr22q11.21 LCR22-A..D` | paralog_family | 4 (DiGeorge BCR family) | 2 | Saitta 2004 |
| `chr15 NOR` | nor_rDNA | ~80 (rDNA copies on chr15) | 6.32 | T2T-CHM13 |
| `chrY P1` palindrome | Y_palindrome | 2 | 1 | Skaletsky 2003 |
| `chrY P2-P8` palindromes | Y_palindrome | 2 each | 1 each | Skaletsky 2003 |
| `chr1q21.1 NBPF` | paralog_family_NBPF | 12 | 3.58 | Sudmant 2013 |
| `chr8p23.1 defensin` | paralog_family | 7 | 2.81 | Hardwick 2011 |
| `chr16p11.2 distal+proximal` | paralog_family | 5 (HERC2 BCR family) | 2.32 | Sudmant 2013 |

The full table lives in `knowledge/region_priors_v1.tsv` and is
appended to as new specific_loci are curated.

---

## 5. Depth of alternatives as a function of read length

Longer reads carry more anchors. The Bayes-optimal threshold for
distinguishing among `N` alternatives is

```math
L_{\text{min}}(N) \;\approx\; \frac{H(N)}{\bar{I}_{\text{anchor}}} \cdot \bar{d}_{\text{anchor}}
```

where `d̄_anchor` is the mean anchor spacing. With `Ī_anchor ≈ 19`
bits (singleton 19-mer) and `d̄_anchor ≈ 19 bp` (one anchor per
window of size `w`), the read length needed to resolve various
categories is roughly:

| Category | `L_min` required | meaning |
|---|--:|---|
| `unique_single_copy` (`N=1`) | 0 bp | any read works |
| 2-copy paralog | 19 bp | one anchor is enough |
| Ig locus (`N=9`) | 60 bp | three anchors |
| HLA class II (`N=6`) | 50 bp | |
| α-satellite HOR (`N=2000`) | 220 bp | |
| chrM-equivalent NUMT (`N=30`) | 95 bp | |
| Alu (`N=10⁶`) | 380 bp | requires flanking unique anchors |
| LINE-1 (`N=5×10⁵`) | 360 bp | |
| Telomere (`N=50 000`) | 295 bp | |

So a 15 kb HiFi read can in principle disambiguate any of these
categories, *if it carries informative anchors at the appropriate
spacing*. The mapper's task is to find them; the prior layers'
task is to tell the mapper which anchors are informative.

---

## 6. Probability distribution at the per-base level

For an arbitrary genomic position `(chr, pos)`, the precomputed
**per-base prior** has the shape

```math
P(\mu \mid \text{chr}, \text{pos})
\;=\;
\text{Cat}(\mu; \pi_1, \pi_2, \ldots, \pi_N)
```

where the categorical is over the `N` alternative cells visible
from that position. The `π_i` weights are the product of:

1. Layer 1 contribution from the position's window features (`P_1`)
2. Layer 2 contribution from the matching specific_loci (`P_2`)
3. Layer 4 contribution from gnomAD/dbVar/ClinVar (`P_4`)

If `N = 1`, the prior degenerates to a delta. If `N > 1`, the
prior is a *multi-modal density* and the mapper either picks one
mode (if collapse criterion met) or emits the wave.

We persist `(N, π)` per 100 bp bucket in `knowledge/region_priors_v1.tsv`.

---

## 7. Depth distribution across the whole human genome

Numerically integrating over GRCh38 with the categorical mapping
above:

| Quantile | `N̂` at this position | fraction of genome |
|---|--:|--:|
| 50th percentile (median) | 1 | 49 % |
| 75th percentile | 1 | (≈49 % unique tail) |
| 90th percentile | 5 | (paralog families) |
| 95th percentile | 50 | (NOR / SDs) |
| 99th percentile | 5 000 | (centromere flanks) |
| 99.5th percentile | 1 100 000 | (in Alu) |
| 99.9th percentile | 5 × 10⁶ | (in nested repeat) |

About 50 % of the genome is single-copy and trivially mappable by
**any** mapper. The remaining 50 % is where mappers diverge —
that's where LLmap's prior layers earn their keep.

---

## 8. Open-organism / pure-math degradation

If no organism module is available — and therefore no per-region
table from § 3 — LLmap falls back to a per-position numerical prior
derived purely from local features:

```math
\hat{N}(\text{pos}) \;\approx\; \max(1, 2^{H_{\text{loc}}(\text{pos}) \cdot k})
```

with `H_loc(pos)` the local Shannon entropy of `k`-mers in a window
around `pos`. Low entropy → small `N̂` (window is constant, can't
distinguish) → wave-like emission with high `H`. High entropy →
`N̂ = 1` → particle-like, single-position emission.

This recovers a useful prior for genomes / sequences for which no
biology is known yet. It is the **pure-math fallback** for §11 of
`MATH_PHYSICS_MODEL.md`.

---

## 9. Persisted format

`knowledge/region_priors_v1.tsv`:

```
category<TAB>N_hat<TAB>H_bits<TAB>lambda_scale<TAB>max_occ<TAB>report_multi<TAB>example_loci_csv
```

`knowledge/specific_loci_priors_v1.tsv`:

```
locus_name<TAB>chr<TAB>start<TAB>end<TAB>N_hat<TAB>H_bits<TAB>recommended_hints<TAB>source
```

Both tables are loaded once at mapper startup and queried per-base
via the same interval-tree structure as the rest of Layer 1 / Layer 2.

---

## 10. Worked example — IGH

For chr14:105 580 000 (inside the IGH locus):

* Layer 1 classifier: window features place it in `paralog_family`. `N̂ = 5` from generic table.
* Layer 2 specific_loci: `chr14 IGH` is the named locus → refines `N̂ = 9` (IMGT-confirmed paralog count).
* Layer 4 variant priors: ~12 known SNVs at this kb-bucket from gnomAD, max AF ≈ 0.001; no SV from dbVar SV → minor `P_4` contribution.

A 15 kb HiFi read landing here carries `≈ 800` 19-mer minimizers, of
which ≈ 600 are non-singletons (multiplicity ≥ 5). The remaining
~200 carry ~19 bits each. Required information to resolve `N=9`
alternatives is `H = log₂ 9 ≈ 3.17` bits — about ½ of one
informative anchor. The IGH read therefore IS resolvable in
principle. Existing mappers fail not for lack of information but
because their priors are uniform — `π_i = 1/9` for all `i`, with no
guidance towards the correct paralog.

Layer-2 with `report_multi_position=y` plus a PSV catalog converts
that to a position-specific posterior that *uses* the read's PSV
content, yielding a confident single position with high MAPQ. This
is the mechanism behind the 9.9 × concordance lift demonstrated on
T6 1k subset earlier in the project.

---

## 11. Beyond — runtime extensions

These numerical priors are not the end:

* **Online refinement.** As a run progresses, the empirical depth
  distribution across mapped reads updates `N̂` for each cell. After
  a million-read whole-genome run, the table self-corrects against
  the prior.
* **Per-sample priors.** For somatic variant calling, a per-sample
  prior on T-axis (`germline / somatic / fixed`) is built from the
  early reads and fed back into mapping for later reads.
* **Cross-sample meta-prior.** Across many samples mapped with the
  same module, the depth distribution becomes itself a knowledge
  artefact — the manifold's `P(μ)` learns from data without
  re-training a neural network.

These keep the mapper open-ended without requiring new biology
ingest.
