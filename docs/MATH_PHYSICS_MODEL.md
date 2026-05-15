# LLmap — formal mathematical & physical model

This document is the working theory of the LLmap mapper. It explains *why*
the implementation is shaped the way it is, in language a physicist /
mathematician / statistician can verify and extend.

Mapping is approached as **Bayesian inference of read origin under
physical constraints**, not max-likelihood pattern matching.

> *A genome is mathematics and physics. We have long modelled it as biology
> with text-search. The text-search era is over.*

The rigid application of one global pattern-matching rule (one `k`, one
`identity_threshold`, one chain-score function) across every region of
every genome of every organism — the assumption made implicitly by every
production short- and long-read mapper — is **no longer adequate** for
the data we now collect (HiFi at 99.9 % accuracy, full-T2T assemblies,
T2T pangenomes, exhaustive population variant catalogs). LLmap is a
re-foundation: parameters vary per base, priors compose, knowledge is
plumbed into the inference rather than glued on afterwards.

---

## 0. The gap in existing mappers

Every classical mapper (BWA, BWA-MEM2, minimap2, winnowmap, pbmm2,
NGMLR, ...) solves

```math
\hat{p} \;=\; \arg\max_{p} \; P(\mathrm{read} \mid p)
```

with an implicit uniform prior `P(p) = const`. They ignore everything
the human community already knows about the reference:

* per-position population variant frequencies (gnomAD, dbSNP, dbVar, DGV, HPRC)
* clinical interpretations (ClinVar)
* segmental-duplication architecture (Bailey 2002, Eichler 2019, T2T 2022)
* known centromere & telomere structures (Willard, Hayden, T2T-CHM13)
* published per-locus biology in literature (PubMed)

The information exists. Mappers don't use it.

LLmap consumes it as **prior knowledge**:

```math
P(p \mid \mathrm{read}) \;\propto\; P(\mathrm{read} \mid p) \cdot P(p)
```

with `P(p)` decomposed into composable layers.

---

## 1. Layered Bayesian posterior

Let `R` be a read and `p` a candidate reference position. Define the
posterior

```math
\log P(p \mid R) \;=\;
  \log P(R \mid p)
  + \log P_{1}(p)
  + \log P_{2}(p)
  + \log P_{3}(p)
  + \log P_{4}(p)
  - \mu \, \mathcal{L}_{\mathrm{novelty}}(R, p)
  - \log Z
```

The four prior layers, by increasing specificity:

| Layer | Source | Determinism |
|---|---|---|
| `P_1` | Universal taxonomy — region-type classification from local features | deterministic |
| `P_2` | Per-locus database — coords + expected architecture of known SDs, centromeres, etc. | deterministic |
| `P_3` | Active LLM agent at checkpoints (tool-use, web fetch, sequence grep) | optional |
| `P_4` | Population variant priors — gnomAD / dbVar / ClinVar / HPRC / 1KG | deterministic |

`L_novelty` is a *negative-evidence* term that **rewards** the mapper for
flagging structures that are not in any layer — keeping the system open
to discovering new biology rather than silently collapsing onto
known-incorrect coordinates. See § 6.

---

## 2. Reads as waves; positions as eigenstates

We treat each read as a complex amplitude over reference positions:

```math
\psi_R(x) \;\in\; \mathbb{C}, \qquad
| \psi_R(x) |^2 \;=\; P(\text{origin}{=}x \mid R)
```

In a unique single-copy region the wavefunction has a sharp peak — the
read is *particle-like*. In a paralog cluster, centromere, or low-complexity
band it is spread across multiple plausible positions — *wave-like*. The
posterior of §1 sets the probability density, but the mapper's **internal
representation** is the full wavefunction, so it can carry interference
between candidates rather than committing too early.

### 2.1 Wavelength = effective k-mer scale

Anchors enter the mapper through minimizers. The effective minimizer
length `k` plays the role of the de-Broglie wavelength:

```math
\lambda(x) \;=\; c \cdot H_{\mathrm{loc}}(x)^{-1/2}
```

where `H_loc(x)` is the local Shannon entropy of 5-mers in a window
around `x`. High-information regions (high `H`) have small `λ` — sharp
localization. Low-information regions (low `H`) have large `λ` —
delocalised but still sensitive.

LLmap implements this directly: the `lambda_scale` field in
`ParamOverride` (declared in `src/annot/annot_types.h`) lets per-region
or per-base annotation set this scale.

### 2.2 Heisenberg-like uncertainty

Position resolution `Δx` and chain-confidence `Δp` are not independent:

```math
\Delta p \cdot \Delta x \;\gtrsim\; \hbar_{\mathrm{map}}
\;:=\; \frac{k \log_2 N}{4}
```

where `N` is the count of distinct k-mers in the reference. This is
not metaphor — it falls out of how a sparse minimizer index can
distinguish positions: the smaller the average `k`, the more anchors
per read, but each carries less information.

The constant is calibrated empirically per organism and stored in the
organism's `MODULE.md`. For human at `k=19`, `N ≈ 2.7×10⁹`, giving
`ℏ_map ≈ 150`. Concretely: a chain with sum-score 600 has position
uncertainty `≈ 0.25 bp`; a chain with sum-score 30 has uncertainty
`≈ 5 bp` and should not be reported as a single point — it should
remain a small wave.

---

## 3. Energy / Action formulation

Define the **alignment action** for a candidate alignment `I = (p, c, …)`:

```math
E(I) \;=\; -\log P(I \mid R)
\;=\; E_{\mathrm{seq}}(I) + \sum_{\ell=1}^{4} \beta_{\ell} \, E_{\ell}(I) + \mu \, E_{\mathrm{novelty}}(I)
```

Each `E_ℓ` is `-log P_ℓ` from § 1.

Mapping is then **minimum-action selection** under physical constraints:

```math
\hat{I} \;=\; \arg\min_{I} \; E(I)
```

This is exactly the chain-DP objective when only `E_seq` is present and
when the per-anchor weights are uniform. Adding the layers makes the
landscape *region-aware*: the saddle points between paralog wells are
raised by `E_2` (per-locus knowledge says "two distinct paralogs live
here, you can't tunnel between them cheaply"), and the floor of a known
polymorphism well is dropped by `E_4` (gnomAD says "this variant is
real, don't penalise its mismatches").

### 3.1 Why minimap2 fails in paralogs

In a paralog cluster the wells of `E_seq` are nearly degenerate
(`< 5%` score gap). A pure-`E_seq` minimiser hops between them
arbitrarily — what looks empirically as **79 % position concordance**
between LLmap and minimap2 on the IGH locus is just two tools
independently rolling dice with different RNGs.

`E_2` breaks the degeneracy: the per-locus database knows the four
IGHG paralogs are at fixed offsets within the IGH locus, and assigns
non-uniform priors based on PSV (paralog-specific-variant) markers
visible in the read. Mapping then becomes deterministic.

---

## 4. Genetic-algorithm dynamics on the wavefunction

The wavefunction is approximated by a finite **population** of
candidate alignments

```math
\mathcal{P}_t = \{ I_1, I_2, \ldots, I_n \}, \qquad
I_i = (p_i, \ell_i, \lambda_i, A_i, \mathcal{S}_i)
```

* `p_i` — candidate position
* `ℓ_i` — effective alignment length
* `λ_i` — wavelength (k-mer scale)
* `A_i` — amplitude (current probability mass)
* `S_i` — anchor set (which minimizers support this candidate)

### 4.1 Operators

* **Mutation** (drift):
  `p_i ← p_i + N(0, σ²(λ_i))` — diffusion proportional to local wavelength
* **λ-mutation** (zoom):
  `λ_i ← λ_i (1 ± ε)` — re-evaluate at a different scale
* **Crossover** (anchor pooling):
  `S_c = S_a ∪ S_b` when positions are compatible — let candidates share evidence
* **Selection** (tournament with diversity penalty):
  `F'(I) = F(I) − γ ⋅ \mathrm{density}(I, \mathcal{P})`
  prevents premature collapse onto one candidate

### 4.2 Fitness

`F(I) = −E(I)`. Selection pressure `T` ("temperature") is annealed across
iterations:

```math
T_t \;=\; T_0 \cdot e^{-t / \tau}
```

Early iterations: high `T`, diverse population, broad search.
Late iterations: low `T`, candidates collapse to local minima.

### 4.3 Collapse criterion

The wavefunction collapses to a single position iff

```math
\frac{\mathrm{Var}(p_i)}{\bar{E}^2} < \tau \quad\text{AND}\quad
\mathrm{eff}(\mathcal{P}) < n_{\min}
```

`eff(P) := 1 / Σᵢ Aᵢ²` is the effective population size — a measure of
how concentrated the amplitude is on one candidate. When collapsed,
LLmap emits a single SAM record. When not collapsed, it emits a
**wave**: multiple alternative alignments with amplitudes summing to
1, written as a `XW:Z:<list>` BAM tag.

This *failure-to-collapse* is biological signal, not failure — it
means the locus genuinely contains multiple equally good explanations
under all currently-known priors.

---

## 5. Stochastic online Bayesian update

Empirical statistics collected during a run feed back into the
priors:

```math
P_{\mathrm{posterior}}(\theta) \;\propto\; P_{\mathrm{prior}}(\theta) \cdot P(\text{empirical} \mid \theta)
```

Concretely: after `n` reads have been mapped, the empirical
identity-distribution in each annotation interval is computed. If the
prior `min_identity = 0.85` for this window but the observed median
is `0.78 ± 0.04`, the threshold is **updated** for subsequent reads in
the same run to `μ − 2σ ≈ 0.70`. Reads that were rejected in an
earlier pass can be re-evaluated in a second pass.

This is invisible to the user and fully deterministic given the read
order. It is implemented as a fifth layer `θ_stochastic` in the
ParamOverride composition (priority 4 in the interval-tree layer
stack, see `src/annot/annot_types.h`).

---

## 6. The novelty channel

A mapper that knows only what's in the databases will collapse novel
biology onto the wrong known position. LLmap reserves a non-zero
*novelty cost* term:

```math
\mathcal{L}_{\mathrm{novelty}}(R, p) \;:=\; - \log P_{\mathrm{novel}}(R, p)
```

where `P_novel` is a penalty against silently accepting an alignment
whose features are inconsistent with **all** of the four layers. When
this term dominates, LLmap emits a BAM record with **`XR:Z:novel_*`
tags**, rather than a single forced position:

* `XR:Z:novel_SV_candidate` — chain confident, but identity tracks a
  large indel pattern not in dbVar/gnomAD-SV → novel candidate SV
* `XR:Z:expanded_locus` — read length exceeds plausible window after
  multi-anchor consistency check → SD expansion candidate
* `XR:Z:divergent_haplotype` — mean identity in this annotation
  interval is systematically below prior expectations across multiple
  reads → novel haplotype in this sample
* `XR:Z:unknown_region` — no layer matched at all → report position
  but flag for human review

Existing mappers don't have this concept. They either accept the best
position with high MAPQ (false positive) or call the read unmapped
(false negative). Having a **third outcome** — *mapped but flagged
novel* — is the actionable biology signal.

---

## 7. Runtime database integration

Layers `P_1`, `P_2`, `P_4` are deterministic files. Layer `P_3` (the
active LLM agent) consults databases at runtime through a controlled
tool set:

| Tool | What it does | Cost guard |
|---|---|---|
| `region_lookup` | lookup in local `specific_loci/*.json` | free, ~1 µs |
| `psv_check` | lookup in local PSV catalog `.tsv` | free, ~1 µs |
| `variant_query` | lookup in local `.priors` table (gnomAD / dbVar / etc.) | free, ~10 µs |
| `local_grep` | grep the reference for a short sequence | free, ~1 ms |
| `web_fetch_ucsc` | HTTP GET against UCSC track server | rate-limited + cached |
| `bash` | whitelisted shell commands inside a sandbox | rate-limited |
| `pubmed_query` | NCBI E-utilities query for keyword / locus | rate-limited + cached |

These are invoked at well-defined **mapping checkpoints**:

* `AmbiguousChain` — top-N chain scores within 5 % of each other
* `UnknownRegion` — no layer-1/2 annotation applies
* `ParalogDisambiguation` — region flagged `require_psv_disambig`
* `SDExpansion` — chain length > window
* `NovelInsertion` — large unexplained soft-clip

The agent's response is **cached on disk** by content hash (read +
context + candidates) so that re-runs are cheap. After the first run
over a dataset, subsequent runs over the same dataset cost essentially
zero LLM calls.

This is also the path by which **new databases become consumable**: we
add a `tool_dbsnp_lookup`, a `tool_alphafold_struct_at_locus`, etc.,
and the agent learns to invoke them. The mapper's knowledge surface
grows continuously as we add data sources.

---

## 8. Modes of operation

| Flag | Behaviour |
|---|---|
| `--llm=off` | Layers 1, 2, 4 only — fully deterministic, fully offline |
| `--llm=auto` | + Layer 3 when available; silent fallback if not |
| `--llm=required` | Fail if Layer 3 not reachable — reproducible review runs |

Across all three modes, **the mathematical model is the same**. What
changes is whether `P_3(p)` is included in the posterior.

---

## 9. Computational complexity

| Phase | Per read | Notes |
|---|---|---|
| Seeding (minimizer query) | `O(k · m)` | `m` = minimizers in read |
| Layer 1 lookup | `O(log W)` | `W` = annotation intervals (≈ ref-size / 1kb) |
| Layer 2 lookup | `O(log L)` | `L` = number of specific loci |
| Layer 4 lookup | `O(log V)` | `V` = variant priors entries |
| Chain DP | `O(m²)` worst, `O(m · max_skip)` typical | unchanged from minimap2-style |
| GA candidates | `O(P · iters)` | `P` = population size (8 by default), `iters` = until collapse |
| Layer 3 (when fired) | `O(1)` calls + cache hit/miss | ~5 % of reads in typical run |

Wall-clock overhead vs minimap2 is dominated by Layers 1-2-4 lookups,
which are O(log) and cache-friendly. End-to-end overhead measured to
date is ~15-20 % on T6 — buying you a posterior-based Bayesian mapping
in exchange.

---

## 10. What this gives you that didn't exist before

* **Mapper output that carries uncertainty.** Wave-output reads,
  multi-position emission, novelty flags — actionable biology rather
  than forced single-position max-likelihood.
* **Region-precise per-base parameter manifold.** Every chain-DP
  parameter (k, w, max_occ, lambda_scale, identity_threshold, gap
  penalties) varies base-by-base. Mapping is no longer one global
  config; it is per-position physics with local knowledge.
* **Composability with knowledge.** Add a new DB → expose a new tool
  → mapper uses it. The knowledge surface is open-ended.
* **Honest novelty discovery.** When a read doesn't fit anything
  known, that's reported, not concealed by a forced position.

---

## 11. Organism-agnostic / pure-math mode

The four prior layers in § 1 are **optional**. Every layer is a `+ log Pᵢ`
term in the action — if no `Pᵢ` is provided, the corresponding term is
zero and the inference still works.

A useful corollary: **LLmap can map sequences for which no biological
knowledge exists at all.** Synthetic oligo pools, designed CRISPR
constructs, sequence from a novel organism we have never seen, fictive
test data, alien DNA — the layered priors degrade gracefully to the
sequence-only term `E_seq`. What you lose is the region-specific
refinement; what you keep is

* the **physical** treatment (wave-particle formulation, Heisenberg
  bound, GA dynamics, energy minimisation),
* the **mathematical** treatment (Bayesian posterior, multi-position
  wave emission, novelty channel),
* the **online stochastic update** of empirical priors during the run.

The field has long conflated *sequence alignment* with *biological
sequence alignment*. They are not the same. The first is a problem in
discrete mathematics and statistical physics; the second is the first
plus a particular set of priors. LLmap separates the two cleanly:

```
mapping(sequence)            ← Pure-math mode, works on anything
mapping(sequence, organism)  ← Add P_taxonomy + P_locus from organism module
mapping(sequence, organism, --llm=auto)
                              ← Plus P_agent at checkpoints
mapping(sequence, organism, --llm=auto, --variants priors)
                              ← Plus P_variant from gnomAD / dbVar / etc.
```

Each step adds capability; none is required. The base case is
sequence-only and is still mathematically principled.

---

## 12. On-premise deployment

A "fully on-premise" LLmap operates with these constraints:

* No network access at runtime — no live API calls, no web fetches
* All prior knowledge present as local files
* No cloud LLM (Layer 3 is `--llm=off` or `--llm=local` with a
  self-hosted model)
* Reproducible: same input + same prior files → bit-identical output

This mode is the default for clinical / forensic / IP-sensitive sites.
LLmap must support it natively.

### 12.1 What ships locally

| Component | On-premise form |
|---|---|
| Reference genome | local FASTA (already standard) |
| Layer 1 (taxonomy) | `knowledge/organisms/<org>/{classifier_rules.json, regions/*.json}` — KB-scale |
| Layer 2 (specific_loci) | `knowledge/organisms/<org>/specific_loci/**.json` — MB-scale |
| Layer 4 (variants) | pre-ingested `.priors` from gnomAD / dbVar / ClinVar / HPRC etc. — typically GB-scale |
| Layer 3 (agent) | one of three options: `--llm=off` (skipped), `--llm=local <endpoint>` (talks to an on-prem LLM gateway), or `--llm=cached` (reuses agent decisions cached from a prior run) |

### 12.2 Air-gapped operation

In a fully air-gapped environment:

* The site obtains the LLmap binary plus a release-tagged `knowledge/`
  tarball and any required `.priors` files via approved channels (USB,
  vetted artefact server).
* `llmap align --llm=off` is the routine command.
* No DNS resolution, no outbound traffic, no telemetry.
* The novelty channel still works — it is intrinsic to the math, not
  to any external service.

### 12.3 Hybrid on-prem / agent

Some sites have an on-prem LLM gateway (e.g., Anthropic on-prem, a
locally-hosted Llama). `--llm=local --endpoint <url>` routes the
agent through that gateway; nothing else changes. The cache layer
makes the agent essentially free on second runs.

### 12.4 Reproducibility audit trail

Every record in the output BAM carries provenance tags:

* `XL:Z:1234567` — local knowledge SHA-256 (matches a specific release of `knowledge/`)
* `XK:Z:0c46…` — variant priors SHA-256 (matches a specific build of the priors file)
* `XA:Z:on|off|local-<endpoint>` — agent mode at the time of mapping
* `XC:Z:<decision_id>` — if Layer 3 fired, the cached decision identifier

A reviewer can later regenerate exact-byte output given the same
inputs and the recorded knowledge SHAs.

---

## 13. Chromosomal positions as a multi-dimensional manifold

A chromosomal "position" in the classical mapper sense is a 1-D
integer: `(chr, base)`. That representation is *too small* to describe
real biology — two reads that both have `chr14:105580000` may belong
to different paralogs, different haplotypes, different somatic clones,
different ancestral haplotype blocks, different gene-conversion
events. The classical 1-D `pos` collapses all of that into a single
number.

LLmap operates on a richer **position manifold**

```math
\mathcal{M} \;=\; \mathbb{R} \;\times\; \mathcal{A} \;\times\; \mathcal{P} \;\times\; \mathcal{S} \;\times\; \mathcal{T} \;\times\; \mathcal{F}
```

with axes:

| Axis | Meaning | Source of values |
|---|---|---|
| `R` | Reference 1-D position `(chr, pos)` | the assembly itself |
| `A` | Allele / haplotype index | population variant catalogs (gnomAD, HPRC) |
| `P` | Paralog identity inside a multi-copy family | `specific_loci/*` + PSV catalogs |
| `S` | Structural context (within centromere, SD block, telomere, NUMT, etc.) | Layer 1 + Layer 2 annotation |
| `T` | Temporal axis (germline vs somatic vs evolutionary fixed) | ClinVar, COSMIC, per-sample history |
| `F` | Population frequency band (singleton, rare, common, fixed) | gnomAD AF distribution |

Every point `μ ∈ M` carries a precomputed local prior `P(μ)`
("blueprint"): how often is something here, what should reads look
like, what variants are expected, what paralog signature should
dominate.

The **mapping problem becomes a projection problem**:

```math
\text{read} \;\xrightarrow{\text{embed}}\; \mathbf{f} \in \mathcal{F}_{\text{feat}}
\;\xrightarrow{\text{collapse}}\; \mu \in \mathcal{M}
```

* `embed` maps the read into a feature space (minimizer signature,
  k-mer rarity, mismatch pattern, length, soft-clip profile).
* `collapse` is the WaveCollapse step: starting from a probability
  cloud `ψ(μ)` over the manifold, apply the four prior layers in
  succession; each layer is a projection operator that contracts the
  cloud along one or more axes:

```math
\psi_0(\mu) \;\xrightarrow{\hat{\Pi}_1}\;
\psi_1(\mu) \;\xrightarrow{\hat{\Pi}_2}\;
\psi_2(\mu) \;\xrightarrow{\hat{\Pi}_3}\;
\psi_3(\mu) \;\xrightarrow{\hat{\Pi}_4}\;
\psi_4(\mu)
```

The cloud either collapses to a single `μ̂` (deterministic mapping)
or remains a wave (multi-position output with amplitudes summing to
one). The novelty channel reports the residue when the read's
features lie far from any well of `P(μ)`.

### 13.1 The blueprint as precompiled prior

The `M` manifold's prior `P(μ)` is **precompiled** from the layered
knowledge — that is the blueprint:

* **Density blueprint** — expected read coverage per `(R, A, P, S)`
  cell, learned from the population coverage distribution. The
  mapper uses this to flag windows that are unexpectedly empty or
  unexpectedly oversaturated.
* **Pattern blueprint** — expected mismatch/soft-clip patterns. A
  read showing the canonical alpha-satellite higher-order-repeat
  shift is *expected* in centromeres and not penalised; the same
  pattern in coding regions triggers the novelty channel.
* **Reunion blueprint** — which reads should *cluster* together at
  what `μ`. Reads from the same paralog `P` should pool; reads from
  different paralogs `P, P'` should stay separated even when their
  `R` is similar.

### 13.2 Three uses of the blueprint at runtime

1. **Pre-sorting.** Before chain DP, each read's coarse manifold cell
   is determined from the minimizer signature alone. Reads landing
   in the same `(R-bin, S-bin)` cell are batched and processed
   together; their decisions share regional context. This is what
   makes the per-base parameter manifold tractable — we don't
   recompute Layer-1/2/4 priors per read, we precompute them per
   cell and reuse.

2. **Pattern recognition.** A *cluster* of reads at the same `μ`
   that exhibit the same off-prior signature is itself a signal: a
   private SV in this sample, a clonal somatic event, a hidden
   paralog copy not in the reference. The mapper reports this
   cluster as a coherent finding rather than as N independent reads.

3. **Targeted reunification.** Reads at neighbouring `μ` (same `P`,
   adjacent `R`) are explicitly grouped for assembly / phasing
   downstream. The mapper's output is a structured manifold-aware
   BAM where the `XP` tag carries the paralog index and the `XT` tag
   carries the temporal-context guess (germline / somatic / fixed).

### 13.3 Mathematical statement

The manifold-aware mapping is

```math
\hat{\mu}(R) \;=\; \arg\max_{\mu \in \mathcal{M}}
\; \log P(R \mid \mu) + \sum_{\ell=1}^4 \beta_\ell \log P_\ell(\mu)
- \mu_{\text{nov}} \mathcal{L}_{\text{novelty}}(R, \mu)
```

subject to the Heisenberg-like bound of § 2.2 — `μ̂` is **emitted as
a point** only if the wavefunction has collapsed; otherwise the
emission is the wavefunction itself.

This is the rigorous statement of what the field has been calling
"alignment". Once one writes it this way, the existing-mapper
behaviour (uniform prior, point-output everywhere, no novelty
channel) is visibly a degenerate special case — not best practice.

---

## 14. Open questions

* Calibration of the `β_ℓ` mixing coefficients. Currently fixed; eventually learned from held-out gold-standard datasets per organism.
* Calibration of `ℏ_map`. Currently derived analytically from `k` and `N`; could be empirical.
* Whether the wavefunction representation is ever provably worse than a single MAP estimate. For confident unique mappings, the population should collapse instantly and the cost should approach zero.
* Multi-organism mode for chimaeric samples (metagenomics, xenografts) — each contig carries its own organism module; the prior layers compose per contig.
