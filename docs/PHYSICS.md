# LLmap Physics — formal math

## 1. The WaveCollapse update

For each read `r` with sparse candidate bucket set `B(r)` (top-K, default K=10 shrinking with level):

```
P_{t+1}(b|r) = (1-γ) · P_t(b|r) + γ · Z⁻¹(r,t) · F(b|r,t)

F(b|r,t) = L(r|b) · λ_t(b) · π_AI(b|r) · π_bio(b) · Σ_{b'∈N(b)} K(b,b') · P_t(b'|r)

Z(r,t) = Σ_{b ∈ B(r)} F(b|r,t)
```

| Symbol | Meaning |
|---|---|
| `P_t(b\|r)` | probability that read r is assigned to bucket b at iteration t |
| `L(r\|b)` | sequence likelihood (path-integrated over alignment trajectories — Forward algorithm, not Viterbi) |
| `λ_t(b)` | coverage prior at iteration t, `λ_t(b) = Σ_r P_t(b\|r)` |
| `π_AI(b\|r)` | AI prior: cosine similarity between read embedding and bucket embedding |
| `π_bio(b)` | biology prior from Claude index-build session (sparse; default 1.0) |
| `K(b,b')` | coverage-coupling kernel (Gaussian over genome distance; learned at each level) |
| `N(b)` | bucket neighbors at the current level |
| `γ` | platform-specific damping (∝ 1/T2; HiFi 0.05, ONT 0.3, Illumina 0.02) |
| `Z(r,t)` | normalizer over `B(r)` |

## 2. Path-integral semantics

`L(r|b)` is **not** the max-likelihood alignment score. It is the **sum over all alignment trajectories**:

```
L(r|b) = Σ_paths exp(-S[path] / kT)
```

where `S[path]` is the action functional accumulated along the alignment path (gap-open + gap-extend + mismatch costs integrated). This is implemented via the Forward algorithm on the pairHMM, not Viterbi. The cost is one matmul-pass per read per candidate, GPU-parallel.

**Why this matters**: max-likelihood alignment commits to one path. Path-integral keeps all paths weighted, which is what makes WaveCollapse lossless at the per-read level — we never discard alternative explanations.

## 3. Symmetry breaking

Sequence-identical paralog copies (e.g. IGHG4 canonical vs. dup-copy with identical exonic sequence) are **degenerate eigenstates** of `L(r|b)`. No per-read evidence breaks the degeneracy.

The coverage-coupling term `K(b,b') · λ_t(b')` is the **symmetry-breaking field**. Asymmetric coverage (e.g. mosaic duplication producing more reads from one paralog) breaks the symmetry collectively, even though no single read could.

**This is the mathematical mechanism behind LLmap's paralog-disambiguation USP.** It is also why per-read disambiguation of sequence-identical exons remains impossible (no symmetry-breaker per read) — only the locus-level coverage CNV recovers in that case.

## 4. Renormalization group flow

The bucket pyramid L0 → L1 → L2 → L3 implements an RG flow. Each level is a coarse-graining of the next:

- L0 (~1k buckets): chromosomes, repeat families, pseudogene clades
- L1 (~600 buckets): 5 MB windows
- L2 (~60k buckets): 50 kb windows
- L3 (continuous): exact position

Reads start at L0. Collapsed reads (`max_b P > τ`) drop out — they have reached a stable mapping. Non-converged reads expand their candidate set to the next level. The fixed points of the RG flow are the stable per-level mappings.

Operators at each level (the `K(b,b')` kernels, the bucket embeddings) are **learned**, not hand-designed. They are trained at Phase 3 of the build via end-to-end optimization on synthetic data with known ground truth.

## 5. Decoherence-T2 → platform damping

Sequencing error is environmental noise that destroys read–reference coherence. Mapped to a quantum-decoherence parameter T2:

| Platform | Q-score | T2 (qualitative) | γ |
|---|---|---|---|
| PacBio HiFi | Q30+ | long | 0.05 |
| ONT R10 | Q15-20 | medium | 0.3 |
| Illumina | Q35+ | very long | 0.02 |

`γ` is the EM damping coefficient. High γ means each iteration commits more aggressively (less smoothing across iterations) — appropriate for noisy data where individual updates carry less signal. Low γ means slow, careful convergence — appropriate for high-quality data where every iteration adds reliable information.

**Per-platform damping is theory-derived, not heuristic-tuned.** This is one of the concrete consequences of taking the physics seriously.

## 6. Collapse criterion

A read `r` is `Mapped` when:

```
max_{b ∈ B(r)} P_t(b|r) > τ
```

Default `τ = 0.99`. The read drops out of the active solver. Its bucket assignment is final at the current level; only the residual WFA2 extension runs on Level 3.

After `max_iter` (default 25) iterations, a read that has not collapsed becomes `Tentative` with its full probability distribution preserved in `tentative_targets`. **No read is ever forced to a primary.**

## 7. Stage 1: Self-Interference (entanglement-like coupling)

Before reads project to the reference, they interfere with each other. Concretely:

- All reads embed via the foundation model (Caduceus-Ph distilled)
- FAISS-GPU sparse k-NN over the embedding space
- Leiden community detection on the resulting similarity graph
- Each cluster runs its own Self-WaveCollapse (intra-cluster EM)
- Cluster representatives emerge

This is mathematically analogous to **entanglement formation**: reads from the same locus develop coupled probability distributions, distinguishable from independent reads via collective coherence. Reads from sequencing errors do not entangle (random noise does not cluster), so they remain as independent low-coupling members.

The result: ~100M raw reads → ~1M coherent representatives. Stage 2 then projects representatives onto the reference, with members inheriting via cheap intra-cluster delta-correction. Compute scales with cluster count, not read count.

## 8. Why this is not just "EM with extras"

Kallisto/Salmon do EM over read-to-transcript assignments. WaveCollapse differs in three ways:

1. **Hierarchical coarse-to-fine** (RG flow), not flat
2. **AI priors** (Foundation-Model embeddings) and **biology priors** (Claude) — not uniform
3. **Path-integral likelihood**, not max-likelihood — preserves alternative explanations as a feature

The combination is novel. To our knowledge, no production mapper implements all three.
