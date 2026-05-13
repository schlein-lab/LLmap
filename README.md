# LLmap

> **Where reads see each other first.**
> Lossless because LLM. Wave-particle because physics. Lossless mapping because mathematics.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: V1.0 Bootstrap](https://img.shields.io/badge/Status-V1.0%20Autonomous%20Build-blue.svg)](#status)
[![Algorithm: WaveCollapse](https://img.shields.io/badge/Algorithm-WaveCollapse-purple.svg)](#the-algorithm-wavecollapse)
[![AI: Claude](https://img.shields.io/badge/AI-Claude%20(Anthropic)-orange.svg)](#claude-as-a-tool-using-agent)
[![Domain](https://img.shields.io/badge/domain-losslessmap.com-green.svg)](https://losslessmap.com)

---

## What is LLmap?

LLmap is a **quantum-inspired, AI-augmented sequence mapper** for long and short reads against any reference (genome, transcriptome, pangenome, custom). It is engineered to replace minimap2 in workflows where **paralogs, pseudogenes, and segmental duplications matter** — IGH, MHC, NPHP1, CLN3, 22q11.2 — while remaining **drop-in compatible** for routine WGS/RNA-seq use.

The name carries two meanings, and we mean both:

| LL | Meaning | What it gives you |
|---|---|---|
| **L**ossless**L**map | No read is silently dropped. Every read produces a probabilistic record. | Mathematically guaranteed via WaveCollapse — the algorithm never forces a measurement it cannot justify. |
| **LL**M-Map | Claude (Anthropic) is a first-class architectural component, not an afterthought. | Index-time biology priors, runtime self-healing CUDA codegen, sample-aware presets. |

---

## Why LLmap exists

Standard mappers (minimap2, winnowmap2, BWA-MEM, strobealign) share three structural failures for paralog-rich loci:

1. **They force a single primary alignment** — for IGH/MHC/SD-regions this means >99% of reads end up `MAPQ=0` or arbitrarily assigned, destroying recoverable information.
2. **Reads do not inform each other** — collective signal that disambiguates paralogs is thrown away by per-read independent mapping.
3. **No biology-aware priors** — every position in the reference has equal a priori weight, ignoring decades of accumulated knowledge about SD regions, recurrent-NAHR loci, and pseudogene families.

LLmap fixes all three at the algorithmic level.

---

## The Algorithm: WaveCollapse

A read is **not a point** to be located. It is a **probability wave** over a hierarchical bucket space. The wavefunction evolves under a Hamiltonian composed of:

- Sequence likelihood `L(r|b)` — what minimap2 computes, but path-integrated over all alignment trajectories
- Coverage prior `λ(b)` — the global signal that reads communicate to each other
- AI prior `π_AI(b|r)` — semantic similarity from a frozen foundation model
- Biology prior `π_bio(b)` — Claude-generated annotation of the reference (Index-Build agent)

Update rule:

```
P_{t+1}(b|r) = (1-γ) P_t(b|r) + γ · Z⁻¹ [
    L(r|b) · λ_t(b) · π_AI(b|r) · π_bio(b) · Σ_{b'∈N(b)} K(b,b') · P_t(b'|r)
]
```

Reads **collapse** to a determined bucket only when `max_b P_t(b|r) > 0.99`. Non-converged reads stay probabilistic in the output. **No read is ever forced to a primary it cannot justify.**

### Two-Stage Pipeline

```
INPUT FASTQ
    │
    ▼
[Foundation Embedder]  ◄── Caduceus-Ph distilled, GPU-batched at 10µs/read
    │
    ▼
┌─ STAGE 1: SELF-INTERFERENCE ────────────────────────────────────┐
│  Reads inform each other BEFORE projecting to the reference.    │
│                                                                 │
│  FAISS-GPU sparse k-NN  →  Leiden clusters  →  intra-cluster    │
│  Self-WaveCollapse  →  ~1M coherent representatives from        │
│  ~100M raw reads                                                │
└─────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─ STAGE 2: REFERENCE WAVECOLLAPSE ───────────────────────────────┐
│  Cluster reps (not raw reads) project onto reference.           │
│                                                                 │
│  Bucket pyramid L0 → L1 → L2 → L3                              │
│  EM iteration with coverage coupling                            │
│  Collapse-dropout per level                                     │
│  WFA2 extension on residual hard reads                          │
└─────────────────────────────────────────────────────────────────┘
    │
    ▼
[Member Propagation]  ◄── cluster members inherit + cheap delta-correction
    │
    ▼
DUAL OUTPUT
    ├─→ BAM-Compat   (samtools/IGV/bcftools drop-in)
    └─→ Probabilistic-Parquet  (full lossless: read × bucket × probability)
```

---

## Physics, not metaphor

The wave-particle analogy is **mathematically isomorphic**, not branding:

| Quantum mechanics | LLmap WaveCollapse |
|---|---|
| Wavefunction ψ(b) | Read probability vector P(b\|r) |
| Hamiltonian | Score functional (Likelihood + Coverage + AI + Biology) |
| Path integrals | Sum over alignment trajectories (Forward algorithm, not Viterbi) |
| Decoherence-T2 | Platform-specific damping γ (HiFi ≪ ONT) |
| Symmetry breaking | Coverage-coupling resolving sequence-identical paralogs |
| Renormalization-group flow | Hierarchical bucket pyramid L0→L1→L2→L3 |
| Entanglement | Read-cluster coupling (Stage 1 Self-Interference) |
| Measurement/collapse | Convergence threshold τ = 0.99 |

See [docs/PHYSICS.md](docs/PHYSICS.md) for the full mapping.

---

## Claude as a Tool-Using Agent

Claude does not act as a per-read voter (which would be prohibitively expensive). Claude is a **tool-using agent** with `Bash`, `Read/Write`, `WebFetch`, and `CUDA Codegen` access, invoked in **four deep sessions** per analysis:

| Session | Trigger | What Claude does | Cost / sample |
|---|---|---|---|
| **A: Index-Build** | Once per reference | Runs `bedtools`, fetches UCSC SD tracks, writes custom Python preprocessors, generates `biology_prior.json` | $5 (amortized) |
| **B: Sample-Init** | Before each run | Reads FASTQ metadata, runs `seqkit stats` + `fastqc`, picks preset + parameters | $1 |
| **C: Diagnostic** | On EM stall | Dumps wave-state, investigates, **writes custom CUDA kernels** if needed, hot-loads via sandboxed compile | $5-15 (only if stalls) |
| **D: Reporter** | Post-run | Generates markdown diagnostic report, updates memoization cache | $2 |

**Total default cost: ~$3 / sample, all async.** Claude never blocks the GPU pipeline — its output is **additive bias** injected when ready.

Modes: `--llm {off, index-only, sample-aware, self-healing, research}`. Default: `sample-aware`.

> **LLmap is the first production bioinformatics tool that lets Claude generate CUDA kernels in the algorithmic hot-path.**

---

## Performance Targets

| Metric | Target vs minimap2 |
|---|---|
| Wallclock (HiFi WGS, 30× coverage) | ≤ 0.55× (~45% speedup) |
| Wallclock (iso-seq with contam) | ≤ 0.45× (~55% speedup) |
| Peak RAM | ≤ 1.2× |
| Read recall (uniquely-mappable) | ≥ 99.5% |
| Paralog assignment accuracy | minimap2 + ≥ 10 percentage points |
| Read recovery in SD regions | ≥ 2× more reads with usable alignment |

Speedup comes from (multiplicative, not additive):

1. **Self-Interference pre-clustering** (~35% — Stage 1 reduces 100M reads to ~1M coherent reps)
2. **Speculative prefetch + kernel fusion + mixed precision** (~20%)
3. **Claude-generated biology priors** (~5-8% from better initial bucket layout)
4. **Memoization across similar samples** (additive 50%+ on replicates/cohorts)

---

## Quickstart

```bash
# Install (target)
brew install llmap

# Build index — runs Claude Session A asynchronously in background
llmap index \
  --fasta GRCh38.fa --gff gencode.v46.gff3 \
  --pangenome HPRC_r2_assemblies/*.fa \
  --llm sample-aware \
  -o GRCh38.llmap.idx

# Map any data type, any reference
llmap align \
  -i GRCh38.llmap.idx -r sample.fastq.gz \
  -x map-hifi --llm sample-aware \
  -o sample.bam -p sample.parquet

# Drop-in minimap2 replacement (auto-translates flags)
llmap-mm GRCh38.fa sample.fastq.gz -x map-hifi > sample.sam
```

---

## What LLmap is NOT

We are honest about limitations:

- **Sequence-identical exonic regions, per-read level**: per-read disambiguation remains impossible. LLmap returns `uninformative` explicitly. Coverage-based CNV inference (Stage 2 by-product) still recovers locus-level dup signal via summed depth.
- **5'-UTR-truncated iso-seq with PSV in UTR**: cannot recover information that was not sequenced. Flagged explicitly, not silently dropped.
- **Phasing generation**: LLmap *consumes* HP-tags from upstream phasers (e.g. [pseudocaller](https://github.com/schlein-lab/pseudocaller)). It does not generate phasing.
- **Foundation model training**: we do not train Caduceus/Evo/Nucleotide-Transformer from scratch. We distill from pre-trained public weights and train only small task heads on our IGHG4 + HPRC corpus.
- **Streaming alignment from sequencer**: V1.0 is batch-only.

---

## Status

LLmap is in **V1.0 autonomous build**. The implementation is being driven by a 96-hour Claude-Code autonomous run on Hummel-2 HPC, beginning May 2026. Phase progress is tracked in [STATE.md](STATE.md) and live-updated by the autonomous driver.

| Phase | Status | Acceptance gate |
|---|---|---|
| 0 — Foundations + Synthetic Data | 🔨 in progress | unit tests + serialization |
| 1 — Foundation Model Integration | ⏳ | 100M reads embed < 2h on H100 |
| 2 — Stage 1 Self-Interference | ⏳ | 100M → 1M clusters in < 30 min |
| 3 — Stage 2 Reference WaveCollapse | ⏳ | synthetic IGH dup-fraction ±2% |
| 4 — Classical Path + WFA2 | ⏳ | CPU fallback runs end-to-end |
| **5 — KILL-SWITCH VALIDATION** ★ | ⏳ | recall ≥ 99.5% minimap2 — go/no-go |
| 6 — Dual Output (BAM + Parquet) | ⏳ | samtools compat + lossless invariant |
| 7 — Claude Agent Integration | ⏳ | biology prior gain + CUDA sandbox |
| 8 — Performance Optimization | ⏳ | ≤ 0.55× minimap2 wallclock |
| 9 — Single-Cell + Paralog Production | ⏳ | HG002 IGHG4 + Kinnex PBMC |

Live progress feed (autonomous commits): https://github.com/schlein-lab/LLmap/commits/main

---

## Citing LLmap

Manuscript in preparation. For now, please cite as:

```bibtex
@software{llmap2026,
  author = {Schlein, Christian},
  title  = {LLmap: Lossless, LLM-augmented, wave-particle sequence mapper},
  year   = {2026},
  url    = {https://github.com/schlein-lab/LLmap},
  note   = {V1.0 autonomous build, May 2026}
}
```

---

## License

MIT. See [LICENSE](LICENSE).

---

## Documentation

- [LLmap_SPEC.md](LLmap_SPEC.md) — full V1.0 specification (700+ lines, agent-feedable)
- [docs/PHYSICS.md](docs/PHYSICS.md) — formal math, path integrals, RG flow, symmetry breaking
- [docs/PHILOSOPHY.md](docs/PHILOSOPHY.md) — the LL double meaning, photon analogy, design ethos
- [docs/CLAUDE_AGENT.md](docs/CLAUDE_AGENT.md) — how Claude is integrated as a tool-using agent
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — module-level architecture
- [docs/DATA_MODEL.md](docs/DATA_MODEL.md) — AlignmentRecord, BucketPyramid, WaveState
- [docs/BENCHMARKS.md](docs/BENCHMARKS.md) — vs minimap2 / winnowmap / strobealign
- [docs/MODELS.md](docs/MODELS.md) — foundation-model catalog, distillation recipes

---

## Acknowledgments

LLmap stands on the shoulders of:

- **minimap2** (Heng Li) — the reference we measure against, and learn from
- **WFA2-lib** (Marco-Sola et al.) — gap-affine extension
- **FAISS** (Meta) — GPU ANN
- **Caduceus**, **Evo**, **Nucleotide Transformer**, **DNABERT-2** — public foundation models
- **HPRC** consortium — pangenome data
- **Anthropic Claude** — the first LLM to write production CUDA kernels in a bioinformatics tool

---

*Reads as photons. Genome as crystal. Mapping as decoherence. Every read accounted for.*
