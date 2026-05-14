# LLmap

A sequence mapper for long and short reads against any reference (genome, transcriptome, pangenome, custom). LLmap targets workflows where **paralogs, pseudogenes, and segmental duplications matter** — IGH, MHC, NPHP1, CLN3, 22q11.2 — while remaining drop-in compatible for routine WGS / RNA-seq.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status: V1.0 Ready](https://img.shields.io/badge/Status-V1.0%20Ready-brightgreen.svg)](#status)
[![Algorithm: WaveCollapse](https://img.shields.io/badge/Algorithm-WaveCollapse-purple.svg)](#algorithm)
[![Domain](https://img.shields.io/badge/domain-losslessmap.com-green.svg)](https://losslessmap.com)

---

## What LLmap does differently

Standard mappers (minimap2, Winnowmap2, BWA-MEM2, strobealign) share three structural failures for paralog-rich loci:

1. **They force a single primary alignment.** For IGH / MHC / SD-regions this means most reads end up at `MAPQ=0` or are assigned arbitrarily, discarding recoverable information.
2. **No probabilistic output.** Uncertainty is collapsed at write-time; downstream tools cannot recover or propagate it.
3. **No biology-aware priors.** Every reference position has equal a priori weight, ignoring known SD regions, recurrent-NAHR loci, and pseudogene families.

LLmap addresses all three at the algorithmic level. Concretely:

- **Lossless output.** Every read produces a probabilistic record. A read whose origin cannot be uniquely determined is not silently dropped to `MAPQ=0` — it is reported as a probability distribution over candidate buckets. Downstream tools may collapse this to a single primary; LLmap itself never discards information.
- **Iterative EM over a hierarchical bucket pyramid (WaveCollapse).** Read placement is updated by an EM iteration that combines sequence likelihood, coverage prior, foundation-model embedding similarity, and biology priors.
- **Biology priors at the index level.** Reference annotation (SD regions, paralog families, recurrent NAHR loci) enters the score functional as `π_bio(b)` so the algorithm starts from the right prior rather than discovering known structure from scratch on every run.

---

## Algorithm

A read is represented as a probability distribution `P(b|r)` over buckets `b` in a hierarchical pyramid `L0 → L1 → L2 → L3`. The EM update rule:

```
P_{t+1}(b|r) = (1-γ) P_t(b|r) + γ · Z⁻¹ [
    L(r|b) · λ_t(b) · π_AI(b|r) · π_bio(b) · Σ_{b'∈N(b)} K(b,b') · P_t(b'|r)
]
```

Terms:

| Term | Meaning |
|---|---|
| `L(r|b)` | Sequence likelihood, path-integrated over alignment trajectories (Forward algorithm, not Viterbi) |
| `λ_t(b)` | Coverage prior — global signal that lets reads inform each other |
| `π_AI(b|r)` | Foundation-model embedding similarity |
| `π_bio(b)` | Biology prior from Claude-generated reference annotation |
| `K(b,b')` | Kernel over neighboring buckets at the same level |
| `γ` | Platform-specific damping (HiFi ≪ ONT) |

A read converges to a determined bucket when `max_b P_t(b|r) > 0.99`. Non-converged reads remain probabilistic in the output.

The formal mapping to path-integral semantics, renormalization-group flow over the bucket pyramid, and platform-dependent decoherence is documented in [docs/PHYSICS.md](docs/PHYSICS.md). This framing is mathematically equivalent to a hierarchical EM; it is used here as an intuition pump for the algorithm structure.

### Pipeline

```
INPUT FASTQ
    │
    ▼
[Foundation Embedder]  (distilled Caduceus, GPU-batched at ~10µs/read)
    │
    ▼
Reference WaveCollapse
  Bucket pyramid L0 → L1 → L2 → L3
  EM iteration with coverage coupling and biology priors
  Collapse-dropout per level
  WFA2 gap-affine extension on residual hard reads
    │
    ▼
Output
  ├─→ BAM / SAM       (samtools / IGV / bcftools-compatible)
  └─→ Parquet         (full lossless: read × bucket × probability)
```

For very large inputs (≥100M reads), an optional embedding-based pre-clustering pass deduplicates near-identical reads down to a smaller set of representatives that are then aligned, with members inheriting the rep alignment plus a delta-correction. This is a throughput optimization, not part of the core algorithm — enable with `--precluster`.

---

## Foundation models and Claude agent

LLmap loads pre-trained foundation models (Caduceus, Evo, Nucleotide Transformer, DNABERT-2) as frozen feature extractors via ONNX Runtime. They contribute the `π_AI(b|r)` term — semantic similarity between read embedding and bucket centroid.

Claude (Anthropic) is invoked as a tool-using agent — not as a per-read voter, which would be prohibitively expensive — at four predefined points:

| Session | Trigger | Role | Default cost / sample |
|---|---|---|---|
| **A: Index-Build** | Once per reference | Generates `biology_prior.json` (SD regions, paralog families, recurrent NAHR loci) using `bedtools`, UCSC tracks, custom preprocessors | $5 (amortized) |
| **B: Sample-Init** | Before each run | Reads FASTQ metadata, runs `seqkit stats` + `fastqc`, picks preset and parameters | $1 |
| **C: Diagnostic** | On EM stall | Inspects wave state; optionally generates a CUDA kernel hot-loaded via sandboxed `nvcc` | $5–15 (only when stalls occur) |
| **D: Reporter** | Post-run | Generates per-sample markdown diagnostic report, updates memoization cache | $2 |

Total typical per-sample cost: **~$3**, all asynchronous. Claude never blocks the alignment pipeline — its output is additive bias applied when ready. CUDA kernels generated by Session C are sandboxed (Bubblewrap network/PID/UTS/IPC isolation, static AST analysis, symbol allow-list) before being loaded.

Modes: `--llm {off, index-only, sample-aware, self-healing, research}`. Default: `sample-aware`. `--llm off` runs a deterministic CPU/GPU-only path with no external API dependencies.

---

## Performance targets

Targets relative to minimap2 on the same hardware. Empirical results are reported by the Phase 11 comparative benchmark campaign in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

| Metric | Target |
|---|---|
| Wallclock (HiFi WGS, 30× coverage) | ≤ 0.55× minimap2 |
| Wallclock (iso-seq with contam)    | ≤ 0.45× minimap2 |
| Peak RAM                           | ≤ 1.2× minimap2 |
| Read recall (uniquely-mappable)    | ≥ 99.5% (Phase 5 kill-switch) |
| Paralog assignment accuracy        | minimap2 + ≥ 10 percentage points |
| Read recovery in SD regions        | ≥ 2× usable alignments |

Speedup sources:

1. Speculative prefetch, kernel fusion, mixed precision (~20%)
2. Biology priors from Session A (~5–8% from better initial bucket layout)
3. Memoization across replicates / cohorts (additive 50%+ on similar samples)
4. Optional `--precluster` pre-deduplication on very large inputs (additional 30–50% when input has high redundancy, e.g. deep WGS)

---

## Quickstart

```bash
# Build index — Session A runs asynchronously in the background
llmap index \
  --fasta GRCh38.fa --gff gencode.v46.gff3 \
  --pangenome HPRC_r2_assemblies/*.fa \
  --llm sample-aware \
  -o GRCh38.llmap.idx

# Map (any data type, any reference)
llmap align \
  -i GRCh38.llmap.idx -r sample.fastq.gz \
  -x map-hifi --llm sample-aware \
  -o sample.bam -p sample.parquet

# Drop-in minimap2 replacement (auto-translates flags)
llmap-mm GRCh38.fa sample.fastq.gz -x map-hifi > sample.sam
```

Presets (`-x`): `map-hifi`, `map-ont`, `sr` (short-read), `splice` (RNA), `asm` (assembly-to-ref).

---

## Limitations

- **Sequence-identical exonic regions, per-read level.** Per-read disambiguation is impossible by information theory. LLmap returns `uninformative` explicitly. Coverage-based CNV inference (Stage 2 by-product) still recovers locus-level dup signal via summed depth.
- **5'-UTR-truncated iso-seq with PSV in UTR.** Information that was not sequenced cannot be recovered. Flagged explicitly, not silently dropped.
- **Phasing generation.** LLmap *consumes* HP-tags from upstream phasers (e.g. [pseudocaller](https://github.com/schlein-lab/pseudocaller)). It does not generate phasing.
- **Foundation model training.** We do not train Caduceus / Evo / Nucleotide Transformer from scratch. We distill from pre-trained public weights and train only small task heads on the IGHG4 + HPRC corpus.
- **Streaming alignment from sequencer.** V1.0 is batch-only.

---

## Status

V1.0 development complete. All 10 phases implemented and tested. 1,433 unit tests passing. GPU validation on Hummel-2 HPC pending. Phase 11 comparative benchmark campaign in progress.

| Phase | Status | Tests |
|---|---|---|
| 0 — Foundations + Synthetic Data       | complete | 56    |
| 1 — Foundation Model Integration       | complete | 110   |
| 2 — Pre-clustering pipeline (optional) | complete | 288   |
| 3 — Reference WaveCollapse             | complete | 493   |
| 4 — Classical Path + WFA2              | complete | 610   |
| **5 — Kill-Switch Validation** ★        | complete (CPU) | 684   |
| 6 — Dual Output (BAM + Parquet)        | complete | 770   |
| 7 — Claude Agent Integration           | complete | 923   |
| 8 — Performance Optimization           | complete | 1,092 |
| 9 — Single-Cell + Paralog Production   | complete | 1,249 |
| 10 — Production Readiness              | complete | 1,433 |
| 11 — Comparative Benchmark Campaign    | in progress | — |

See [CHANGELOG.md](CHANGELOG.md) for detailed feature list and [STATE.md](STATE.md) for build history.

---

## Citing LLmap

Manuscript in preparation. For now:

```bibtex
@software{llmap2026,
  author = {Schlein, Christian},
  title  = {LLmap: Lossless mapping with hierarchical EM for paralog-rich regions},
  year   = {2026},
  url    = {https://github.com/schlein-lab/LLmap}
}
```

---

## License

MIT. See [LICENSE](LICENSE).

---

## Documentation

- [LLmap_SPEC.md](LLmap_SPEC.md) — full V1.0 specification
- [docs/PHYSICS.md](docs/PHYSICS.md) — formal math (path integrals, RG flow, symmetry breaking)
- [docs/CLAUDE_AGENT.md](docs/CLAUDE_AGENT.md) — Claude agent integration
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — module-level architecture
- [docs/DATA_MODEL.md](docs/DATA_MODEL.md) — AlignmentRecord, BucketPyramid, WaveState
- [docs/BENCHMARKS.md](docs/BENCHMARKS.md) — comparative benchmark results
- [docs/MODELS.md](docs/MODELS.md) — foundation model catalog

---

## Acknowledgments

- **minimap2** (Heng Li) — the reference we measure against, and learn from
- **WFA2-lib** (Marco-Sola et al.) — gap-affine extension
- **FAISS** (Meta) — GPU approximate nearest-neighbor
- **Caduceus**, **Evo**, **Nucleotide Transformer**, **DNABERT-2** — pre-trained foundation models
- **HPRC** consortium — pangenome data
- **Anthropic Claude** — tool-using agent integration
