# LLmap — The Lossless, LLM-Augmented, Wave-Particle Mapper

**Tagline**: *Where reads see each other first. Lossless because LLM. Wave-particle because physics.*

**Status**: Greenfield V1.0 specification. Feed to a coding agent and start.

**Domain**: [losslessmap.com](https://losslessmap.com)

---

## 0. Brand & Philosophy

### The double meaning of "LL"
- **LosslessMap** — every input read produces a probabilistic record; no read is silently dropped, ever. Lossless **by construction** via WaveCollapse (the algorithm never forces a measurement it cannot justify).
- **LLM-Map** — Claude (Anthropic API) is a first-class architectural component, not an afterthought. Claude is a *tool-using agent* (Bash, Code, Fetch, Read), not a per-read voter.

### The physics
Reads are not points to be located — they are **probability waves over a hierarchical bucket space**. The wavefunction evolves under a Hamiltonian composed of sequence likelihood, coverage coupling, AI prior, and Claude-generated biology prior. Reads collapse to a determined position only when the probability mass concentrates above threshold. Unconverged reads stay probabilistic in the output.

This is *not metaphor*. The math is isomorphic to many-body quantum-inspired EM with path-integral semantics over alignment trajectories. See `docs/PHYSICS.md`.

### Non-negotiable principles
1. **Lossless by construction** — every read produces an `AlignmentRecord` with status ∈ {Mapped, Tentative, Unmapped}. Unmapped records carry a `RejectionReason`.
2. **Speed never sacrificed to AI** — Claude runs asynchronously in background sessions. Foundation-model inference is GPU-batched. Claude output is *additive bias* injected when ready; pipeline never blocks on it.
3. **GPU+AI are V1.0 architectural constraints** — without them the tool runs ~10× slower with reduced paralog accuracy, but remains functional (CPU classical fallback).
4. **Phase-5 is an empirical kill-switch** — if WaveCollapse cannot beat minimap2 on a recall gate, the project pivots or stops before Phase 6-9 build out.
5. **Drop-in compatibility** — BAM-Compat output by default; `samtools`/`bcftools`/IGV must work without modification.

---

## 1. Why this exists

Current mappers (minimap2, winnowmap, BWA-MEM, strobealign) share three structural failures:

1. **Forced measurement** — a single primary alignment must be chosen at MAPQ-time. For IGH/MHC/SD/paralog loci this destroys >99% of recoverable information.
2. **No read-to-read information sharing** — each read is mapped independently. The collective signal that disambiguates paralogs (coverage coupling, cluster coherence) is lost.
3. **No biology-aware priors** — every position in the genome is treated as having equal a priori mapping likelihood. Known SD-regions, pseudogene families, recurrent-NAHR loci get no special treatment.

LLmap fixes all three at the algorithmic level (WaveCollapse) and adds an LLM-agent layer for biology-aware reasoning, anomaly detection, and runtime self-improvement.

---

## 2. Algorithm: WaveCollapse

### 2.1 Formal update

For each read `r` with candidate bucket set `B(r)` (sparse, top-K):

```
P_{t+1}(b|r) = (1-γ) P_t(b|r) + γ · Z⁻¹ [
    L(r|b) · λ_t(b) · π_AI(b|r) · π_bio(b) · Σ_{b'∈N(b)} K(b,b') · P_t(b'|r)
]
```

where:
- `L(r|b)` — sequence likelihood (minimizer-seed + WFA2 alignment score, exp-transformed)
- `λ_t(b)` — coverage prior at iteration t (∑_r P_t(b|r))
- `π_AI(b|r)` — AI prior: cosine similarity between read embedding and bucket embedding (Foundation Model)
- `π_bio(b)` — biology prior from Claude Index-Build Session (`biology_prior.json`)
- `K(b,b')` — coverage-coupling smoothing kernel (Gaussian over genome distance, learned)
- `N(b)` — bucket neighbors (1MB radius at L1, 50kb at L2, etc.)
- `γ` — platform-specific damping constant (∝ Q-score, ~ Decoherence-T2⁻¹)
- `Z` — normalization over `B(r)`

### 2.2 Collapse criterion

Read `r` collapses at iteration t if:
```
max_b P_t(b|r) > τ_collapse  (default τ = 0.99)
```
Collapsed reads drop out of the active solver — they are now determined. Their bucket assignment is final; only the residual WFA-extension runs on Level 3.

### 2.3 Non-convergence

If a read has not collapsed after `max_iter` iterations:
- Status = `Tentative`
- `tentative_targets` = top-K buckets with their final probabilities
- `RejectionReason::DidNotConverge` with full distribution preserved

**No read is forced to a primary.** This is the lossless guarantee at the type level.

### 2.4 Bucket Pyramid (Hierarchical RG-Flow)

| Level | Granularity | Bucket count (human genome) | Purpose |
|---|---|---|---|
| L0 | Chromosomes + repeat families + pseudogene clades | ~1,000 | Coarse routing, contam rejection |
| L1 | 5 MB windows | ~600 | Region localization |
| L2 | 50 kb windows | ~60,000 | Fine localization |
| L3 | Exact position | continuous | WFA2 extension |

Reads start at L0; converged reads drop out, non-converged refine to next level by expanding their bucket candidate set.

### 2.5 Two-Stage Pipeline

```
INPUT FASTQ
    │
    ▼
[Foundation Embedder: Caduceus-Ph distilled, GPU]
    │
    ▼
[STAGE 1: Self-Interference] ◄─── reads inform each other BEFORE projecting onto reference
    │
    │   1.1  FAISS-GPU sparse k-NN over read embeddings
    │   1.2  Leiden community detection on similarity graph
    │   1.3  Self-WaveCollapse: clusters refine each other
    │   1.4  Output: ~1M cluster representatives from ~100M reads
    │
    ▼
[STAGE 2: Reference WaveCollapse] ◄─── now only ~1M reps project to reference
    │
    │   2.1  Bucket pyramid L0 → L1 → L2 → L3
    │   2.2  EM iteration per level with coverage coupling
    │   2.3  Collapse-dropout per level
    │   2.4  WFA2 extension at L3 for residual hard reads
    │
    ▼
[STAGE 3: Member Propagation] ◄─── each cluster member inherits + delta-correction
    │
    ▼
DUAL OUTPUT
    ├─→ BAM-Compat (primary=argmax, alternatives in XA, drop-in for samtools)
    └─→ Probabilistic-Parquet (full lossless: per-read × per-bucket × probability)
```

---

## 3. Tech Stack

- **C++23** core (LLVM/Clang 18+, GCC 13+)
- **CUDA 12.3+** kernels (Compute Capability 8.0+ baseline; Hopper/Blackwell preferred)
- **ONNX Runtime 1.18+** with CUDA Execution Provider for Foundation-Model inference
- **TensorRT 10+** optional acceleration for hot kernels
- **FAISS 1.8+ GPU** for read-vs-read ANN
- **WFA2-lib** for gap-affine extension (C FFI)
- **htslib** for BAM/SAM IO
- **Apache Arrow + Parquet C++** for probabilistic-output
- **RocksDB** for memoization-cache (cross-dataset)
- **libcurl + anthropic-sdk-cpp** for Claude API
- **llama.cpp** optional, only for `--llm self-healing` offline mode
- **CMake 3.28+** build system
- **GoogleTest** unit testing
- **Google Benchmark** + custom harness for performance regression

---

## 4. Project Layout

```
LLmap/
├── CMakeLists.txt
├── README.md                           # tagline + install + quickstart
├── LLmap_SPEC.md                       # this file
├── docs/
│   ├── PHYSICS.md                      # path-integral, RG, symmetry-breaking, T2
│   ├── PHILOSOPHY.md                   # LLM-Map double meaning, photon analogy
│   ├── ARCHITECTURE.md
│   ├── DATA_MODEL.md
│   ├── CLAUDE_AGENT.md                 # how Claude is integrated as tool-user
│   ├── BENCHMARKS.md                   # vs minimap2/winnowmap/strobealign
│   └── MODELS.md                       # foundation-model catalog + training recipes
├── src/
│   ├── core/
│   │   ├── alignment_record.{h,cpp}        # lossless type
│   │   ├── bucket_pyramid.{h,cpp}          # L0..L3 hierarchy
│   │   ├── wave_state.{h,cpp}              # sparse Read×Bucket probability matrix
│   │   ├── em_iterator.{h,cu}              # CUDA: P-update, λ-update, K-smoothing
│   │   ├── coverage_kernel.{h,cu}          # K(b,b') Gaussian smoothing
│   │   ├── collapse_check.{h,cu}           # convergence test + dropout
│   │   ├── refinement.{h,cu}               # coarse→fine bucket expansion
│   │   └── path_integral.{h,cpp}           # multi-path likelihood
│   ├── self_interference/                  # STAGE 1
│   │   ├── faiss_wrapper.{h,cpp}           # FAISS-GPU ANN
│   │   ├── similarity_graph.{h,cpp}        # sparse edge list
│   │   ├── leiden_clustering.{h,cu}        # GPU Leiden
│   │   ├── self_wavecollapse.{h,cu}        # intra-cluster EM
│   │   └── cluster_rep.{h,cpp}             # representative selection
│   ├── reference_collapse/                 # STAGE 2
│   │   ├── stage2_pipeline.{h,cpp}         # orchestrator
│   │   └── member_propagation.{h,cu}       # cluster→member delta-correction
│   ├── ai/
│   │   ├── foundation_embedder.{h,cpp}     # Caduceus-Ph distilled, ONNX
│   │   ├── bucket_embedder.{h,cpp}         # Evo-1.5B distilled, ONNX
│   │   ├── prior_combiner.{h,cpp}          # π_AI(b|r) generation
│   │   ├── nn_triage_head.{h,cpp}          # trivial/repeat/paralog/reject
│   │   ├── nn_paralog_head.{h,cpp}         # PSV→P(paralog)
│   │   ├── nn_confidence_head.{h,cpp}      # convergence predictor
│   │   ├── nn_reject_head.{h,cpp}          # rejection-reason classifier
│   │   └── memoization_store.{h,cpp}       # RocksDB cache
│   ├── claude_agent/
│   │   ├── anthropic_client.{h,cpp}        # API wrapper, async
│   │   ├── session_a_index_builder.{h,cpp} # biology_prior.json generation
│   │   ├── session_b_sample_init.{h,cpp}   # preset + parameter selection
│   │   ├── session_c_diagnostic.{h,cpp}    # stalled-convergence handler + CUDA-gen
│   │   ├── session_d_reporter.{h,cpp}      # post-run report
│   │   ├── cuda_codegen_sandbox.{h,cpp}    # secure compile + load of generated kernels
│   │   └── async_bias_injector.{h,cpp}     # threadsafe bias-term updates
│   ├── classical/
│   │   ├── minimizer_index.{h,cpp}         # for L(r|b)
│   │   ├── chain.{h,cpp}                   # collinear chaining
│   │   ├── wfa2_extend.{h,cpp}             # WFA2 FFI
│   │   └── classical_argmax.{h,cpp}        # CPU-only fallback path
│   ├── psv/                                # paralog-specific
│   │   ├── psv_catalog.{h,cpp}
│   │   ├── inter_paralog_disambig.{h,cpp}
│   │   ├── intra_paralog_disambig.{h,cpp}
│   │   └── informativeness.{h,cpp}
│   ├── singlecell/
│   │   ├── cb_preservation.{h,cpp}         # CB/UB/UMI through pipeline
│   │   └── per_cell_paralog.{h,cpp}        # cell × paralog matrix
│   ├── output/
│   │   ├── bam_compat_writer.{h,cpp}       # htslib-based, drop-in
│   │   ├── probabilistic_parquet.{h,cpp}   # Arrow-based, lossless
│   │   └── coverage_track.{h,cpp}          # probabilistic BedGraph
│   ├── perf/
│   │   ├── mixed_precision.{h,cu}          # FP16/BF16 kernels
│   │   ├── kernel_fusion.{h,cu}            # P-update + λ + K fused
│   │   ├── speculative_prefetch.{h,cpp}    # NN-predicted reference preload
│   │   └── io_uring_reader.{h,cpp}         # async FASTQ IO
│   └── cli/
│       ├── llmap_index.cpp
│       ├── llmap_align.cpp
│       ├── llmap_validate.cpp
│       ├── llmap_bench.cpp
│       └── llmap_cnv.cpp
├── models/
│   ├── README.md                           # versioning + provenance + sha256
│   ├── foundation/
│   │   ├── caduceus_ph_distilled_v1.onnx   # ~50MB
│   │   └── evo_distilled_v1.onnx           # ~200MB
│   ├── heads/
│   │   ├── triage_v1.onnx                  # ~2MB
│   │   ├── paralog_v1.onnx                 # ~20MB
│   │   ├── confidence_v1.onnx              # ~0.5MB
│   │   └── reject_classifier_v1.onnx       # ~4MB
│   └── orchestrator/
│       └── README.md                       # Claude API endpoint config
├── training/
│   ├── README.md
│   ├── distill_caduceus.py                 # Foundation → student
│   ├── distill_evo.py
│   ├── train_triage_head.py
│   ├── train_paralog_head.py
│   ├── train_confidence_head.py
│   ├── train_reject_head.py
│   └── data/
│       └── synthetic_generators/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── synthetic/
│   └── golden/
├── benchmarks/
│   ├── vs_minimap2/
│   ├── vs_winnowmap/
│   └── vs_strobealign/
└── data/
    ├── igh_constants_multi_mrna.fa         # bundled
    └── grch38_chr14_igh.fa                 # bundled, ~2MB
```

---

## 5. Data Model

### 5.1 AlignmentRecord (lossless container)

```cpp
struct AlignmentRecord {
    std::string read_id;
    uint32_t read_len;
    AlignmentStatus status;                   // Mapped | Tentative | Unmapped
    std::optional<AlignmentHit> primary;      // Some when collapsed
    std::vector<AlignmentHit> alternatives;   // all hits ≥ 0.9 × primary score
    std::vector<TentativeTarget> tentative_targets;  // sub-threshold candidates
    std::vector<float> confidence_scores;
    AnchorEvidence anchor_resolution;
    std::optional<ParalogCall> paralog_assignment;
    std::optional<std::string> cell_barcode;
    std::optional<uint8_t> haplotype;
    std::optional<RejectionReason> rejection_reason;

    // WaveCollapse-specific
    uint32_t collapsed_at_iteration;          // 0 if not collapsed
    uint8_t collapsed_at_level;               // L0..L3
    uint32_t cluster_id;                      // Stage-1 cluster membership
    bool is_cluster_representative;
};

enum class AlignmentStatus { Mapped, Tentative, Unmapped };

struct TentativeTarget {
    std::string target_id;
    uint64_t approx_start, approx_end;
    uint32_t n_seeds;
    int32_t partial_chain_score;
    float sequence_identity_estimate;
    float final_probability;                  // P_T(b|r) at termination
};

enum class RejectionReason {
    NoSeeds,
    LowSeedDensity,
    ChainScoreBelowThreshold,
    AmbiguousNoAnchor,
    FailedWfa2Extension,
    DidNotConverge,                           // WaveCollapse-specific
    HostContamination,                        // Triage classifier verdict
    LowComplexity
};
```

### 5.2 Bucket pyramid layout (memory-contiguous)

```cpp
struct BucketPyramid {
    // L0: ~1k buckets — chromosomes, repeat families, pseudogene clades
    std::vector<L0Bucket> l0;

    // L1: ~600 buckets — 5MB windows
    std::vector<L1Bucket> l1;
    std::vector<uint32_t> l1_to_l0;           // parent index

    // L2: ~60k buckets — 50kb windows
    std::vector<L2Bucket> l2;
    std::vector<uint32_t> l2_to_l1;

    // L3: continuous positions, lazy-allocated

    // Embeddings (frozen at index-time)
    arma::Mat<float> l0_embeddings;           // (1k, 256)
    arma::Mat<float> l1_embeddings;           // (600, 256)
    arma::Mat<float> l2_embeddings;           // (60k, 256)

    // Biology prior from Claude Session A
    std::unordered_map<uint32_t, BiologyHint> biology_prior;
};

struct BiologyHint {
    float prior_weight;
    std::optional<std::string> annotation;    // e.g. "IGH-Constants, expect dup"
    std::optional<uint32_t> paralog_partner_bucket;
    float expected_coverage_multiplier;       // 2.0 for known dup-regions
};
```

### 5.3 Wave-state (sparse)

```cpp
// Per read: at most top-K buckets retained (K=10 default, shrinks with level)
struct WaveState {
    // Compressed Sparse Row format for GPU efficiency
    std::vector<uint32_t> read_offsets;       // size N_reads + 1
    std::vector<uint32_t> bucket_indices;     // size = sum(top_K per read)
    std::vector<float> probabilities;         // matched to bucket_indices
    std::vector<uint8_t> current_level;       // L0..L3 per read
    std::vector<bool> collapsed;              // dropout flag

    // GPU mirror
    cudaStream_t stream;
    thrust::device_vector<uint32_t> d_read_offsets;
    thrust::device_vector<uint32_t> d_bucket_indices;
    thrust::device_vector<float> d_probabilities;
};
```

---

## 6. CLI Surface

```bash
# Reference index build (with optional Claude Index-Build agent)
llmap index \
  --fasta GRCh38_plus_alt.fa \
  --gff gencode.v46.gff3 \
  --pangenome HPRC_r2_assemblies/*.fa \
  --bucket-granularity 5000000:50000 \
  --llm sample-aware \
  --foundation-model caduceus-ph-distilled \
  --bucket-model evo-distilled \
  --output GRCh38.llmap.idx

# Alignment — universal across short/long/RNA/DNA, all reference types
llmap align \
  --idx GRCh38.llmap.idx \
  --reads in.{fastq,bam} \
  --preset {sr,map-hifi,map-ont,map-pb,splice,asm5,asm10,asm20,custom} \
  --llm {off,index-only,sample-aware,self-healing,research} \
  --threads 32 \
  --gpu 0,1 \
  --keep-all-alt 0.9 \
  --anchor-resolve \
  --keep-cb-tag \
  --emit-tentative \
  -o out.bam \
  -p out.parquet \
  --unmapped out.unmapped.bam

# Coverage-based CNV (works when per-read disambiguation impossible)
llmap cnv-from-coverage \
  --bam in.bam \
  --parquet in.parquet \
  --regions igh_paralogs.bed \
  --baseline chr14:50M-100M \
  --bootstrap 1000 \
  -o cnv.tsv

# All-vs-all read mapping (Stage 1 exposed as standalone command)
llmap allpair \
  --reads-a iso.bam --reads-b hifi.bam \
  --foundation-model caduceus-ph-distilled \
  --min-cluster-size 5 \
  -o iso_vs_hifi.parquet

# PSV-based paralog assignment
llmap assign-paralog \
  --bam in.bam --parquet in.parquet \
  --psv-catalog psv.json \
  --intra-paralog-mode \
  --phased-bam phaser.bam \
  -o assignment.tsv

# Single-cell paralog matrix
llmap sc-paralog-matrix \
  --parquet in.parquet \
  --assignment assignment.tsv \
  -o cells.h5ad

# Cross-tool benchmark harness
llmap bench \
  --reference GRCh38.fa \
  --synthetic synth_data/ \
  --tools llmap,minimap2,winnowmap,strobealign \
  --metrics recall,precision,paralog-acc,wallclock,peak-ram \
  -o bench_report.html
```

### CLI compatibility shims for minimap2 drop-in

```bash
# minimap2 invocation: minimap2 -ax map-hifi ref.fa reads.fq > out.sam
# LLmap equivalent: llmap-mm ref.fa reads.fq -x map-hifi -o out.sam
# llmap-mm is a thin wrapper that maps minimap2 flags to llmap flags
```

---

## 7. Claude Agent Integration

### 7.1 Sessions overview

| Session | Trigger | Tools Used | Calls (incl. internal) | Cost / sample | Latency-sensitive |
|---|---|---|---|---|---|
| A: Index-Build | Once per reference | Bash, WebFetch, WebSearch, Write, Read | ~50 | $5 amortized | No (ahead of run) |
| B: Sample-Init | Once per sample, before run | Bash, Read | ~20 | $1 | Async to first 5min of run |
| C: Diagnostic | On EM stall trigger | Bash, Read, Write, CUDA-Sandbox | ~100 | $5-15 | Async, additive bias |
| D: Reporter | Once per sample, post-run | Bash, Read, Write | ~30 | $2 | No (after run) |

**Total V1.0 default mode (`sample-aware`)**: ~$3 / sample, all async.

### 7.2 Session A: Index-Build Agent

**Goal**: produce `biology_prior.json` annotated with bucket-level hints.

**Input**: reference FASTA, GFF, pangenome assemblies.

**Claude's tool-use loop**:
1. `Bash: bedtools intersect` reference vs RepeatMasker → identify repeat-dense regions
2. `Bash: samtools faidx` extract candidate paralog regions
3. `WebFetch` UCSC SD-track for current GRCh38/CHM13
4. `WebSearch` recent papers on paralog catalogs for major loci
5. `Write` custom Python script for region-specific bucket boundary logic (e.g. IGH-Constants get 10kb sub-buckets instead of 50kb)
6. `Bash: python custom_bucketize.py` → emits regional bucket overrides
7. `Write biology_prior.json` with annotations + prior weights + expected paralog mappings

**Output schema**:
```json
{
  "version": "1.0",
  "reference_sha256": "...",
  "buckets": {
    "12345": {
      "level": "L2",
      "annotation": "IGH-Constants C-region",
      "prior_weight": 1.2,
      "paralog_partner_bucket": 12346,
      "expected_coverage_multiplier": 2.0,
      "claude_rationale": "Known duplication +650kb in HPRC samples..."
    }
  },
  "regional_overrides": {
    "chr14:105M-107M": {
      "sub_bucket_granularity_kb": 10,
      "max_iter": 25,
      "convergence_threshold": 0.95
    }
  }
}
```

### 7.3 Session B: Sample-Init Agent

**Goal**: choose preset + tune parameters before main run.

**Input**: FASTQ header, first 100k reads (sampled).

**Claude's tool-use loop**:
1. `Bash: seqkit stats` read length distribution
2. `Bash: fastqc` quality profile
3. `Read` sample-sheet metadata if present
4. Reason: "iso-seq from PBMC B-cells → expect 5'-truncation, class-switched IGH dominant"
5. `Write sample_params.json` — preset, threshold, expected coverage profile, Foundation-Model choice

### 7.4 Session C: Diagnostic Agent (with CUDA codegen)

**Trigger**: EM stalls (convergence rate < 10% per iteration for ≥ 3 iterations) OR > 5% of reads non-converged on Level 2.

**Claude's tool-use loop**:
1. `Bash: llmap dump-wave-state --batch 42 -o stalled.json`
2. `Read stalled.json` → analyze bucket distribution
3. `Write investigate.py` for custom diagnostic plot
4. `Bash: python investigate.py` → identify pattern (e.g. bimodal P over two paralogs)
5. `WebFetch` literature on the genomic region
6. **(V1.0 capability)** If a custom CUDA kernel would resolve the issue: `Write custom_kernel.cu` with a specialized convergence-step (e.g. paralog-aware update for sequence-identical exons)
7. `Bash: nvcc --shared -o custom_kernel.so custom_kernel.cu`
8. Pipeline hot-loads `custom_kernel.so` via the **CUDA Sandbox** (see 7.6)
9. Resumed EM iteration uses the custom kernel for the stalled batch

**Self-healing log**: every Session-C invocation generates a learning record in `llmap-experience/` — patterns of which kernel-types resolve which stall-types accumulate over runs.

### 7.5 Session D: Reporter Agent

**Goal**: per-sample diagnostic markdown report + memoization-cache update.

**Claude's tool-use loop**:
1. `Bash: samtools flagstat`, `mosdepth`, etc.
2. `Read` final wave-state summary
3. Reason over patterns: "IGHG4 region shows 4.2× coverage with bimodal cluster — consistent with mosaic duplication"
4. `Write` markdown report
5. `Bash: rocksdb_put` cache generalizable findings

### 7.6 CUDA Sandbox (security)

Claude-generated CUDA code is a security/stability concern. V1.0 mitigation:

1. **Static analysis pass** on generated `.cu` files (no system calls, no file IO, no network) using a vendored AST checker
2. **Compile in container** — nvcc invocation runs inside a rootless `bubblewrap` sandbox with no network and read-only filesystem except `/tmp`
3. **Symbol allow-list** — generated kernels can only call into LLmap's CUDA-helper namespace
4. **Resource limits** — generated kernel must accept a `kernel_budget_ns` parameter; if exceeded, abort and fall back
5. **Audit log** — every generated kernel is signed by Claude session ID + retained for post-hoc review

---

## 8. Foundation Models Catalog

All frozen at runtime. Distillation pipelines live in `training/`.

| Model | Source | Role | Distilled Size | Inference Cost |
|---|---|---|---|---|
| **Caduceus-Ph** (Cornell) | HuggingFace `kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16` | Read embedder | ~50 MB | <10µs/read GPU |
| **Evo 1.5B** (Arc Institute) | HuggingFace `togethercomputer/evo-1.5-8k-base` | Bucket embedder (index-time) | ~200 MB | one-time per ref |
| **Nucleotide-Transformer-500M** (InstaDeep) | HuggingFace `InstaDeepAI/nucleotide_transformer_v2_500m_multi_species` | Cross-species fallback embedder | ~200 MB | <50µs/read GPU |
| **DNABERT-2** (Tsinghua) | HuggingFace `zhihan1996/DNABERT-2-117M` | CPU-only fallback embedder | ~100 MB | <5ms/read CPU |

### Trained-by-us Task Heads

| Head | Backbone | Params | Training Data Source | Trained On |
|---|---|---|---|---|
| Triage | Caduceus-features → 4-layer MLP | 500k | HG002 + HPRC simulated | trivial/repeat/paralog/reject |
| Paralog disambig | Caduceus-features → 4-layer Transformer | 5M | pseudocaller IGHG4↔IGHGP (155 PSVs) + synthetic | inter+intra-paralog probability |
| Confidence | Wave-state features → logistic regression | 10k | EM-convergence logs from HG002 run | will-converge-in-N-iters |
| Reject classifier | Read features → 3-layer MLP | 1M | Manually labeled ~5k unmapped reads | host/repeat/novel/lowqual/chimeric |

**Training budget**: ~2 GPU-weeks total on Hummel-2 (post-distillation). All trainable from public data + your existing HPRC + pseudocaller assets.

---

## 9. Physics Foundation (excerpt from `docs/PHYSICS.md`)

### 9.1 Path-integral semantics

`L(r|b)` is not the max-likelihood alignment score — it is the sum over all alignment trajectories:
```
L(r|b) = Σ_trajectories exp(-S[trajectory] / kT)
```
where `S` is the action functional (gap-penalty + mismatch-penalty integrated along path). This is implemented via Forward algorithm in `path_integral.cpp`, replacing standard Viterbi.

### 9.2 Symmetry breaking

Sequence-identical paralog copies are **degenerate eigenstates** of `L(r|b)`. The coverage-coupling term `K(b,b') · λ(b')` is the **symmetry-breaking field**: collective behavior (mosaic coverage asymmetry) breaks the symmetry that no single read can break alone. This is the mathematical mechanism behind LLmap's paralog-disambiguation USP.

### 9.3 Renormalization group

Bucket pyramid L0→L1→L2→L3 is an RG flow. Each level is a coarse-graining of the next. Fixed points of the flow are stable mappings (collapsed reads). Operators at each level are learned (not hand-designed) during Phase 3 of training.

### 9.4 Decoherence-T2 → platform damping

Sequencing error rate maps to a quantum-decoherence parameter T2:
- PacBio HiFi: Q30+ → long T2 → γ ≈ 0.05
- ONT R10: Q15-20 → short T2 → γ ≈ 0.3
- Illumina: Q35+ → very long T2 → γ ≈ 0.02

Per-platform damping is theory-derived, not heuristic-tuned.

---

## 10. Phase-Ordered Tasks

**Each phase has acceptance criteria. Do not advance until passed. Phase 5 is a kill-switch.**

### Phase 0 — Repo + Foundations + Synthetic Data

- [ ] CMake workspace, C++23, CUDA 12.3, ONNX-Runtime, FAISS-GPU all linking cleanly
- [ ] `AlignmentRecord` type with `Mapped|Tentative|Unmapped` triad enforced at type level
- [ ] `BucketPyramid` data structure with serialization
- [ ] `WaveState` sparse CSR with CPU↔GPU sync primitives
- [ ] Synthetic IGH-locus generator: planted PSVs, controlled mosaic-dup-fraction {0,5,10,30,50,100}%, plus sequence-identical-exon edge case
- [ ] PacBio HiFi simulator (pbsim3 wrapper) + iso-seq simulator (with 5'-truncation profile) + Illumina simulator (ART wrapper)

**Acceptance**: full unit-test suite passes; synthetic data deterministically reproducible by seed; round-trip equality of in-memory and on-disk bucket pyramid.

### Phase 1 — Foundation Model Integration

- [ ] Caduceus-Ph distillation pipeline → 50MB ONNX
- [ ] Evo-1.5B distillation pipeline → 200MB ONNX
- [ ] ONNX-Runtime CUDA-EP wired up; TensorRT optimization for hot models
- [ ] Foundation-Embedder C++ wrapper with GPU batch streaming (target: 10µs/read at batch=10k)
- [ ] Bucket-Embedder runs at index-time, persists embeddings into `.idx`
- [ ] CPU fallback path via DNABERT-2 (slow but functional)

**Acceptance**: 100M reads embed in <2h on 1× H100; embedding-quality validated (read-vs-bucket cosine recovers known mapping ≥ 95% on synthetic).

### Phase 2 — Stage 1: Self-Interference (Read-vs-Read)

- [ ] FAISS-GPU IndexIVFFlat wrapper for read embeddings (10M reads/batch)
- [ ] Sparse k-NN extraction (k=50) → similarity graph in CSR
- [ ] Leiden community detection GPU kernel
- [ ] Self-WaveCollapse: intra-cluster EM, iterative refinement
- [ ] Cluster representative selection (medoid via embedding distance)
- [ ] `llmap allpair` CLI command

**Acceptance**: 100M reads → <1M coherent clusters in <30 min on 1× H100; recall vs naive all-vs-all (on 1M subset) ≥ 99%; clusters split correctly along PSV signatures on synthetic IGHG1/2/3/4 mix.

### Phase 3 — Stage 2: Reference WaveCollapse

- [ ] Bucket pyramid L0..L3 builder
- [ ] CUDA kernel: `em_iteration` — P-update, λ-update, K-smoothing fused
- [ ] CUDA kernel: `collapse_check` — convergence test, dropout
- [ ] CUDA kernel: `refinement` — coarse→fine bucket expansion
- [ ] Member propagation: cluster-rep mapping → member alignment via cheap intra-cluster banded-WFA delta
- [ ] CPU fallback (`classical_argmax`) for no-GPU hosts

**Acceptance**: synthetic IGH locus → reads from canonical and dup-copy correctly separated, dup-fraction recovered within ±2% across mosaic spectrum; convergence in ≤ 15 EM iterations on 95% of reads.

### Phase 4 — Classical Path + WFA2 Extension

- [ ] Minimizer index (matches minimap2 defaults initially)
- [ ] Chain extraction keeping all chains ≥ 0.9× best
- [ ] WFA2-lib FFI for gap-affine extension at Level 3
- [ ] CPU-only fallback runs end-to-end (slow but produces correct output)

**Acceptance**: with `--gpu off --llm off`, LLmap produces output equivalent to minimap2 within ±0.5% recall on HG002 chr14; wallclock ~10× minimap2 on CPU-only path (expected, documented).

### Phase 5 ★ — VALIDATION KILL SWITCH

**This phase is a go/no-go gate. If it fails, project pivots or stops before Phase 6-9.**

- [ ] HG002 HiFi WGS full run: recall ≥ 99.5% of minimap2 on uniquely-mappable regions
- [ ] HG002 IGH-locus targeted run: paralog assignment accuracy ≥ minimap2 baseline + 10 percentage points
- [ ] Synthetic mosaic-dup datasets: dup-fraction recovered within ±2% across {5,10,30,50}% spectrum
- [ ] Stage 1 vs Stage 2 ablation: each stage's contribution measurable and positive
- [ ] No silent read loss: `count(input) == count(AlignmentRecord)` exactly

**KILL CRITERIA — if any of these fail, halt project advancement to Phase 6:**
- Recall < 99% of minimap2 (regardless of paralog gains)
- Paralog accuracy ≤ minimap2 baseline
- Any silent read drop detected
- Convergence rate < 80% on standard HiFi data

If Phase 5 fails: spend ≤ 2 weeks on algorithmic remediation; if still failing, pivot to "minimap2 + lossless wrapper" lower-ambition project.

### Phase 6 — Dual Output (BAM + Parquet)

- [ ] BAM-Compat writer via htslib; primary = argmax, alternatives in `XA` and `om:Z` tags; full samtools/IGV/bcftools compatibility
- [ ] Probabilistic Parquet writer via Arrow C++; schema `(read_id, bucket_id, probability, confidence, level, collapsed_at_iter)`
- [ ] Both files emitted from same `AlignmentRecord` source
- [ ] Round-trip: read Parquet → reconstruct AlignmentRecord → re-emit → byte-identical
- [ ] Lossless invariant test: 1M input reads → 1M Parquet records (no silent drops)

**Acceptance**: `samtools view` and `bcftools call` work unmodified on LLmap BAMs; Parquet readable by `pyarrow`, `polars`, `duckdb`; integration test passes lossless invariant.

### Phase 7 — Claude Agent Integration

- [ ] `anthropic_client.cpp` — async wrapper with token-bucket rate limiting + retry
- [ ] Session A (Index-Builder): produces `biology_prior.json`, integrated into bucket pyramid
- [ ] Session B (Sample-Init): produces `sample_params.json`, applied before run
- [ ] Session C (Diagnostic): triggered on stall; CUDA codegen + sandbox + hot-load functional
- [ ] Session D (Reporter): post-run markdown report
- [ ] CUDA Sandbox: bubblewrap-based; static analyzer rejects unsafe code; symbol allow-list enforced
- [ ] `--llm` flag with all 5 modes (off, index-only, sample-aware, self-healing, research)
- [ ] Async bias-injector: Claude output flows into next EM iteration if ready, never blocks pipeline

**Acceptance**: on HG002 IGH run with `--llm sample-aware`, biology_prior contributes measurable accuracy gain (≥ 3pp paralog); Claude calls do not slow wallclock by more than 5%; CUDA sandbox passes adversarial-input tests (no escape, no resource exhaustion).

### Phase 8 — Performance Optimization

- [ ] Mixed precision: FP16 for embedding similarity, BF16 for EM-update, FP32 for convergence check
- [ ] Kernel fusion: P-update + λ-update + K-smoothing as single CUDA kernel
- [ ] Speculative prefetch: NN-triage predicts target regions → GPU memory preloaded ahead of EM
- [ ] io_uring async FASTQ reader; zero-copy mmap path for indexed FASTA
- [ ] Bench harness vs minimap2/winnowmap/strobealign on 5 datasets (HG002 HiFi, HG002 ONT, HPRC iso-seq, single-cell Kinnex, Illumina WGS)

**Acceptance**: wallclock vs minimap2 on HG002 HiFi: ≤ 0.55× (45%+ speedup); peak RAM ≤ 1.2× minimap2; all metrics published in `docs/BENCHMARKS.md` with reproducible scripts.

### Phase 9 — Single-Cell + Paralog Production

- [ ] CB/UB/UMI preservation across `bam → align → bam` pipeline (no `samtools fastq` round-trip loss)
- [ ] Per-cell paralog matrix → AnnData `.h5ad`
- [ ] PSV catalog builder from HPRC pangenome assemblies
- [ ] Inter-paralog disambiguation (NN paralog-head integrated)
- [ ] Intra-paralog (dup vs canonical) disambiguation
- [ ] Phasing integration (HP-tag propagation from pseudocaller)
- [ ] HG002 IGH end-to-end reproduces published canonical/dup ratios within 5%
- [ ] BioIVT PBMC Kinnex single-cell test: cell × paralog matrix has correct dimensions, B-cell IgG class-switching pattern recovered

**Acceptance**: full IGHG4 use case (the original motivation) demonstrably solved on real HG002 + HPRC data; reproduces your published findings; ready for Belios paper use.

---

## 11. V1.0 Release Acceptance

1. All 10 phases (0-9) pass acceptance criteria
2. **Phase 5 kill-switch passed** — recall + paralog gates met
3. Cross-tool benchmark publicly reproducible; LLmap dominant on SD-region read-recovery, ≥ 45% wallclock speedup average
4. Real HPRC IGH analysis reproduces canonical/dup ratios within 5%
5. Single-cell Kinnex PBMC produces sensible cell × paralog matrix
6. `pip install llmap` works (Python wrapper around C++ binary); brew formula exists
7. `docs/PHYSICS.md`, `docs/PHILOSOPHY.md`, `docs/CLAUDE_AGENT.md` complete
8. Tool paper draft ready (target: Nature Methods or Bioinformatics)
9. Anthropic case study draft ready ("First production bioinformatics tool with Claude in the algorithmic hot-path")

---

## 12. What LLmap explicitly does NOT do

- **Sequence-identical exonic regions, per-read level**: still impossible. LLmap returns `uninformative` per-read but coverage-based CNV via probabilistic depth (Stage 2 by-product) still recovers locus-level CNV.
- **5'-UTR-truncated iso-seq with discriminating PSV in UTR**: cannot recover info that wasn't sequenced. Phase 5d (`readqc`) flags these explicitly.
- **Phasing generation**: LLmap reads HP-tags from upstream phasers (pseudocaller-compatible). It does not generate phasing.
- **Foundation model training**: we do not train Caduceus/Evo/NT from scratch. We distill from pre-trained weights and train only small task heads.
- **Streaming alignment from sequencer**: V1.0 is batch-only. Real-time integration deferred to V1.2+.

---

## 13. Quickstart (target user experience)

```bash
# Install
brew install llmap

# Build index (once per reference, runs Claude Session A if --llm)
llmap index \
  --fasta GRCh38.fa --gff gencode.v46.gff3 \
  --llm sample-aware \
  -o GRCh38.llmap.idx
# ~15 min on H100 + ~5min Claude background; produces GRCh38.llmap.idx (~12 GB)

# Map any data type, any reference
llmap align \
  -i GRCh38.llmap.idx -r sample.fastq.gz \
  -x map-hifi --llm sample-aware \
  -o out.bam -p out.parquet

# Compare to minimap2 (drop-in via wrapper)
llmap-mm GRCh38.fa sample.fastq.gz -x map-hifi > out.sam
```

---

**END OF SPECIFICATION**

*Reads as photons. Genome as crystal. Mapping as decoherence. Every read accounted for.*
