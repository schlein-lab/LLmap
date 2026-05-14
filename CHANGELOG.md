# Changelog

All notable changes to LLmap are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-05-14

### Added

#### Core Framework (Phase 0)
- `AlignmentRecord` type with full probabilistic state (bucket assignments, confidences, iterations)
- `BucketPyramid` hierarchical spatial index (L0: 1kb, L1: 10kb, L2: 100kb buckets)
- `WaveState` sparse CSR format for efficient probability propagation
- Synthetic IGH-locus generator for validation with configurable mosaic/duplication patterns
- pbsim3/ART simulator wrappers for realistic read generation

#### Foundation Model Integration (Phase 1)
- ONNX Runtime integration with CUDA-EP and TensorRT acceleration
- `FoundationEmbedder` with PIMPL pattern for clean API separation
- `BucketEmbedder` for batch embedding with GPU batching
- Support for CPU, CUDA, and TensorRT execution providers

#### Stage 1: Self-Interference (Phase 2)
- FAISS wrapper with IndexIVFFlat for GPU-accelerated ANN search
- `SimilarityGraph` sparse k-NN extraction in CSR format
- Leiden community detection for read clustering
- `SelfWaveCollapse` intra-cluster EM refinement
- Cluster representative selection via medoid-based algorithm
- `llmap allpair` CLI command for self-interference analysis

#### Stage 2: Reference WaveCollapse (Phase 3)
- `ReferenceIndex` with Builder pattern and efficient serialization
- EM iteration kernel (CPU fallback) with configurable convergence
- Collapse check with entropy-based dropout
- Coarse-to-fine refinement (bucket pyramid traversal)
- Member propagation from cluster representatives
- `Stage2Pipeline` orchestrator wiring all components

#### Classical Path + WFA2 (Phase 4)
- `MinimizerIndex` with configurable k-mer/window parameters
- Chain extraction with colinear DP and 0.9x best filter
- WFA2-lib FFI with Gotoh gap-affine DP fallback
- `ClassicalPipeline` seed-chain-extend end-to-end orchestrator
- Full CPU-only fallback mode

#### Kill-Switch Validation (Phase 5)
- `KillSwitchValidator` framework for lossless/position/origin checks
- End-to-end synthetic validation with configurable presets
- Real reference integration infrastructure (BED/BAM parsing)
- SLURM job management for HPC validation
- `llmap validate-real` CLI command
- `llmap generate-synth` CLI for synthetic data generation
- GPU validation SLURM script (`scripts/slurm_phase55_validation.sh`)

#### Dual Output (Phase 6)
- `BamWriter` for SAM/BAM output with full tag support
- CIGAR generation, statistics, and validation utilities
- WaveCollapse tags: XC (confidence), XL (level), XI (iterations), XP (probability)
- Paralog tags: PD (distance), PC (copy number), PP (paralog probability)
- `ParquetWriter` for probabilistic output (Arrow/Parquet with CSV fallback)
- `ParquetReader` with round-trip validation
- `llmap align` CLI with --bam/--sam/--parquet flags
- Integration tests for full align workflow

#### Claude Agent Integration (Phase 7)
- `AgentTypes` with 4 session types: IndexBuild, SampleInit, Diagnostic, Reporter
- `AnthropicClient` with PIMPL pattern and TokenBucket rate limiter
- `AgentSession` management for async Claude interactions
- `BiologyPrior` JSON serialization for reference annotations
- `CudaSandbox` for secure agent-generated kernel execution
- Static AST analysis for forbidden patterns (syscalls, file I/O, network)
- Bubblewrap containerized compilation (--unshare-net/pid/uts/ipc)
- CUmodule loading with budget enforcement and audit logging
- `PipelineAgent` integration with stall detection (NoProgress, Oscillation, HighUncertainty)
- `--llm` flag family for `llmap align` command

#### Performance Optimization (Phase 8)
- `Profiler` infrastructure with ScopedTimer and ManualTimer
- `Arena` bump allocator for zero-allocation hot paths
- `ScratchBuffer<T>` reusable vector with thread-local `ScratchSpace`
- SIMD optimization (SSE4.2/AVX2) for minimizer extraction
- `EncodeBases`, `PackKmers`, `HashKmersBatch` vectorized operations
- `MmapFastaReader` with lazy index building and zero-copy access
- Memory advice hints (AdviseSequential/Random/WillNeed/DontNeed)
- `ThreadPool` with work-stealing and ParallelFor/ParallelMap APIs
- `BatchProcessor<In,Out>` with progress callbacks
- Cache-friendly SoA layouts: `MinimizerSoA`, `AnchorSoA`, `DPStateSoA`
- `TiledProcessor` for cache-oblivious 2D DP
- Measured 1.56x speedup for sequential field access

#### Single-Cell + Paralog Production (Phase 9)
- `SingleCellTags` struct (CB/CR/CY/UB/UR/UY/RG/BC/QT tags)
- `CellBarcodeWhitelist` with Hamming-distance correction
- Platform presets (10x Genomics, Parse Biosciences, Kinnex)
- `CellParalogMatrix` with aggregation methods (Mean/Max/Sum/Weighted)
- CSV/TSV/dense CSV/H5AD export support
- `llmap sc-paralog-matrix` CLI command
- PSV (Paralog-Specific Variant) module:
  - `PsvCatalog` with BED/VCF I/O and informativeness scoring
  - `PsvAssigner` with Bayesian log-likelihood and posterior normalization
- `llmap align` PSV integration (--psv-catalog/--psv-weight flags)
- Single-cell QC reporting:
  - `CellQcMetrics`, `ParalogQcMetrics`, `GlobalQcSummary`
  - Entropy, dominance, confidence distribution computation
  - JSON/TSV export with quality filtering
- `llmap sc-qc-report` CLI command

#### Production Readiness (Phase 10)
- `Logger` singleton with thread-safe logging
- LogLevel (Trace/Debug/Info/Warn/Error/Fatal/Off) and LogFormat (Text/Json)
- LLMAP_LOG_* macros with source location capture
- Environment configuration (LLMAP_LOG_LEVEL, LLMAP_LOG_FORMAT)
- `Result<T,E>` type for explicit error handling
- `ErrorCode` enum with IO/Parse/Config/Validate/Resource/System/Algo/External categories
- `LLmapError` class with context and source location
- Monadic operations (map, and_then, or_else, inspect) and LLMAP_TRY macro
- `LLmapConfig` struct with TOML parsing via `ConfigParser`
- Config file search paths (./llmap.toml, ~/.config/llmap/config.toml, /etc/llmap/config.toml)
- Environment overrides and CLI flag application
- Version string with --version CLI flag
- `llmap index` CLI command for building minimizer indices
- `llmap check` CLI for V1.0 readiness validation (46 checks, 11 categories)

#### Comparative Benchmark Campaign (Phase 11)
- Benchmark specification vs minimap2/BWA-MEM2/Winnowmap2/STAR/Bowtie2
- Synthetic-truth dataset generator for T1 (WGS) and T2 (paralog stress) tasks
- Tool installation manifest with version verification gate
- Per-tool runners with smoke test framework
- Metrics collector: mapping rate, MAPQ histogram, ground truth evaluation
- Concordance analysis: pairwise BAM comparison
- SLURM orchestrator with dry-run mode, task/tool filtering, replicate management
- Report generator: per-task README, cross-tool tables, matplotlib plots
- T1/T2 benchmark results documented in `docs/BENCHMARKS.md`
- SLURM submission scripts for T3-T6 (real reference tasks)

#### Critical Fixes (Phase A)
- Chain threshold tuning: `min_chain_score` 20→10, `min_score_fraction` 0.9→0.5
- WFA2 extension wiring: `AlignGap()` + `ExtendChain()` for base-accurate CIGAR
- Chain end extension: leftward/rightward extension using WFA2Aligner

#### Performance Improvements (Phase B)
- Parallel alignment: `AlignReadsParallel()` using ThreadPool
- Zero-allocation chaining: thread-local `ChainScratch` buffers
- Index caching: `-i/--index` flag for pre-built index loading

#### Polish + Precision (Phase C)
- Identity filter: default 0.70→0.80 for improved precision
- Filtered alignment tracking (`filtered_by_identity`, `filtered_by_length` stats)
- `-x` presets: `map-hifi` (k=19, w=19, identity=0.90, chain=40), `map-ont`/`map-pb` (k=15, w=10, identity=0.70, chain=20), `sr` (k=21, w=11, identity=0.95, chain=50)
- Probability-based MAPQ calculation using score gap formula
- `--classical-only` mode for pure seed-chain-extend without foundation models

### Architecture

- **1,509 unit tests** with 100% pass rate
- **Modular design**: no source file exceeds 400 LOC
- **Clear dependency graph**: core -> classical -> ai -> reference_collapse -> output -> cli
- **PIMPL pattern** for third-party includes (ONNX Runtime, FAISS)
- **One CMake target per logical unit**
- **Test files mirror source files 1:1**

### CLI Commands

| Command | Description |
|---------|-------------|
| `llmap align` | Full alignment pipeline with BAM/SAM/Parquet output |
| `llmap allpair` | Stage 1 self-interference analysis |
| `llmap index` | Build minimizer index from FASTA |
| `llmap generate-synth` | Generate synthetic validation data |
| `llmap validate-real` | Run kill-switch validation |
| `llmap sc-paralog-matrix` | Build single-cell paralog matrix |
| `llmap sc-qc-report` | Generate single-cell QC report |
| `llmap check` | V1.0 readiness validation |
| `llmap --version` | Show version information |

### Dependencies

- C++23 standard
- CMake 3.20+
- ONNX Runtime 1.17+ (optional, for foundation model)
- FAISS (optional, for GPU-accelerated ANN)
- CUDA 12.0+ (optional, for GPU acceleration)
- Arrow/Parquet (optional, for probabilistic output)

### Notes

- GPU validation on Hummel-2 HPC pending manual SLURM submission
- CPU validation passes with 100% recall on synthetic IGH locus data
- Build tested with GCC 13.2 and Clang 17

## [Unreleased]

### Planned
- Streaming alignment from sequencer
- Phase integration (consume HP-tags)
- Pangenome index format
- Multi-sample joint calling mode

[1.0.0]: https://github.com/schlein-lab/LLmap/releases/tag/v1.0.0
[Unreleased]: https://github.com/schlein-lab/LLmap/compare/v1.0.0...HEAD
