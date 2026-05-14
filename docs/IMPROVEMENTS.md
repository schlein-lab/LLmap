# LLmap Improvement Targets

This document identifies performance gaps revealed by Phase 11 benchmarks (T1/T2 synthetic tests) and prioritizes improvement targets for future development iterations.

**Generated**: 2026-05-14
**Baseline Version**: V1.0

---

## Executive Summary

T1/T2 benchmarks reveal three critical performance gaps vs. minimap2:

| Gap | LLmap | minimap2 | Impact |
|-----|-------|----------|--------|
| **Mapping Rate** | 46% | 92% | 50% of reads unmapped |
| **Recall** | 38% | 100% | Missing correct positions |
| **Wallclock** | 32s | 9s | 3.5× slower |

Root cause analysis indicates these stem from overly stringent filtering thresholds and an incomplete extension implementation, not fundamental algorithmic limitations.

---

## Issue 1: Low Mapping Rate (Critical)

### Symptoms
- T1: 46.2% mapping rate vs. 91.8% (minimap2)
- T2: 45.0% mapping rate vs. 91.9% (minimap2)
- ~50% of reads produce no alignment

### Root Cause Analysis

**Primary cause: Chain score threshold too aggressive**

Current `ChainConfig` defaults (`src/classical/chain.h:27`):
```cpp
int32_t min_chain_score = 20;       // Too high for short chains
uint32_t min_chain_anchors = 3;     // Reasonable
float min_score_fraction = 0.9f;    // Extremely aggressive secondary filter
```

Reads with valid but modest chain scores (10-19) are discarded entirely. The 90% score fraction filter aggressively removes chains that might be valid alignments.

**Secondary cause: Minimizer suppression too aggressive**

Current `MinimizerConfig` defaults (`src/classical/minimizer_index.h:30`):
```cpp
size_t max_occ = 500;  // May suppress useful repetitive minimizers
```

While 500 is reasonable, repetitive regions may need higher thresholds for mapping.

### Recommended Fixes

| Priority | Fix | File | Change |
|----------|-----|------|--------|
| P0 | Reduce `min_chain_score` to 10 | `chain.h` | `int32_t min_chain_score = 10;` |
| P0 | Reduce `min_score_fraction` to 0.5 | `chain.h` | `float min_score_fraction = 0.5f;` |
| P1 | Add `--min-chain-score` CLI flag | `cmd_align.cpp` | User-tunable threshold |
| P2 | Increase `max_occ` to 1000 for repeat-rich refs | `minimizer_index.h` | Configurable |

### Expected Impact
- Mapping rate: 46% → 80-90%
- May slightly increase false positives (address via P1 MAPQ tuning)

---

## Issue 2: Low Recall (Critical)

### Symptoms
- T1: 38.7% recall vs. 100% (minimap2)
- T2: 36.0% recall vs. 100% (minimap2)
- Mapped reads often placed incorrectly (within 10bp tolerance)

### Root Cause Analysis

**Primary cause: Extension generates approximate CIGAR, not base-accurate**

In `classical_pipeline.cpp:207-258`, the `ExtendChain` function builds CIGAR by interpolating between anchor positions rather than performing actual sequence alignment:

```cpp
// Build a simple CIGAR from chain anchors
// For a proper implementation, we would do full WFA2 extension
// between anchors. For now, we approximate with M operations
// connecting the anchors.
```

This means:
1. Alignment coordinates are approximate (anchor-derived, not base-aligned)
2. CIGAR strings are synthetic, not computed from actual alignment
3. Position errors accumulate between anchors

**Secondary cause: WFA2 aligner not invoked**

The `Wfa2Aligner` module exists (`src/classical/wfa2_aligner.cpp`) but is not called during extension. The pipeline relies on anchor interpolation.

### Recommended Fixes

| Priority | Fix | File | Change |
|----------|-----|------|--------|
| P0 | Invoke WFA2 between adjacent anchors | `classical_pipeline.cpp` | Call `aligner_.Extend()` |
| P0 | Compute base-accurate CIGAR | `classical_pipeline.cpp` | Replace interpolation logic |
| P1 | Add left/right extension from chain ends | `classical_pipeline.cpp` | Call `aligner_.ExtendLeft/Right()` |
| P2 | Implement anchor-skipping for gap bridging | `classical_pipeline.cpp` | Handle > max_gap cases |

### Implementation Sketch

Replace `ExtendChain()` with:
```cpp
// Between each pair of adjacent anchors:
//   1. Extract ref substring: ref_seq[prev_end : curr_start]
//   2. Extract query substring: query_seq[prev_qend : curr_qstart]
//   3. Call aligner_.Align(query_sub, ref_sub)
//   4. Append resulting CIGAR
// Extend from first anchor leftward, last anchor rightward
```

### Expected Impact
- Recall: 38% → 95%+ (limited by seeding)
- Precision maintained or improved
- F1: 48% → 97%+

---

## Issue 3: Slow Wallclock Time (High)

### Symptoms
- T1: 31.9s vs. 9.1s (3.5× slower)
- T2: 16.2s vs. 5.6s (2.9× slower)

### Root Cause Analysis

**Primary causes:**

1. **No index caching**: Index rebuilt per alignment batch in CLI
2. **Sequential read processing**: Single-threaded in `AlignReads()` (`classical_pipeline.cpp:288-325`)
3. **Heap allocations in hot path**: Despite arena work, chaining still allocates

**Profile breakdown (estimated):**
- Index query: 15%
- Chaining: 35%
- Extension (interpolation): 10%
- I/O: 40%

### Recommended Fixes

| Priority | Fix | File | Change |
|----------|-----|------|--------|
| P1 | Parallelize `AlignReads` with ThreadPool | `classical_pipeline.cpp` | Use `ParallelFor` |
| P1 | Cache built index between batches | `cmd_align.cpp` | Load once, align many |
| P2 | Use `ExtractChainsFromAnchorsWithScratch` | `classical_pipeline.cpp` | Zero-alloc chaining |
| P2 | Batch I/O with larger read chunks | `cmd_align.cpp` | Buffer reads |
| P3 | Vectorize minimizer extraction (SIMD) | `minimizer_index_extract.cpp` | Use AVX2 from `simd.h` |

### Implementation Sketch

```cpp
// In AlignReads():
auto& pool = core::GetGlobalThreadPool();
results.resize(query_names.size());
pool.ParallelFor(0, query_names.size(), [&](size_t i) {
    results[i] = AlignRead(query_names[i], query_seqs[i]);
}, /*chunk_size=*/32);
```

### Expected Impact
- With parallelization: 32s → 8-12s on 4 cores
- With all optimizations: 32s → 5-8s (competitive with minimap2)

---

## Issue 4: Precision Gap (Medium)

### Symptoms
- T1: 62.3% precision (mapped reads often incorrect)
- T2: 58.5% precision

### Root Cause Analysis

Reads with weak chains are still reported. Without proper extension, even "mapped" reads may be positioned incorrectly.

### Recommended Fixes

| Priority | Fix | File | Change |
|----------|-----|------|--------|
| P1 | Add identity filter after extension | `classical_pipeline.cpp` | Filter `identity < 0.8` |
| P2 | Implement proper MAPQ calculation | `classical_pipeline.cpp` | Based on alignment score |
| P2 | Add `-x` presets (map-hifi, map-ont) | `cmd_align.cpp` | Tuned thresholds |

---

## Issue 5: Memory Overhead (Low)

### Symptoms
- T1: 305 MB vs. 252 MB (21% overhead)

### Root Cause Analysis

Probabilistic framework structures allocated even for classical pipeline path.

### Recommended Fixes

| Priority | Fix | File | Change |
|----------|-----|------|--------|
| P3 | Lazy allocation of Wave structures | `wave_state.cpp` | Allocate on use |
| P3 | Add `--classical-only` mode | `cmd_align.cpp` | Skip probabilistic init |

---

## Prioritized Implementation Plan

### Phase A: Critical Fixes (Next 2 Iterations)

1. **A.1**: Adjust chain thresholds (`min_chain_score=10`, `min_score_fraction=0.5`)
2. **A.2**: Wire WFA2 extension into `ExtendChain()`
3. **A.3**: Add left/right extension for chain ends

**Goal**: Mapping rate ≥85%, Recall ≥90%, F1 ≥90%

### Phase B: Performance (Following 2 Iterations)

4. **B.1**: Parallelize `AlignReads()` with ThreadPool
5. **B.2**: Use zero-allocation chaining (`ChainScratch`)
6. **B.3**: Index caching in CLI

**Goal**: Wallclock ≤12s on T1 (within 30% of minimap2)

### Phase C: Polish (Following 2 Iterations)

7. **C.1**: Add CLI flags for tunable thresholds
8. **C.2**: Implement `-x` presets
9. **C.3**: Proper MAPQ calculation
10. **C.4**: Add `--classical-only` mode

**Goal**: Production-quality alignment with competitive metrics

---

## Validation Criteria

After implementing Phase A fixes, re-run T1/T2 benchmarks and verify:

| Metric | Current | Target |
|--------|---------|--------|
| Mapping Rate | 46% | ≥85% |
| Recall | 38% | ≥90% |
| Precision | 62% | ≥95% |
| F1 | 48% | ≥90% |
| Wallclock | 32s | ≤40s (no regression) |

After Phase B:
- Wallclock ≤12s (T1)
- Memory ≤280 MB (T1)

---

## Technical Notes

### Why minimap2 achieves 100% F1 on synthetic data

Synthetic truth datasets (T1/T2) simulate reads with exact known origin. Deterministic mappers like minimap2 can perfectly recover these origins because:
1. No sequencing errors in simulation
2. No ambiguous regions in synthetic references
3. Perfect chain recovery with tuned parameters

Real-world data (T3-T6) will show more realistic performance gaps where LLmap's probabilistic framework may provide advantages.

### WaveCollapse framework not engaged

T1/T2 benchmarks use only the classical seed-chain-extend pipeline. The foundation-model embeddings and WaveCollapse probabilistic refinement are designed for paralog disambiguation in ambiguous regions. These features will be evaluated in T6 (IGH locus) benchmarks.

---

## References

- `src/classical/chain.h` - Chain configuration
- `src/classical/classical_pipeline.cpp` - Pipeline orchestration
- `src/classical/wfa2_aligner.cpp` - WFA2 alignment (unused)
- `benchmarks/reports/T1/README.md` - T1 results
- `benchmarks/reports/T2/README.md` - T2 results
