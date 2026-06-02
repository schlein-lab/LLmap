# Worktree changes — SegDup mapper fallback hierarchy

Branch: `feat/nahr-flanks-3318-linter`
Scope: `src/mapper/fallback*`, `tests/mapper/test_fallback*`
Status: implemented, builds clean, all unit tests pass

## TL;DR

a. **Class hierarchy + API** — abstract `FallbackStage` + 5 concrete stages + `FallbackChain` orchestrator. Probes are injectable so tests do not pull in `llmap_classical`.
b. **Test results** — 9/9 GoogleTest pass (one positive test per stage + skip/disable variants + factory + name round-trip).
c. **LOC** — largest .cpp = 212 (fallback_stages.cpp), well under the 400-LOC soft cap.
d. **Catalog-loader interface** — `StageParams` struct in `src/mapper/fallback_types.h` is the agreed-upon hand-off shape. The catalog agent's loader should produce `std::vector<StageParams>` from the JSON `fallback_chain[]` array.

## Files written

### Source — `src/mapper/`

| File | LOC | Purpose |
|---|---:|---|
| `fallback_types.h` | 135 | `FallbackResult`, `StageId`, `StageParams`, `MappingContext`, `CandidatePlacement`, `FallbackOutcome` |
| `fallback_types.cpp` | 18 | `StageIdName()` |
| `fallback_stage.h` | 53 | Abstract `FallbackStage` base + `MakeStage()` factory decl |
| `fallback_stages.h` | 113 | Concrete-stage decls + `StageProbes` injection seam |
| `fallback_stages.cpp` | 212 | 5 concrete `apply()` impls + global probe registry + factory |
| `fallback_chain.h` | 50 | `FallbackChain` container |
| `fallback_chain.cpp` | 54 | `FallbackChain::FromParams` + `run()` |
| `CMakeLists.txt` | 17 | Builds `llmap_mapper` static lib, links `llmap_core` |

### Tests — `tests/mapper/`

| File | LOC | Purpose |
|---|---:|---|
| `test_fallback_stages.cpp` | 319 | 9 GoogleTest cases, all stages + edge cases |
| `CMakeLists.txt` | 13 | Wires `test_fallback_stages` exec |

### Wiring edits (out of scope, but unavoidable for cmake to find the new targets)

- `src/CMakeLists.txt` — one-line `add_subdirectory(mapper)` added after `igh_locus` block.
- `tests/CMakeLists.txt` — one-line `add_subdirectory(mapper)` added at the end.

These are the only two edits outside the `*fallback*` scope. Both are pure additions next to the parallel `classify` module's analogous lines, so the merge surface is minimal.

## Architecture

### Public API

```cpp
namespace llmap::mapper::fallback {

// (a) — what the caller hands in
struct MappingContext {
    std::string_view read_id;
    std::span<const uint8_t> read_seq;
    bool llm_fallback_enabled = false;     // CLI flag --llm-fallback
    std::string_view locus_id;
    // outputs
    std::vector<CandidatePlacement> placements;
    std::vector<std::string> flags;
    std::string last_error;
    StageId last_stage_attempted{};
    StageId resolved_by{};
};

// (b) — per-stage signal
enum class FallbackResult : uint8_t {
    success,                // stop the chain, use ctx.placements
    fail_try_next_stage,    // continue to next stage
    abort_with_warning,     // terminal failure (stage 5)
};

// (c) — what every rung implements
class FallbackStage {
  public:
    virtual FallbackResult apply(MappingContext&) = 0;
    virtual StageId id() const noexcept = 0;
    const char* name() const noexcept;
    const StageParams& params() const noexcept;
};

// (d) — the orchestrator
class FallbackChain {
  public:
    static FallbackChain FromParams(const std::vector<StageParams>&);
    void Push(std::unique_ptr<FallbackStage>);
    FallbackOutcome run(MappingContext&) const;
};

} // namespace
```

### 5 concrete stages

| # | Class | `apply()` behaviour |
|--:|---|---|
| 1 | `RelaxedMismatchStage` | Calls `relaxed_mismatch_probe(read, kmer_size, max_mismatch)`. Catalog default mm=4. |
| 2 | `ChainOnlyStage` | Calls `chain_only_probe(read)`, marks all placements `extension_applied=false`. |
| 3 | `MultiPositionStage` | If `params.report_all_top_k`, calls `multi_position_probe(read, top_k, tie_delta)`. Else fail-through. |
| 4 | `LlmCheckpointStage` (stub) | If `params.require_llm_flag && !ctx.llm_fallback_enabled` → skip. Else call `llm_checkpoint_probe`. Default impl returns empty → fall through. |
| 5 | `NovelHaplotypeFlagStage` | Pushes `warning_tag` into `ctx.flags`, returns `abort_with_warning`. |

### Probe injection seam

Each stage delegates the actual mapping work to a `std::function` in the global `StageProbes` table. Production code installs real probes (talking to `llmap_classical`); tests install mocks. This decouples the chain from the mapper internals — `llmap_mapper` only depends on `llmap_core`, and tests do not need to link any heavy classical/WFA2 path.

```cpp
struct StageProbes {
    std::function<std::vector<CandidatePlacement>(span, uint32_t, uint32_t)>
        relaxed_mismatch_probe;
    std::function<std::vector<CandidatePlacement>(span)>
        chain_only_probe;
    std::function<std::vector<CandidatePlacement>(span, uint32_t, float)>
        multi_position_probe;
    std::function<std::vector<CandidatePlacement>(span, string_view)>
        llm_checkpoint_probe;
};
void SetGlobalStageProbes(StageProbes);
const StageProbes& GlobalStageProbes() noexcept;
StageProbes DefaultStageProbes();   // all no-op
```

## Test results

```
[==========] Running 9 tests from 1 test suite.
[ RUN      ] FallbackChainTest.Stage1RelaxedMismatchSucceeds                       OK
[ RUN      ] FallbackChainTest.Stage2ChainOnlySucceedsAfterStage1Fail              OK
[ RUN      ] FallbackChainTest.Stage3MultiPositionReportsTiedChains                OK
[ RUN      ] FallbackChainTest.Stage3SkippedWhenCatalogDisablesMultiPosition       OK
[ RUN      ] FallbackChainTest.Stage4SkippedWithoutLlmFallbackFlag                 OK
[ RUN      ] FallbackChainTest.Stage4InvokedWithLlmFallbackFlagButStubFails        OK
[ RUN      ] FallbackChainTest.Stage5FlagEmittedWhenAllPriorFail                   OK
[ RUN      ] FallbackChainTest.MakeStageBuildsEveryKnownId                         OK
[ RUN      ] FallbackChainTest.StageNamesAreStable                                 OK
[==========] 9 tests from 1 test suite ran. (0 ms total)
[  PASSED  ] 9 tests.
```

Coverage matrix:

| Spec requirement | Test |
|---|---|
| Stage 1: 3 mismatches in 25mer → success | `Stage1RelaxedMismatchSucceeds` |
| Stage 2: no-extension fallback after stage 1 fail | `Stage2ChainOnlySucceedsAfterStage1Fail` |
| Stage 3: 2 tied chains after stages 1+2 fail | `Stage3MultiPositionReportsTiedChains` |
| Stage 4: skipped without `--llm-fallback` | `Stage4SkippedWithoutLlmFallbackFlag` |
| Stage 4 invoked when flag set (stub) | `Stage4InvokedWithLlmFallbackFlagButStubFails` |
| Stage 5: all prior fail → flag emitted | `Stage5FlagEmittedWhenAllPriorFail` |
| Catalog disables stage 3 → fall-through still works | `Stage3SkippedWhenCatalogDisablesMultiPosition` |
| Factory contract | `MakeStageBuildsEveryKnownId` |
| Name round-trip (catalog JSON ↔ enum) | `StageNamesAreStable` |

## Catalog-loader interface (for the parallel `src/catalog/` agent)

The loader should produce, per locus, a `std::vector<llmap::mapper::fallback::StageParams>` in catalog order. Field mapping from JSON `fallback_chain[]`:

| JSON key | StageParams field | Notes |
|---|---|---|
| `name` | `StageParams::id` | Map via `StageIdName()` (string ↔ enum), see test `StageNamesAreStable` |
| `kmer_size` | `kmer_size` | stage 1 |
| `max_mismatch` | `max_mismatch` | stage 1; default 4 if absent |
| `use_extension` | `use_extension` | stage 2 |
| `report_all_top_k` | `report_all_top_k` | stage 3 |
| `k` | `top_k` | stage 3; default 4 |
| `opt_in_flag` | `opt_in_flag` | stage 4 |
| (implicit) | `require_llm_flag = true` | stage 4 when `opt_in_flag` present |
| `emit_warning` | `emit_warning` | stage 5 |
| `rationale` | `rationale` | optional, free-form |

Then:
```cpp
auto chain = llmap::mapper::fallback::FallbackChain::FromParams(loader_output);
```

`FromParams` silently drops entries whose `id` does not parse to a known stage — the loader is responsible for validation against the JSON schema.

## What's stubbed (intentional)

- **LLM checkpoint probe** (`LlmCheckpointStage`): interface only. Returns `fail_try_next_stage` unless a real `llm_checkpoint_probe` is installed via `SetGlobalStageProbes`. Production wiring goes through `llmap_checkpoint`/`llmap_claude_agent` — out of scope here.
- **Production `DefaultStageProbes()`**: returns no-op probes. The CLI/pipeline that owns the upstream mapper will install real probes that route into `llmap_classical` (minimizer + chain DP + WFA2). That wiring lives in the CLI/integration layer.

## Build verification

```
cd /home/<user>/llmap-local/build
cmake .. -DLLMAP_ENABLE_CUDA=OFF -DLLMAP_ENABLE_FAISS=OFF \
         -DLLMAP_ENABLE_FOUNDATION=OFF -DLLMAP_ENABLE_CLAUDE=OFF
cmake --build . --target test_fallback_stages
./tests/mapper/test_fallback_stages    # 9/9 PASSED
```

## Constraints honored

- No `git commit`, no `git push`, no `git checkout`.
- All new files under `src/mapper/` and `tests/mapper/`.
- Wiring touches (cmake `add_subdirectory`) are 1 line in `src/CMakeLists.txt` and 1 line in `tests/CMakeLists.txt` — additive next to the parallel `classify` module's lines.
- No heavy compute beyond cmake + unit tests (sub-second total runtime).
- Namespace `llmap::mapper::fallback::` consistent across all files.
- C++23 (matches project standard); cxx_std_23 declared on `llmap_mapper`.
- Doxygen-style block comments on every public type and method.
