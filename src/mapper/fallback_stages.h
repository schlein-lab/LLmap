// LLmap — SegDup mapper fallback hierarchy: concrete stages.
//
// Five rungs, in catalog order:
//   1. RelaxedMismatchStage   — retry with max_mismatch 2->4
//   2. ChainOnlyStage         — skip extension; report best chain
//   3. MultiPositionStage     — report top-k tied chains
//   4. LlmCheckpointStage     — opt-in LLM consult (stub here)
//   5. NovelHaplotypeFlagStage — emit warning, surface to curation
//
// Each stage exposes an injectable "probe" function. The probe is the
// part that actually talks to the upstream mapper (k-mer index +
// chain + WFA2 extension). Default probes are linked against the real
// pipeline; tests inject mocks via the StageProbes adapter so the unit
// tests do not pull in llmap_classical.

#pragma once

#include "mapper/fallback_stage.h"

#include <functional>
#include <memory>
#include <vector>

namespace llmap::mapper::fallback {

/// A "probe" abstracts the upstream mapping call each stage performs.
/// Production code wires these to llmap_classical (minimizer + chain +
/// WFA2). Tests inject mocks.
///
/// Return: list of candidate placements found (possibly empty).
struct StageProbes {
    /// Stage 1: re-run mapping with relaxed mismatch budget.
    /// args: read_seq, kmer_size, max_mismatch
    std::function<std::vector<CandidatePlacement>(
        std::span<const uint8_t>, uint32_t, uint32_t)>
        relaxed_mismatch_probe;

    /// Stage 2: produce best chain without extension.
    /// args: read_seq
    std::function<std::vector<CandidatePlacement>(
        std::span<const uint8_t>)>
        chain_only_probe;

    /// Stage 3: report top-k tied chains.
    /// args: read_seq, k, tie_delta_score
    std::function<std::vector<CandidatePlacement>(
        std::span<const uint8_t>, uint32_t, float)>
        multi_position_probe;

    /// Stage 4: LLM checkpoint. Returns empty in stub.
    /// args: read_seq, locus_id
    std::function<std::vector<CandidatePlacement>(
        std::span<const uint8_t>, std::string_view)>
        llm_checkpoint_probe;
};

/// Install the production probes (calls into llmap_classical / llmap_checkpoint).
/// Returns a StageProbes ready to hand to MakeStage(...).
/// NOTE: real wiring lives outside this scope; here we only provide
/// the seam so the chain compiles independently.
StageProbes DefaultStageProbes();

/// Override the probes for the current process. Used by tests; thread-
/// unsafe by design (set once, fork tests as needed).
void SetGlobalStageProbes(StageProbes probes);

/// Read back whatever was last installed (or DefaultStageProbes() if
/// nothing was set).
const StageProbes& GlobalStageProbes() noexcept;

// ---------------------------------------------------------------------------
// Concrete stages
// ---------------------------------------------------------------------------

class RelaxedMismatchStage final : public FallbackStage {
public:
    using FallbackStage::FallbackStage;
    FallbackResult apply(MappingContext& ctx) override;
    StageId id() const noexcept override { return StageId::relaxed_mismatch; }
};

class ChainOnlyStage final : public FallbackStage {
public:
    using FallbackStage::FallbackStage;
    FallbackResult apply(MappingContext& ctx) override;
    StageId id() const noexcept override { return StageId::chain_only; }
};

class MultiPositionStage final : public FallbackStage {
public:
    using FallbackStage::FallbackStage;
    FallbackResult apply(MappingContext& ctx) override;
    StageId id() const noexcept override { return StageId::multi_position; }
};

/// Stub. Always returns fail_try_next_stage unless the --llm-fallback
/// flag is set in the MappingContext AND a non-null
/// llm_checkpoint_probe is installed.
class LlmCheckpointStage final : public FallbackStage {
public:
    using FallbackStage::FallbackStage;
    FallbackResult apply(MappingContext& ctx) override;
    StageId id() const noexcept override { return StageId::llm_checkpoint; }
};

class NovelHaplotypeFlagStage final : public FallbackStage {
public:
    using FallbackStage::FallbackStage;
    FallbackResult apply(MappingContext& ctx) override;
    StageId id() const noexcept override { return StageId::novel_haplotype_flag; }
};

}  // namespace llmap::mapper::fallback
