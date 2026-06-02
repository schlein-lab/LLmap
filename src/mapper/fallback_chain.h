// LLmap — SegDup mapper fallback hierarchy: chain container.
//
// FallbackChain is the orchestrator. It owns an ordered list of stages
// (built once from a catalog locus) and applies them to each read
// until one returns success or the chain is exhausted.

#pragma once

#include "mapper/fallback_stage.h"
#include "mapper/fallback_types.h"

#include <memory>
#include <vector>

namespace llmap::mapper::fallback {

/// Ordered chain of fallback stages. Construct from a vector of
/// StageParams (one per stage, in execution order); call run() per
/// read.
class FallbackChain {
public:
    FallbackChain() = default;

    /// Build from a flat vector of StageParams (catalog order).
    /// Unknown stage ids are silently dropped — the caller is
    /// responsible for validating before construction.
    static FallbackChain FromParams(const std::vector<StageParams>& params);

    /// Append a pre-built stage. Mostly for tests.
    void Push(std::unique_ptr<FallbackStage> stage);

    /// Number of stages in the chain.
    std::size_t size() const noexcept { return stages_.size(); }
    bool empty() const noexcept { return stages_.empty(); }

    /// Read-only access — used by diagnostics and tests.
    const FallbackStage& at(std::size_t i) const { return *stages_.at(i); }

    /// Run the chain against `ctx`. Returns an aggregate outcome.
    /// Side-effects on `ctx`:
    ///   - ctx.placements:   filled by whichever stage succeeded
    ///   - ctx.flags:        accumulated across all stages
    ///   - ctx.resolved_by:  set when a stage returns success
    FallbackOutcome run(MappingContext& ctx) const;

private:
    std::vector<std::unique_ptr<FallbackStage>> stages_;
};

}  // namespace llmap::mapper::fallback
