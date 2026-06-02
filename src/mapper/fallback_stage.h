// LLmap — SegDup mapper fallback hierarchy: abstract stage interface.
//
// Each stage is a pure pass over a MappingContext. Stages MUST be
// stateless across reads — any state lives in the context. This makes
// the chain trivially thread-safe (one chain per thread) and lets us
// swap in mock stages for tests without inheritance gymnastics.

#pragma once

#include "mapper/fallback_types.h"

#include <memory>
#include <string>

namespace llmap::mapper::fallback {

/// Abstract base for a single rung of the SegDup fallback ladder.
///
/// Lifetime model: stages are owned by a FallbackChain (unique_ptr).
/// Construct once from catalog params, call apply() many times.
class FallbackStage {
public:
    explicit FallbackStage(StageParams params) : params_(std::move(params)) {}
    virtual ~FallbackStage() = default;

    FallbackStage(const FallbackStage&) = delete;
    FallbackStage& operator=(const FallbackStage&) = delete;

    /// Run this stage against `ctx`. Implementations:
    ///   - may write to ctx.placements / ctx.flags / ctx.last_error
    ///   - MUST set ctx.last_stage_attempted to id()
    ///   - return one of {success, fail_try_next_stage, abort_with_warning}
    virtual FallbackResult apply(MappingContext& ctx) = 0;

    /// Stable identifier. Used by the chain to record which stage
    /// resolved the read and for structured logging.
    virtual StageId id() const noexcept = 0;

    /// Human-readable name. Mirrors catalog JSON `name`.
    const char* name() const noexcept { return StageIdName(id()); }

    const StageParams& params() const noexcept { return params_; }

protected:
    StageParams params_;
};

/// Convenience factory: instantiate the correct concrete stage given
/// StageParams.id. Returns nullptr for unknown ids (callers should
/// treat that as a catalog-load error).
std::unique_ptr<FallbackStage> MakeStage(StageParams params);

}  // namespace llmap::mapper::fallback
