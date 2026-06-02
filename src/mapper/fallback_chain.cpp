// LLmap — SegDup mapper fallback hierarchy: chain orchestrator.

#include "mapper/fallback_chain.h"

#include <utility>

namespace llmap::mapper::fallback {

FallbackChain FallbackChain::FromParams(
    const std::vector<StageParams>& params) {
    FallbackChain chain;
    chain.stages_.reserve(params.size());
    for (const auto& p : params) {
        auto stage = MakeStage(p);
        if (stage) {
            chain.stages_.push_back(std::move(stage));
        }
    }
    return chain;
}

void FallbackChain::Push(std::unique_ptr<FallbackStage> stage) {
    if (stage) {
        stages_.push_back(std::move(stage));
    }
}

FallbackOutcome FallbackChain::run(MappingContext& ctx) const {
    FallbackOutcome out;
    for (const auto& stage : stages_) {
        ++out.stages_run;
        const auto r = stage->apply(ctx);
        if (r == FallbackResult::success) {
            out.success = true;
            out.resolved_by = stage->id();
            out.flags = ctx.flags;
            return out;
        }
        if (r == FallbackResult::abort_with_warning) {
            // Terminal stage emitted a warning. Stop the chain but
            // surface the flag set for downstream consumers.
            out.success = false;
            out.flags = ctx.flags;
            out.last_error = ctx.last_error;
            return out;
        }
        // fail_try_next_stage: continue.
    }
    out.flags = ctx.flags;
    out.last_error = ctx.last_error;
    return out;
}

}  // namespace llmap::mapper::fallback
