// LLmap — SegDup mapper fallback hierarchy: stage implementations.
//
// Each apply() method:
//   1. Marks ctx.last_stage_attempted with its own id().
//   2. Calls the matching probe in the global StageProbes table.
//   3. Translates the probe's empty/non-empty result into a
//      FallbackResult.
//
// The probes are injectable so unit tests can drive deterministic
// scenarios without pulling in llmap_classical.

#include "mapper/fallback_stages.h"

#include <string>
#include <utility>

namespace llmap::mapper::fallback {

namespace {

/// Process-wide probe table. Set via SetGlobalStageProbes().
StageProbes& MutableProbes() {
    static StageProbes probes = DefaultStageProbes();
    return probes;
}

}  // namespace

StageProbes DefaultStageProbes() {
    // Production wiring lives in another translation unit (cli or a
    // dedicated mapper_pipeline.cpp). Here we provide no-op defaults
    // so the chain still links and behaves safely if the host forgot
    // to install probes — every stage will simply fail_try_next_stage.
    StageProbes p;
    p.relaxed_mismatch_probe = [](std::span<const uint8_t>, uint32_t, uint32_t) {
        return std::vector<CandidatePlacement>{};
    };
    p.chain_only_probe = [](std::span<const uint8_t>) {
        return std::vector<CandidatePlacement>{};
    };
    p.multi_position_probe = [](std::span<const uint8_t>, uint32_t, float) {
        return std::vector<CandidatePlacement>{};
    };
    p.llm_checkpoint_probe = [](std::span<const uint8_t>, std::string_view) {
        return std::vector<CandidatePlacement>{};
    };
    return p;
}

void SetGlobalStageProbes(StageProbes probes) {
    MutableProbes() = std::move(probes);
}

const StageProbes& GlobalStageProbes() noexcept {
    return MutableProbes();
}

// ---------------------------------------------------------------------------
// Stage 1: relaxed_mismatch
// ---------------------------------------------------------------------------

FallbackResult RelaxedMismatchStage::apply(MappingContext& ctx) {
    ctx.last_stage_attempted = id();
    const auto& probes = GlobalStageProbes();
    if (!probes.relaxed_mismatch_probe) {
        ctx.last_error = "relaxed_mismatch_probe not installed";
        return FallbackResult::fail_try_next_stage;
    }

    auto hits = probes.relaxed_mismatch_probe(
        ctx.read_seq, params_.kmer_size, params_.max_mismatch);

    if (hits.empty()) {
        ctx.last_error = "no hits at max_mismatch=" +
                         std::to_string(params_.max_mismatch);
        return FallbackResult::fail_try_next_stage;
    }

    ctx.clear_placements();
    ctx.placements = std::move(hits);
    ctx.resolved_by = id();
    return FallbackResult::success;
}

// ---------------------------------------------------------------------------
// Stage 2: chain_only
// ---------------------------------------------------------------------------

FallbackResult ChainOnlyStage::apply(MappingContext& ctx) {
    ctx.last_stage_attempted = id();
    const auto& probes = GlobalStageProbes();
    if (!probes.chain_only_probe) {
        ctx.last_error = "chain_only_probe not installed";
        return FallbackResult::fail_try_next_stage;
    }

    auto hits = probes.chain_only_probe(ctx.read_seq);
    if (hits.empty()) {
        ctx.last_error = "chain-only: no chain";
        return FallbackResult::fail_try_next_stage;
    }

    for (auto& h : hits) {
        h.extension_applied = false;
    }
    ctx.clear_placements();
    ctx.placements = std::move(hits);
    ctx.resolved_by = id();
    return FallbackResult::success;
}

// ---------------------------------------------------------------------------
// Stage 3: multi_position
// ---------------------------------------------------------------------------

FallbackResult MultiPositionStage::apply(MappingContext& ctx) {
    ctx.last_stage_attempted = id();
    if (!params_.report_all_top_k) {
        // Catalog says don't allow multi-position for this locus.
        ctx.last_error = "multi_position disabled by catalog";
        return FallbackResult::fail_try_next_stage;
    }

    const auto& probes = GlobalStageProbes();
    if (!probes.multi_position_probe) {
        ctx.last_error = "multi_position_probe not installed";
        return FallbackResult::fail_try_next_stage;
    }

    auto hits = probes.multi_position_probe(
        ctx.read_seq, params_.top_k, params_.tie_delta_score);

    if (hits.empty()) {
        ctx.last_error = "multi_position: no tied chains";
        return FallbackResult::fail_try_next_stage;
    }

    ctx.clear_placements();
    ctx.placements = std::move(hits);
    ctx.resolved_by = id();
    return FallbackResult::success;
}

// ---------------------------------------------------------------------------
// Stage 4: llm_checkpoint (stub)
// ---------------------------------------------------------------------------

FallbackResult LlmCheckpointStage::apply(MappingContext& ctx) {
    ctx.last_stage_attempted = id();

    // Honor the catalog's opt-in flag. If the CLI did not pass
    // --llm-fallback, we skip silently.
    if (params_.require_llm_flag && !ctx.llm_fallback_enabled) {
        ctx.last_error = "llm_checkpoint skipped (no " +
                         params_.opt_in_flag + ")";
        return FallbackResult::fail_try_next_stage;
    }

    const auto& probes = GlobalStageProbes();
    if (!probes.llm_checkpoint_probe) {
        ctx.last_error = "llm_checkpoint_probe not installed (stub mode)";
        return FallbackResult::fail_try_next_stage;
    }

    auto hits = probes.llm_checkpoint_probe(ctx.read_seq, ctx.locus_id);
    if (hits.empty()) {
        ctx.last_error = "llm_checkpoint: inconclusive";
        return FallbackResult::fail_try_next_stage;
    }
    ctx.clear_placements();
    ctx.placements = std::move(hits);
    ctx.resolved_by = id();
    return FallbackResult::success;
}

// ---------------------------------------------------------------------------
// Stage 5: novel_haplotype_flag
// ---------------------------------------------------------------------------

FallbackResult NovelHaplotypeFlagStage::apply(MappingContext& ctx) {
    ctx.last_stage_attempted = id();
    if (params_.emit_warning) {
        ctx.flags.push_back(params_.warning_tag);
    }
    // Terminal stage: chain stops regardless. We report
    // abort_with_warning so callers can distinguish "no resolution"
    // from "no stage applicable".
    ctx.last_error = "novel_haplotype: no canonical/dup PSV pattern";
    return FallbackResult::abort_with_warning;
}

// ---------------------------------------------------------------------------
// Factory
// ---------------------------------------------------------------------------

std::unique_ptr<FallbackStage> MakeStage(StageParams params) {
    switch (params.id) {
        case StageId::relaxed_mismatch:
            return std::make_unique<RelaxedMismatchStage>(std::move(params));
        case StageId::chain_only:
            return std::make_unique<ChainOnlyStage>(std::move(params));
        case StageId::multi_position:
            return std::make_unique<MultiPositionStage>(std::move(params));
        case StageId::llm_checkpoint:
            return std::make_unique<LlmCheckpointStage>(std::move(params));
        case StageId::novel_haplotype_flag:
            return std::make_unique<NovelHaplotypeFlagStage>(std::move(params));
    }
    return nullptr;
}

}  // namespace llmap::mapper::fallback
