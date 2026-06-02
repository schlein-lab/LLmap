// LLmap — SegDup mapper fallback hierarchy: shared types.
//
// This header carries the small POD types shared by every fallback stage:
//   - FallbackResult: per-stage status (success / try-next / abort)
//   - StageParams:    catalog-driven knobs (one struct per stage)
//   - MappingContext: per-read mutable state passed through the chain
//   - FallbackOutcome: aggregate result the chain returns to the caller
//
// Catalog loader interface (parallel agent in src/catalog/):
//   The catalog loader is expected to produce a `MappingStrategy` object
//   that contains a `std::vector<StageParams>` with one entry per stage
//   defined in the JSON `fallback_chain` array. We define the canonical
//   StageParams here so both sides agree on field names and types.
//
//   Concretely:
//     llmap::catalog::CuratedLocus::mapping_strategy.fallback ->
//         std::vector<llmap::mapper::fallback::StageParams>
//   The catalog side fills in the params; the mapper side consumes them.

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::mapper::fallback {

/// Per-stage outcome semantics.
///
/// success                — read mapped, stop the chain
/// fail_try_next_stage    — stage produced no usable mapping; continue
/// abort_with_warning     — terminal failure (only the novel-haplotype
///                          flag stage uses this; everything else uses
///                          fail_try_next_stage so the chain can run on)
enum class FallbackResult : uint8_t {
    success = 0,
    fail_try_next_stage = 1,
    abort_with_warning = 2,
};

/// Stable identity tag for each stage. Mirrors the JSON `name` field
/// in catalog/curated/*.json `fallback_chain[].name`.
enum class StageId : uint8_t {
    relaxed_mismatch     = 1,
    chain_only           = 2,
    multi_position       = 3,
    llm_checkpoint       = 4,
    novel_haplotype_flag = 5,
};

const char* StageIdName(StageId id) noexcept;

/// Catalog-driven knobs. Every stage reads from this struct; unused
/// fields are simply ignored. This keeps the catalog->mapper interface
/// flat and trivially serializable.
struct StageParams {
    StageId id{};

    // Stage 1 (relaxed_mismatch)
    uint32_t kmer_size = 25;
    uint32_t max_mismatch = 4;

    // Stage 2 (chain_only)
    bool use_extension = true;

    // Stage 3 (multi_position)
    bool report_all_top_k = false;
    uint32_t top_k = 4;
    float tie_delta_score = 0.05f;   // fractional Δscore for tie set

    // Stage 4 (llm_checkpoint)
    bool require_llm_flag = true;    // skip unless --llm-fallback is set
    std::string opt_in_flag = "--llm-fallback";

    // Stage 5 (novel_haplotype_flag)
    bool emit_warning = true;
    std::string warning_tag = "novel_haplotype";

    // Optional free-form note carried over from JSON (rationale string).
    std::string rationale;
};

/// A single candidate placement reported by the chain. The chain is
/// agnostic to upstream chain/anchor types — it only needs the minimal
/// fields downstream callers (BAM writer, classifier) consume.
struct CandidatePlacement {
    uint32_t ref_id = 0;
    uint32_t ref_start = 0;
    uint32_t ref_end = 0;
    int32_t  score = 0;
    bool     forward = true;
    uint8_t  mapq = 0;
    std::string copy_tag;          // e.g. "canonical", "dup", "dupVar"
    bool extension_applied = true; // false when stage 2 short-circuited
};

/// Per-read state threaded through the chain. Implementations of
/// FallbackStage::apply() mutate `placements`, `flags`, and `last_error`
/// to record what happened; they MUST NOT touch the input span.
struct MappingContext {
    // Inputs (read-only for stages)
    std::string_view read_id;
    std::span<const uint8_t> read_seq;     // 2-bit or ASCII; opaque here
    bool llm_fallback_enabled = false;     // CLI flag --llm-fallback

    // Catalog locus identity (for diagnostics / logging)
    std::string_view locus_id;

    // Mutable output state
    std::vector<CandidatePlacement> placements;
    std::vector<std::string> flags;        // free-form warning tags
    std::string last_error;                // populated by failing stages
    StageId last_stage_attempted{};        // most recent stage that ran
    StageId resolved_by{};                 // stage that returned success

    /// Convenience: did any stage succeed?
    bool resolved() const noexcept { return !placements.empty(); }

    /// Reset placements but keep flags (used between stages).
    void clear_placements() noexcept { placements.clear(); }
};

/// Aggregate result for the whole chain.
struct FallbackOutcome {
    bool success = false;
    StageId resolved_by{};
    uint32_t stages_run = 0;
    std::vector<std::string> flags;
    std::string last_error;
};

}  // namespace llmap::mapper::fallback
