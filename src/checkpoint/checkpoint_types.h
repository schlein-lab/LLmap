// LLmap -- Layer 3 active LLM checkpoint types.
//
// A "checkpoint" is a point inside the mapping pipeline where the chain DP /
// extension may consult a live Claude agent because deterministic scoring is
// not confident enough. The dispatcher routes these requests through cache +
// agent and returns an AgentDecision that the caller can fold back into its
// mapping decision (multi-position wave, parameter override, free-form
// "special finding" annotation).
//
// Layer priority -- see annot/annot_types.h -- AgentDecision overrides
// SpecificLocus which overrides Taxonomy which overrides Default.

#pragma once

#include "annot/annot_types.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace llmap::checkpoint {

// Why was the agent invoked at this point in the pipeline?
enum class CheckpointType : uint8_t {
    AmbiguousChain        = 0,  // top-N chains within 5% of best
    UnknownRegion         = 1,  // window has no annotation match
    ParalogDisambiguation = 2,  // region flagged require_psv_disambig
    SDExpansion           = 3,  // read length > window size by > 10%
    NovelInsertion        = 4   // big soft-clip after extension
};

constexpr int kCheckpointTypeCount = 5;

constexpr const char* CheckpointTypeName(CheckpointType t) {
    switch (t) {
        case CheckpointType::AmbiguousChain:        return "ambiguous-chain";
        case CheckpointType::UnknownRegion:         return "unknown-region";
        case CheckpointType::ParalogDisambiguation: return "paralog-disambig";
        case CheckpointType::SDExpansion:           return "sd-expansion";
        case CheckpointType::NovelInsertion:        return "novel-insertion";
    }
    return "unknown";
}

// Three operating modes mirroring `--llm=auto|off|required` from
// knowledge/EXTENDING.md.
enum class LlmMode : uint8_t {
    Off      = 0,  // never invoke agent, fully deterministic
    Auto     = 1,  // try cache + agent, fall back silently
    Required = 2   // cache + agent, throw on failure
};

constexpr const char* LlmModeName(LlmMode m) {
    switch (m) {
        case LlmMode::Off:      return "off";
        case LlmMode::Auto:     return "auto";
        case LlmMode::Required: return "required";
    }
    return "unknown";
}

// Context handed to the dispatcher at a checkpoint. Keep this serialisable
// so the cache key is deterministic.
struct CheckpointContext {
    std::string read_id;
    std::string read_seq;          // optional, <= 1 kb truncation for prompt size
    uint32_t ref_id = 0;
    std::vector<std::pair<uint32_t, int32_t>> candidate_positions;  // (pos, chain_score)
    std::string region_name;       // from annotation if available
    std::string user_note;         // free-form, optional
};

// Result of the consultation. consulted=false means "didn't ask anyone" --
// either Off mode or cache-only path produced nothing. fallback_used=true
// means classical scoring should make the decision.
struct AgentDecision {
    bool consulted = false;
    bool fallback_used = false;
    std::string reasoning;         // human-readable explanation
    std::vector<std::pair<uint32_t, float>> wave;  // (pos, amplitude)
    std::string special_finding;   // "novel SD", "VDJ rearrangement", etc.
    annot::ParamOverride override;
};

}  // namespace llmap::checkpoint
