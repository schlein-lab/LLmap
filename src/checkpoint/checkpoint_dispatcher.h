// LLmap -- the Layer 3 consult dispatcher.
//
// The chain DP / extension stage calls Consult() at each checkpoint. The
// dispatcher implements the LlmMode policy:
//
//   Off      -> short-circuit, return consulted=false, fallback_used=true
//   Auto     -> cache -> agent; on agent failure, silent fallback
//   Required -> cache -> agent; on agent failure, throw runtime_error
//
// The dispatcher owns its statistics counters so the CLI can print an
// end-of-run summary.

#pragma once

#include "checkpoint/checkpoint_cache.h"
#include "checkpoint/checkpoint_tools.h"
#include "checkpoint/checkpoint_types.h"

#include <cstddef>
#include <stdexcept>

namespace llmap::claude_agent { class PipelineAgent; }

namespace llmap::checkpoint {

class CheckpointDispatcher {
public:
    struct Stats {
        size_t consulted = 0;
        size_t cache_hits = 0;
        size_t cache_misses = 0;
        size_t fallback_used = 0;
        size_t by_type[kCheckpointTypeCount] = {0};
    };

    // The dispatcher does not take ownership of cache or agent. `agent`
    // may be nullptr; that's how callers run in Off mode without paying
    // the cost of constructing a PipelineAgent.
    CheckpointDispatcher(LlmMode mode,
                         CheckpointCache* cache,
                         claude_agent::PipelineAgent* agent);

    // The hot-path entry point. Always returns a usable AgentDecision; the
    // caller checks `fallback_used` to decide whether to fall back to
    // classical scoring. May throw std::runtime_error in Required mode if
    // the agent isn't reachable.
    AgentDecision Consult(CheckpointType type, const CheckpointContext& ctx);

    const Stats& GetStats() const { return stats_; }
    LlmMode Mode() const { return mode_; }
    CheckpointToolRegistry& Tools() { return tools_; }
    const CheckpointToolRegistry& Tools() const { return tools_; }

private:
    AgentDecision FallbackDecision(const std::string& reason);

    // Try to consult the agent. Returns nullopt if the agent isn't
    // reachable / didn't reply. Implementation lives in the .cpp so the
    // header doesn't pull in PipelineAgent's full definition.
    std::optional<AgentDecision> AskAgent(CheckpointType type,
                                          const CheckpointContext& ctx);

    LlmMode mode_;
    CheckpointCache* cache_;
    claude_agent::PipelineAgent* agent_;
    CheckpointToolRegistry tools_;
    Stats stats_;
};

// Exception type so callers can distinguish a Required-mode failure from
// any other runtime error.
class CheckpointAgentUnavailable : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

}  // namespace llmap::checkpoint
