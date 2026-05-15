// LLmap -- checkpoint dispatcher implementation.

#include "checkpoint/checkpoint_dispatcher.h"
#include "checkpoint/checkpoint_prompts.h"
#include "claude_agent/pipeline_agent.h"

namespace llmap::checkpoint {

CheckpointDispatcher::CheckpointDispatcher(LlmMode mode,
                                           CheckpointCache* cache,
                                           claude_agent::PipelineAgent* agent)
    : mode_(mode), cache_(cache), agent_(agent), tools_() {}

AgentDecision CheckpointDispatcher::FallbackDecision(const std::string& reason) {
    AgentDecision d;
    d.consulted = false;
    d.fallback_used = true;
    d.reasoning = reason;
    ++stats_.fallback_used;
    return d;
}

std::optional<AgentDecision> CheckpointDispatcher::AskAgent(
    CheckpointType type, const CheckpointContext& ctx) {

    if (!agent_) return std::nullopt;

    // The PipelineAgent today only exposes diagnostic / reporter sessions.
    // Wiring a synchronous "consult" call against the Anthropic client is
    // a follow-up: for the scaffold, we treat "agent present but no
    // checkpoint-shaped session yet" as a controlled fallback. When the
    // agent gains a checkpoint session (e.g. claude_agent::CheckpointSession)
    // this branch is the one to extend.
    //
    // The prompt is built here so it is exercised in tests and so the
    // hook-in point is unmistakable.
    (void)BuildPrompt(type, ctx);
    (void)SystemPrompt();
    return std::nullopt;
}

AgentDecision CheckpointDispatcher::Consult(CheckpointType type,
                                            const CheckpointContext& ctx) {
    ++stats_.by_type[static_cast<size_t>(type)];

    if (mode_ == LlmMode::Off) {
        return FallbackDecision("llm-mode=off");
    }

    // Cache lookup
    std::string key;
    if (cache_) {
        key = MakeCacheKey(type, ctx);
        if (auto hit = cache_->Lookup(key)) {
            ++stats_.cache_hits;
            ++stats_.consulted;
            return *hit;
        }
        ++stats_.cache_misses;
    }

    // Agent path
    auto agent_decision = AskAgent(type, ctx);
    if (agent_decision) {
        ++stats_.consulted;
        agent_decision->consulted = true;
        agent_decision->fallback_used = false;
        if (cache_) {
            cache_->Store(key, *agent_decision);
        }
        return *agent_decision;
    }

    if (mode_ == LlmMode::Required) {
        throw CheckpointAgentUnavailable(
            "llm-mode=required but agent unreachable for checkpoint "
            + std::string(CheckpointTypeName(type)));
    }

    // Auto mode: fall back silently.
    return FallbackDecision("agent unreachable / not yet wired");
}

}  // namespace llmap::checkpoint
