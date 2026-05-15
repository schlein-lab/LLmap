// LLmap -- checkpoint dispatcher implementation.

#include "checkpoint/checkpoint_dispatcher.h"
#include "checkpoint/checkpoint_prompts.h"
#include "claude_agent/pipeline_agent.h"
#include "claude_agent/anthropic_client.h"

#include <iostream>

namespace llmap::checkpoint {

namespace {

// Minimal JSON-ish response parser tailored to the prompt template in
// checkpoint_prompts.cpp. Caller has freshly-built BuildPrompt() text and
// gets back the agent's reply; we look for "decision": "<pos>" and a
// "reasoning": "<text>" pair plus optional "wave": [...] amplitudes.
AgentDecision ParseAgentReply(std::string_view text) {
    AgentDecision d;
    d.consulted = true;
    d.fallback_used = false;

    auto find_field = [&](std::string_view key) -> std::string {
        auto pos = text.find(key);
        if (pos == std::string_view::npos) return {};
        pos = text.find(':', pos);
        if (pos == std::string_view::npos) return {};
        ++pos;
        while (pos < text.size() && (text[pos] == ' ' || text[pos] == '"')) ++pos;
        auto end = text.find('"', pos);
        if (end == std::string_view::npos) end = text.find(',', pos);
        if (end == std::string_view::npos) end = text.size();
        return std::string(text.substr(pos, end - pos));
    };

    d.reasoning = find_field("\"reasoning\"");
    d.special_finding = find_field("\"finding\"");
    auto decision_str = find_field("\"decision\"");
    // The decision string format is "chr:pos:amplitude,chr:pos:amplitude,...".
    size_t start = 0;
    while (start < decision_str.size()) {
        auto end = decision_str.find(',', start);
        if (end == std::string::npos) end = decision_str.size();
        auto seg = decision_str.substr(start, end - start);
        auto c1 = seg.find(':');
        auto c2 = seg.find(':', c1 == std::string::npos ? 0 : c1 + 1);
        if (c1 != std::string::npos && c2 != std::string::npos) {
            try {
                uint32_t pos = std::stoul(seg.substr(c1 + 1, c2 - c1 - 1));
                float amp = std::stof(seg.substr(c2 + 1));
                d.wave.push_back({pos, amp});
            } catch (...) {}
        }
        start = end + 1;
    }
    return d;
}

}  // namespace

CheckpointDispatcher::CheckpointDispatcher(LlmMode mode,
                                           CheckpointCache* cache,
                                           claude_agent::PipelineAgent* agent)
    : mode_(mode), cache_(cache), agent_(agent), client_(nullptr), tools_() {}

CheckpointDispatcher::CheckpointDispatcher(LlmMode mode,
                                           CheckpointCache* cache,
                                           claude_agent::AnthropicClient* client,
                                           DirectClientTag)
    : mode_(mode), cache_(cache), agent_(nullptr), client_(client), tools_() {}

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

    auto user_prompt = BuildPrompt(type, ctx);
    auto system = SystemPrompt();

    // Path 1: direct AnthropicClient (preferred for synchronous consults).
    if (client_ && client_->HasApiKey()) {
        try {
            std::vector<claude_agent::Message> msgs;
            claude_agent::Message m;
            m.role = "user";
            m.content = user_prompt;
            msgs.push_back(std::move(m));
            auto turn = client_->Send(system, std::move(msgs));
            // ConversationTurn.messages contains the assistant's reply
            // appended to the input messages. Grab the last assistant
            // message content.
            std::string reply;
            for (auto it = turn.messages.rbegin(); it != turn.messages.rend(); ++it) {
                if (it->role == "assistant") { reply = it->content; break; }
            }
            return ParseAgentReply(reply);
        } catch (const std::exception& e) {
            std::cerr << "[checkpoint] agent error: " << e.what() << "\n";
            return std::nullopt;
        }
    }

    // Path 2: PipelineAgent (long-running diagnostics, not yet suitable for
    // per-checkpoint sync calls — left as a fallback for forward compat).
    if (agent_) {
        // The PipelineAgent exposes only diagnostic/reporter sessions; a
        // CheckpointSession is the natural extension point.
        (void)agent_;
    }

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
