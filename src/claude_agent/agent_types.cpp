#include "agent_types.h"

#include <algorithm>
#include <cctype>

namespace llmap::claude_agent {

std::optional<AgentMode> ParseAgentMode(std::string_view s) {
    std::string lower;
    lower.reserve(s.size());
    for (char c : s) {
        lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }

    if (lower == "off" || lower == "0") {
        return AgentMode::Off;
    }
    if (lower == "index-only" || lower == "indexonly" || lower == "1") {
        return AgentMode::IndexOnly;
    }
    if (lower == "sample-aware" || lower == "sampleaware" || lower == "2") {
        return AgentMode::SampleAware;
    }
    if (lower == "self-healing" || lower == "selfhealing" || lower == "3") {
        return AgentMode::SelfHealing;
    }
    if (lower == "research" || lower == "4") {
        return AgentMode::Research;
    }
    return std::nullopt;
}

}  // namespace llmap::claude_agent
