// LLmap -- prompt templates per CheckpointType.
//
// Each checkpoint kind has its own structured prompt that pins the agent
// to a JSON output shape the dispatcher can parse back into an
// AgentDecision. Prompts live inline in the .cpp; no external file
// dependency yet.

#pragma once

#include "checkpoint/checkpoint_types.h"

#include <string>

namespace llmap::checkpoint {

// Build the user-facing prompt for a checkpoint. Includes context fields
// and the JSON schema we expect the agent to respond with.
std::string BuildPrompt(CheckpointType type, const CheckpointContext& ctx);

// Shared system prompt used for every checkpoint consultation. Spells out
// the agent's role, available tools, and the JSON shape it must emit.
std::string SystemPrompt();

// Parse the agent's JSON response into an AgentDecision. Returns a decision
// with `fallback_used=true` if parsing fails so the caller can degrade
// gracefully.
AgentDecision ParseAgentResponse(const std::string& response);

}  // namespace llmap::checkpoint
