#pragma once

#include "agent_types.h"

#include <chrono>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace llmap::claude_agent {

struct Message {
    std::string role;
    std::string content;
};

struct ToolCall {
    std::string id;
    std::string name;
    std::string input;
};

struct ToolResult {
    std::string tool_use_id;
    std::string output;
    bool is_error{false};
};

struct ConversationTurn {
    std::vector<Message> messages;
    std::vector<ToolCall> tool_calls;
    std::string stop_reason;
    size_t input_tokens{0};
    size_t output_tokens{0};
};

using ToolHandler = std::function<std::string(const std::string& name, const std::string& input)>;

class AnthropicClient {
public:
    explicit AnthropicClient(AgentConfig config);
    ~AnthropicClient();

    AnthropicClient(const AnthropicClient&) = delete;
    AnthropicClient& operator=(const AnthropicClient&) = delete;
    AnthropicClient(AnthropicClient&&) noexcept;
    AnthropicClient& operator=(AnthropicClient&&) noexcept;

    void SetToolHandler(ToolHandler handler);

    std::future<ConversationTurn> SendAsync(
        std::string_view system_prompt,
        std::vector<Message> messages
    );

    ConversationTurn Send(
        std::string_view system_prompt,
        std::vector<Message> messages
    );

    ConversationTurn RunConversation(
        std::string_view system_prompt,
        std::string_view initial_user_message,
        size_t max_turns = 100
    );

    double EstimateCostUsd() const;
    size_t TotalToolCalls() const;

    bool HasApiKey() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class TokenBucket {
public:
    TokenBucket(size_t tokens_per_minute, size_t burst = 0);

    bool TryAcquire(size_t tokens = 1);
    void Acquire(size_t tokens = 1);
    std::chrono::milliseconds WaitTime(size_t tokens = 1) const;

private:
    void Refill();

    mutable std::mutex mutex_;
    size_t tokens_per_minute_;
    size_t burst_;
    double tokens_;
    std::chrono::steady_clock::time_point last_refill_;
};

std::string GetApiKeyFromEnv();

}  // namespace llmap::claude_agent
