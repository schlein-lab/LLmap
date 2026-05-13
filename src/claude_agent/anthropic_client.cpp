#include "anthropic_client.h"

#include <atomic>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <thread>

namespace llmap::claude_agent {

std::string GetApiKeyFromEnv() {
    const char* key = std::getenv("ANTHROPIC_API_KEY");
    return key ? std::string(key) : std::string{};
}

TokenBucket::TokenBucket(size_t tokens_per_minute, size_t burst)
    : tokens_per_minute_(tokens_per_minute)
    , burst_(burst > 0 ? burst : tokens_per_minute)
    , tokens_(static_cast<double>(burst_))
    , last_refill_(std::chrono::steady_clock::now()) {
}

void TokenBucket::Refill() {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        now - last_refill_).count();
    double tokens_to_add = (static_cast<double>(elapsed) / 60000.0)
                           * static_cast<double>(tokens_per_minute_);
    tokens_ = std::min(static_cast<double>(burst_), tokens_ + tokens_to_add);
    last_refill_ = now;
}

bool TokenBucket::TryAcquire(size_t tokens) {
    std::lock_guard<std::mutex> lock(mutex_);
    Refill();
    if (tokens_ >= static_cast<double>(tokens)) {
        tokens_ -= static_cast<double>(tokens);
        return true;
    }
    return false;
}

void TokenBucket::Acquire(size_t tokens) {
    while (!TryAcquire(tokens)) {
        auto wait = WaitTime(tokens);
        std::this_thread::sleep_for(wait);
    }
}

std::chrono::milliseconds TokenBucket::WaitTime(size_t tokens) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (tokens_ >= static_cast<double>(tokens)) {
        return std::chrono::milliseconds(0);
    }
    double needed = static_cast<double>(tokens) - tokens_;
    double minutes_needed = needed / static_cast<double>(tokens_per_minute_);
    return std::chrono::milliseconds(
        static_cast<int64_t>(minutes_needed * 60000.0) + 1);
}

struct AnthropicClient::Impl {
    AgentConfig config;
    ToolHandler tool_handler;
    TokenBucket rate_limiter;
    std::atomic<size_t> total_input_tokens{0};
    std::atomic<size_t> total_output_tokens{0};
    std::atomic<size_t> total_tool_calls{0};

    explicit Impl(AgentConfig cfg)
        : config(std::move(cfg))
        , rate_limiter(cfg.rate_limit_rpm) {
    }

    ConversationTurn CallApi(
        std::string_view system_prompt,
        const std::vector<Message>& messages
    ) {
        rate_limiter.Acquire();

        ConversationTurn turn{};

        if (config.api_key.empty()) {
            turn.stop_reason = "error";
            return turn;
        }

        // Build CLI command for claude subprocess
        std::ostringstream cmd;
        cmd << "claude -p ";

        // Escape system prompt for shell
        std::string escaped_prompt;
        for (char c : system_prompt) {
            if (c == '\'') {
                escaped_prompt += "'\\''";
            } else {
                escaped_prompt += c;
            }
        }

        // Build full prompt with conversation history
        std::string full_prompt = std::string(system_prompt);
        full_prompt += "\n\n";
        for (const auto& msg : messages) {
            full_prompt += "[" + msg.role + "]: " + msg.content + "\n";
        }

        // For now, return a simulated response
        // Real implementation would use libcurl or subprocess
        turn.stop_reason = "end_turn";
        turn.input_tokens = full_prompt.size() / 4;  // rough estimate
        turn.output_tokens = 100;
        turn.messages.push_back({"assistant", "Agent session completed."});

        total_input_tokens += turn.input_tokens;
        total_output_tokens += turn.output_tokens;

        return turn;
    }
};

AnthropicClient::AnthropicClient(AgentConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {
}

AnthropicClient::~AnthropicClient() = default;

AnthropicClient::AnthropicClient(AnthropicClient&&) noexcept = default;
AnthropicClient& AnthropicClient::operator=(AnthropicClient&&) noexcept = default;

void AnthropicClient::SetToolHandler(ToolHandler handler) {
    impl_->tool_handler = std::move(handler);
}

std::future<ConversationTurn> AnthropicClient::SendAsync(
    std::string_view system_prompt,
    std::vector<Message> messages
) {
    return std::async(std::launch::async, [this, sp = std::string(system_prompt),
                                           msgs = std::move(messages)]() {
        return impl_->CallApi(sp, msgs);
    });
}

ConversationTurn AnthropicClient::Send(
    std::string_view system_prompt,
    std::vector<Message> messages
) {
    return impl_->CallApi(system_prompt, messages);
}

ConversationTurn AnthropicClient::RunConversation(
    std::string_view system_prompt,
    std::string_view initial_user_message,
    size_t max_turns
) {
    std::vector<Message> messages;
    messages.push_back({"user", std::string(initial_user_message)});

    ConversationTurn last_turn;
    for (size_t i = 0; i < max_turns; ++i) {
        last_turn = Send(system_prompt, messages);

        if (last_turn.stop_reason == "end_turn"
            || last_turn.stop_reason == "error") {
            break;
        }

        if (last_turn.stop_reason == "tool_use" && impl_->tool_handler) {
            for (const auto& tc : last_turn.tool_calls) {
                ++impl_->total_tool_calls;
                std::string result = impl_->tool_handler(tc.name, tc.input);
                messages.push_back({"assistant", "Tool call: " + tc.name});
                messages.push_back({"user", "Tool result: " + result});
            }
        } else {
            for (const auto& msg : last_turn.messages) {
                messages.push_back(msg);
            }
            break;
        }
    }

    return last_turn;
}

double AnthropicClient::EstimateCostUsd() const {
    constexpr double kInputPricePerMToken = 3.0;
    constexpr double kOutputPricePerMToken = 15.0;

    double input_cost = (static_cast<double>(impl_->total_input_tokens)
                         / 1000000.0) * kInputPricePerMToken;
    double output_cost = (static_cast<double>(impl_->total_output_tokens)
                          / 1000000.0) * kOutputPricePerMToken;
    return input_cost + output_cost;
}

size_t AnthropicClient::TotalToolCalls() const {
    return impl_->total_tool_calls;
}

bool AnthropicClient::HasApiKey() const {
    return !impl_->config.api_key.empty();
}

}  // namespace llmap::claude_agent
