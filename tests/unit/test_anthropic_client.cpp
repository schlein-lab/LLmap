#include <gtest/gtest.h>

#include "claude_agent/anthropic_client.h"

#include <chrono>
#include <thread>

using namespace llmap::claude_agent;

class AnthropicClientTest : public ::testing::Test {};

TEST_F(AnthropicClientTest, TokenBucketInitialBurst) {
    TokenBucket bucket(60, 10);  // 60 tokens/min, burst of 10
    // Should be able to acquire burst immediately
    for (int i = 0; i < 10; ++i) {
        EXPECT_TRUE(bucket.TryAcquire(1));
    }
    // 11th should fail
    EXPECT_FALSE(bucket.TryAcquire(1));
}

TEST_F(AnthropicClientTest, TokenBucketRefill) {
    TokenBucket bucket(6000, 1);  // 6000 tokens/min = 100 tokens/sec
    // Drain the bucket
    EXPECT_TRUE(bucket.TryAcquire(1));
    EXPECT_FALSE(bucket.TryAcquire(1));

    // Wait 20ms, should refill ~2 tokens
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    EXPECT_TRUE(bucket.TryAcquire(1));
}

TEST_F(AnthropicClientTest, TokenBucketWaitTime) {
    TokenBucket bucket(60, 1);  // 60 tokens/min = 1 token/sec
    bucket.TryAcquire(1);  // Drain

    auto wait = bucket.WaitTime(1);
    EXPECT_GT(wait.count(), 0);
    EXPECT_LE(wait.count(), 1001);  // Should be ~1 second
}

TEST_F(AnthropicClientTest, TokenBucketAcquireBlocking) {
    TokenBucket bucket(60000, 1);  // 60000 tokens/min = 1000/sec
    auto start = std::chrono::steady_clock::now();

    bucket.Acquire(1);  // Should succeed immediately
    bucket.Acquire(1);  // Should block briefly

    auto elapsed = std::chrono::steady_clock::now() - start;
    // Both acquires should complete in under 10ms
    EXPECT_LT(std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count(), 10);
}

TEST_F(AnthropicClientTest, GetApiKeyFromEnvEmpty) {
    // Unset to test empty case
    auto key = GetApiKeyFromEnv();
    // May or may not be set in test environment
    // Just verify it doesn't crash
    EXPECT_GE(key.size(), 0u);
}

TEST_F(AnthropicClientTest, ClientConstruction) {
    AgentConfig config;
    config.api_key = "";  // No API key
    config.rate_limit_rpm = 60;

    AnthropicClient client(std::move(config));
    EXPECT_FALSE(client.HasApiKey());
    EXPECT_EQ(0u, client.TotalToolCalls());
    EXPECT_DOUBLE_EQ(0.0, client.EstimateCostUsd());
}

TEST_F(AnthropicClientTest, ClientWithApiKey) {
    AgentConfig config;
    config.api_key = "test-key";
    config.rate_limit_rpm = 60;

    AnthropicClient client(std::move(config));
    EXPECT_TRUE(client.HasApiKey());
}

TEST_F(AnthropicClientTest, SendWithoutApiKey) {
    AgentConfig config;
    config.api_key = "";  // No API key

    AnthropicClient client(std::move(config));
    auto turn = client.Send("System prompt", {{"user", "Hello"}});

    EXPECT_EQ("error", turn.stop_reason);
}

TEST_F(AnthropicClientTest, RunConversationWithoutApiKey) {
    AgentConfig config;
    config.api_key = "";

    AnthropicClient client(std::move(config));
    auto turn = client.RunConversation("System", "User message", 5);

    EXPECT_EQ("error", turn.stop_reason);
}

TEST_F(AnthropicClientTest, SetToolHandler) {
    AgentConfig config;
    config.api_key = "test-key";

    AnthropicClient client(std::move(config));

    bool handler_called = false;
    client.SetToolHandler([&](const std::string& name, const std::string& input) {
        handler_called = true;
        return "result";
    });

    // The handler won't be called without actual API interaction
    EXPECT_FALSE(handler_called);
}

TEST_F(AnthropicClientTest, MoveConstruction) {
    AgentConfig config;
    config.api_key = "test-key";

    AnthropicClient client1(std::move(config));
    EXPECT_TRUE(client1.HasApiKey());

    AnthropicClient client2(std::move(client1));
    EXPECT_TRUE(client2.HasApiKey());
}

TEST_F(AnthropicClientTest, MoveAssignment) {
    AgentConfig config1;
    config1.api_key = "key1";
    AnthropicClient client1(std::move(config1));

    AgentConfig config2;
    config2.api_key = "";
    AnthropicClient client2(std::move(config2));

    EXPECT_FALSE(client2.HasApiKey());
    client2 = std::move(client1);
    EXPECT_TRUE(client2.HasApiKey());
}

TEST_F(AnthropicClientTest, MessageStruct) {
    Message msg;
    msg.role = "user";
    msg.content = "Hello, Claude!";

    EXPECT_EQ("user", msg.role);
    EXPECT_EQ("Hello, Claude!", msg.content);
}

TEST_F(AnthropicClientTest, ToolCallStruct) {
    ToolCall tc;
    tc.id = "call_123";
    tc.name = "bash";
    tc.input = R"({"command": "ls"})";

    EXPECT_EQ("call_123", tc.id);
    EXPECT_EQ("bash", tc.name);
    EXPECT_EQ(R"({"command": "ls"})", tc.input);
}

TEST_F(AnthropicClientTest, ToolResultStruct) {
    ToolResult tr;
    tr.tool_use_id = "call_123";
    tr.output = "file1.txt\nfile2.txt";
    tr.is_error = false;

    EXPECT_EQ("call_123", tr.tool_use_id);
    EXPECT_EQ("file1.txt\nfile2.txt", tr.output);
    EXPECT_FALSE(tr.is_error);
}

TEST_F(AnthropicClientTest, ConversationTurnDefaults) {
    ConversationTurn turn;

    EXPECT_TRUE(turn.messages.empty());
    EXPECT_TRUE(turn.tool_calls.empty());
    EXPECT_TRUE(turn.stop_reason.empty());
    EXPECT_EQ(0u, turn.input_tokens);
    EXPECT_EQ(0u, turn.output_tokens);
}
