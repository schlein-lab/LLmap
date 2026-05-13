#include <gtest/gtest.h>

#include "ai/foundation_embedder.h"

#include <filesystem>
#include <string>
#include <vector>

namespace llmap::ai {
namespace {

class FoundationEmbedderTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default config for tests
        config_.model_path = "/nonexistent/model.onnx";
        config_.provider = ExecutionProvider::CPU;
        config_.embedding_dim = 256;
        config_.max_sequence_length = 512;
        config_.batch_size = 1024;
    }

    EmbedderConfig config_;
};

// Configuration tests

TEST_F(FoundationEmbedderTest, DefaultConfigValues) {
    EmbedderConfig cfg;
    EXPECT_EQ(cfg.provider, ExecutionProvider::CPU);
    EXPECT_EQ(cfg.device_id, 0);
    EXPECT_EQ(cfg.batch_size, 1024);
    EXPECT_EQ(cfg.embedding_dim, 256);
    EXPECT_EQ(cfg.max_sequence_length, 512);
    EXPECT_TRUE(cfg.enable_memory_pattern);
    EXPECT_EQ(cfg.intra_op_threads, 0);
    EXPECT_EQ(cfg.inter_op_threads, 0);
}

TEST_F(FoundationEmbedderTest, ConfigCanBeModified) {
    config_.provider = ExecutionProvider::CUDA;
    config_.device_id = 1;
    config_.batch_size = 2048;
    config_.embedding_dim = 512;
    config_.max_sequence_length = 1024;

    EXPECT_EQ(config_.provider, ExecutionProvider::CUDA);
    EXPECT_EQ(config_.device_id, 1);
    EXPECT_EQ(config_.batch_size, 2048);
    EXPECT_EQ(config_.embedding_dim, 512);
    EXPECT_EQ(config_.max_sequence_length, 1024);
}

// Factory tests

TEST_F(FoundationEmbedderTest, CreateReturnsNullptrForNonexistentModel) {
    auto embedder = FoundationEmbedder::Create(config_);
    // Should return nullptr because model file doesn't exist
    // (or ONNX Runtime isn't available)
    EXPECT_EQ(embedder, nullptr);
}

TEST_F(FoundationEmbedderTest, CreateWithEmptyPathReturnsNullptr) {
    config_.model_path = "";
    auto embedder = FoundationEmbedder::Create(config_);
    EXPECT_EQ(embedder, nullptr);
}

// Utility function tests

TEST_F(FoundationEmbedderTest, IsOnnxRuntimeAvailableReturnsBool) {
    // Just verify it returns a valid bool without crashing
    bool available = IsOnnxRuntimeAvailable();
    EXPECT_TRUE(available || !available);  // Always passes, just tests no crash
}

TEST_F(FoundationEmbedderTest, IsCudaEPAvailableReturnsBool) {
    bool available = IsCudaEPAvailable();
    EXPECT_TRUE(available || !available);
}

TEST_F(FoundationEmbedderTest, IsTensorRTEPAvailableReturnsBool) {
    bool available = IsTensorRTEPAvailable();
    EXPECT_TRUE(available || !available);
}

TEST_F(FoundationEmbedderTest, ListAvailableProvidersReturnsNonEmpty) {
    auto providers = ListAvailableProviders();
    // Should always have at least CPU provider (or stub message)
    EXPECT_FALSE(providers.empty());
}

// ExecutionProvider enum tests

TEST_F(FoundationEmbedderTest, ExecutionProviderEnumValues) {
    EXPECT_EQ(static_cast<int>(ExecutionProvider::CPU), 0);
    EXPECT_EQ(static_cast<int>(ExecutionProvider::CUDA), 1);
    EXPECT_EQ(static_cast<int>(ExecutionProvider::TensorRT), 2);
}

// EmbeddingResult struct tests

TEST_F(FoundationEmbedderTest, EmbeddingResultDefaultConstruction) {
    EmbeddingResult result{};
    EXPECT_TRUE(result.embedding.empty());
    EXPECT_EQ(result.inference_time_us, 0.0f);
}

TEST_F(FoundationEmbedderTest, EmbeddingResultCanHoldData) {
    EmbeddingResult result{
        .embedding = {0.1f, 0.2f, 0.3f, 0.4f},
        .inference_time_us = 123.45f,
    };
    EXPECT_EQ(result.embedding.size(), 4);
    EXPECT_FLOAT_EQ(result.embedding[0], 0.1f);
    EXPECT_FLOAT_EQ(result.inference_time_us, 123.45f);
}

// BatchEmbeddingResult struct tests

TEST_F(FoundationEmbedderTest, BatchEmbeddingResultDefaultConstruction) {
    BatchEmbeddingResult result{};
    EXPECT_TRUE(result.embeddings.empty());
    EXPECT_EQ(result.total_inference_time_us, 0.0f);
    EXPECT_EQ(result.avg_time_per_read_us, 0.0f);
}

TEST_F(FoundationEmbedderTest, BatchEmbeddingResultCanHoldMultipleEmbeddings) {
    BatchEmbeddingResult result{
        .embeddings = {
            {0.1f, 0.2f, 0.3f},
            {0.4f, 0.5f, 0.6f},
            {0.7f, 0.8f, 0.9f},
        },
        .total_inference_time_us = 300.0f,
        .avg_time_per_read_us = 100.0f,
    };
    EXPECT_EQ(result.embeddings.size(), 3);
    EXPECT_EQ(result.embeddings[0].size(), 3);
    EXPECT_FLOAT_EQ(result.embeddings[1][1], 0.5f);
    EXPECT_FLOAT_EQ(result.total_inference_time_us, 300.0f);
}

// Config path tests

TEST_F(FoundationEmbedderTest, ModelPathSupportsFilesystemPath) {
    std::filesystem::path model_path = "/opt/models/caduceus_ph_v1.onnx";
    config_.model_path = model_path;
    EXPECT_EQ(config_.model_path.filename().string(), "caduceus_ph_v1.onnx");
    EXPECT_EQ(config_.model_path.extension().string(), ".onnx");
}

TEST_F(FoundationEmbedderTest, ModelPathSupportsRelativePath) {
    config_.model_path = "models/foundation/caduceus.onnx";
    EXPECT_EQ(config_.model_path.filename().string(), "caduceus.onnx");
}

// Thread configuration tests

TEST_F(FoundationEmbedderTest, ThreadConfigurationDefaults) {
    EmbedderConfig cfg;
    EXPECT_EQ(cfg.intra_op_threads, 0);  // 0 = auto-detect
    EXPECT_EQ(cfg.inter_op_threads, 0);
}

TEST_F(FoundationEmbedderTest, ThreadConfigurationCanBeSet) {
    config_.intra_op_threads = 4;
    config_.inter_op_threads = 2;
    EXPECT_EQ(config_.intra_op_threads, 4);
    EXPECT_EQ(config_.inter_op_threads, 2);
}

// Memory pattern configuration

TEST_F(FoundationEmbedderTest, MemoryPatternDefaultsToEnabled) {
    EmbedderConfig cfg;
    EXPECT_TRUE(cfg.enable_memory_pattern);
}

TEST_F(FoundationEmbedderTest, MemoryPatternCanBeDisabled) {
    config_.enable_memory_pattern = false;
    EXPECT_FALSE(config_.enable_memory_pattern);
}

}  // namespace
}  // namespace llmap::ai
