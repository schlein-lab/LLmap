#include <gtest/gtest.h>

#include "ai/foundation_embedder.h"

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

namespace llmap::ai {
namespace {

// Helper to get test model path relative to project root
std::filesystem::path GetTestModelPath() {
    // Try different relative paths depending on where tests are run from
    std::vector<std::filesystem::path> candidates = {
        "tests/data/test_dna_embedder.onnx",
        "../tests/data/test_dna_embedder.onnx",
        "../../tests/data/test_dna_embedder.onnx",
    };
    for (const auto& p : candidates) {
        if (std::filesystem::exists(p)) {
            return p;
        }
    }
    return "tests/data/test_dna_embedder.onnx";
}

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

// Tests that require ONNX Runtime to be available
// These are conditional on runtime availability

class RealModelEmbedderTest : public ::testing::Test {
protected:
    void SetUp() override {
        model_path_ = GetTestModelPath();
        onnx_available_ = IsOnnxRuntimeAvailable() &&
                          std::filesystem::exists(model_path_);
    }

    std::filesystem::path model_path_;
    bool onnx_available_ = false;
};

TEST_F(RealModelEmbedderTest, CreateWithValidModelReturnsEmbedder) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    if (!embedder) {
        // Debug info for failure
        FAIL() << "Failed to create embedder. Model path: " << model_path_
               << " exists: " << std::filesystem::exists(model_path_)
               << " ONNX Runtime available: " << IsOnnxRuntimeAvailable();
    }
    EXPECT_TRUE(embedder->IsReady());
}

TEST_F(RealModelEmbedderTest, EmbedSingleSequenceProducesCorrectDimension) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result = embedder->Embed("ACGTACGTACGT");
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->embedding.size(), 256);
    EXPECT_GT(result->inference_time_us, 0.0f);
}

TEST_F(RealModelEmbedderTest, EmbeddingIsNormalized) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result = embedder->Embed("ACGTACGTACGT");
    ASSERT_TRUE(result.has_value());

    // The test model normalizes embeddings to unit length
    float norm_sq = 0.0f;
    for (float v : result->embedding) {
        norm_sq += v * v;
    }
    float norm = std::sqrt(norm_sq);
    EXPECT_NEAR(norm, 1.0f, 0.01f);
}

TEST_F(RealModelEmbedderTest, DifferentSequencesProduceDifferentEmbeddings) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result1 = embedder->Embed("AAAAAAAAAA");
    auto result2 = embedder->Embed("TTTTTTTTTT");

    ASSERT_TRUE(result1.has_value());
    ASSERT_TRUE(result2.has_value());

    // Embeddings should differ
    float diff_sum = 0.0f;
    for (size_t i = 0; i < result1->embedding.size(); ++i) {
        diff_sum += std::abs(result1->embedding[i] - result2->embedding[i]);
    }
    EXPECT_GT(diff_sum, 0.1f);  // Non-trivial difference
}

TEST_F(RealModelEmbedderTest, BatchEmbeddingWorks) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::vector<std::string_view> sequences = {
        "ACGTACGTACGT",
        "TTTTTTTTTTTT",
        "GGGGGGGGGGGG",
        "CCCCCCCCCCCC",
    };

    auto result = embedder->EmbedBatch(sequences);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->embeddings.size(), 4);

    for (const auto& emb : result->embeddings) {
        EXPECT_EQ(emb.size(), 256);
    }

    EXPECT_GT(result->total_inference_time_us, 0.0f);
    EXPECT_GT(result->avg_time_per_read_us, 0.0f);
}

TEST_F(RealModelEmbedderTest, LongSequenceIsTruncated) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    // Create a sequence longer than max_sequence_length
    std::string long_seq(1000, 'A');
    auto result = embedder->Embed(long_seq);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->embedding.size(), 256);
}

TEST_F(RealModelEmbedderTest, EmbedBatchIntoWorks) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::vector<std::string_view> sequences = {
        "ACGTACGT",
        "TTTTTTTT",
    };

    std::vector<float> output(2 * 256);
    bool success = embedder->EmbedBatchInto(sequences, output);
    EXPECT_TRUE(success);

    // Verify first embedding is not all zeros
    float sum1 = 0.0f;
    for (size_t i = 0; i < 256; ++i) {
        sum1 += std::abs(output[i]);
    }
    EXPECT_GT(sum1, 0.0f);
}

TEST_F(RealModelEmbedderTest, ModelInfoAccessors) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    EXPECT_EQ(embedder->ModelName(), "test_dna_embedder");
    EXPECT_FALSE(embedder->ModelVersion().empty());
    EXPECT_FALSE(embedder->IsGpuEnabled());
    EXPECT_EQ(embedder->ActiveProvider(), ExecutionProvider::CPU);
}

TEST_F(RealModelEmbedderTest, InvalidSequenceReturnsNullopt) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    // Sequence with invalid characters
    auto result = embedder->Embed("ACGT123XYZ");
    EXPECT_FALSE(result.has_value());
}

TEST_F(RealModelEmbedderTest, ConsistentEmbeddingsForSameSequence) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    EmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;

    auto embedder = FoundationEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result1 = embedder->Embed("ACGTACGTACGT");
    auto result2 = embedder->Embed("ACGTACGTACGT");

    ASSERT_TRUE(result1.has_value());
    ASSERT_TRUE(result2.has_value());

    // Same sequence should produce identical embeddings
    for (size_t i = 0; i < result1->embedding.size(); ++i) {
        EXPECT_FLOAT_EQ(result1->embedding[i], result2->embedding[i]);
    }
}

}  // namespace
}  // namespace llmap::ai
