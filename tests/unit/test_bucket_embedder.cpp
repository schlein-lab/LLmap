#include <gtest/gtest.h>

#include "ai/bucket_embedder.h"

#include <cmath>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace llmap::ai {
namespace {

// Helper to get test model path
std::filesystem::path GetTestModelPath() {
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

// Generate random DNA sequence
std::string GenerateRandomSequence(size_t length, uint32_t seed = 42) {
    static const char nucleotides[] = "ACGT";
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += nucleotides[dist(gen)];
    }
    return seq;
}

class BucketEmbedderTest : public ::testing::Test {
protected:
    void SetUp() override {
        config_.model_path = "/nonexistent/model.onnx";
        config_.provider = ExecutionProvider::CPU;
        config_.embedding_dim = 256;
        config_.max_sequence_length = 8192;
        config_.batch_size = 64;
        config_.chunk_size = 1024;
        config_.chunk_overlap = 128;
    }

    BucketEmbedderConfig config_;
};

// Configuration tests

TEST_F(BucketEmbedderTest, DefaultConfigValues) {
    BucketEmbedderConfig cfg;
    EXPECT_EQ(cfg.provider, ExecutionProvider::CPU);
    EXPECT_EQ(cfg.device_id, 0);
    EXPECT_EQ(cfg.batch_size, 64);
    EXPECT_EQ(cfg.embedding_dim, 256);
    EXPECT_EQ(cfg.max_sequence_length, 8192);
    EXPECT_EQ(cfg.chunk_size, 1024);
    EXPECT_EQ(cfg.chunk_overlap, 128);
    EXPECT_TRUE(cfg.enable_memory_pattern);
    EXPECT_EQ(cfg.intra_op_threads, 0);
    EXPECT_EQ(cfg.inter_op_threads, 0);
    EXPECT_TRUE(cfg.normalize_output);
}

TEST_F(BucketEmbedderTest, ConfigCanBeModified) {
    config_.provider = ExecutionProvider::CUDA;
    config_.device_id = 1;
    config_.batch_size = 128;
    config_.embedding_dim = 512;
    config_.max_sequence_length = 16384;
    config_.normalize_output = false;

    EXPECT_EQ(config_.provider, ExecutionProvider::CUDA);
    EXPECT_EQ(config_.device_id, 1);
    EXPECT_EQ(config_.batch_size, 128);
    EXPECT_EQ(config_.embedding_dim, 512);
    EXPECT_EQ(config_.max_sequence_length, 16384);
    EXPECT_FALSE(config_.normalize_output);
}

// Factory tests

TEST_F(BucketEmbedderTest, CreateReturnsNullptrForNonexistentModel) {
    auto embedder = BucketEmbedder::Create(config_);
    EXPECT_EQ(embedder, nullptr);
}

TEST_F(BucketEmbedderTest, CreateWithEmptyPathReturnsNullptr) {
    config_.model_path = "";
    auto embedder = BucketEmbedder::Create(config_);
    EXPECT_EQ(embedder, nullptr);
}

// BucketEmbeddingStats tests

TEST_F(BucketEmbedderTest, StatsDefaultConstruction) {
    BucketEmbeddingStats stats{};
    EXPECT_EQ(stats.total_buckets, 0);
    EXPECT_EQ(stats.total_nucleotides, 0);
    EXPECT_EQ(stats.total_chunks, 0);
    EXPECT_EQ(stats.total_time_ms, 0.0f);
    EXPECT_EQ(stats.avg_time_per_bucket_ms, 0.0f);
    EXPECT_EQ(stats.throughput_buckets_per_sec, 0.0f);
    EXPECT_EQ(stats.throughput_mb_per_sec, 0.0f);
}

TEST_F(BucketEmbedderTest, StatsCanHoldData) {
    BucketEmbeddingStats stats{
        .total_buckets = 100,
        .total_nucleotides = 5000000,
        .total_chunks = 500,
        .total_time_ms = 1000.0f,
        .avg_time_per_bucket_ms = 10.0f,
        .throughput_buckets_per_sec = 100.0f,
        .throughput_mb_per_sec = 5.0f,
    };
    EXPECT_EQ(stats.total_buckets, 100);
    EXPECT_EQ(stats.total_nucleotides, 5000000);
    EXPECT_EQ(stats.total_chunks, 500);
    EXPECT_FLOAT_EQ(stats.total_time_ms, 1000.0f);
    EXPECT_FLOAT_EQ(stats.throughput_mb_per_sec, 5.0f);
}

// Utility function tests

TEST_F(BucketEmbedderTest, CosineSimilarityIdentical) {
    std::vector<float> a = {1.0f, 0.0f, 0.0f, 0.0f};
    std::vector<float> b = {1.0f, 0.0f, 0.0f, 0.0f};
    float sim = CosineSimilarity(a, b);
    EXPECT_NEAR(sim, 1.0f, 1e-6f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityOrthogonal) {
    std::vector<float> a = {1.0f, 0.0f, 0.0f, 0.0f};
    std::vector<float> b = {0.0f, 1.0f, 0.0f, 0.0f};
    float sim = CosineSimilarity(a, b);
    EXPECT_NEAR(sim, 0.0f, 1e-6f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityOpposite) {
    std::vector<float> a = {1.0f, 0.0f};
    std::vector<float> b = {-1.0f, 0.0f};
    float sim = CosineSimilarity(a, b);
    EXPECT_NEAR(sim, -1.0f, 1e-6f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityNonUnit) {
    std::vector<float> a = {3.0f, 4.0f};  // magnitude 5
    std::vector<float> b = {6.0f, 8.0f};  // magnitude 10, same direction
    float sim = CosineSimilarity(a, b);
    EXPECT_NEAR(sim, 1.0f, 1e-6f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityEmpty) {
    std::vector<float> a = {};
    std::vector<float> b = {};
    float sim = CosineSimilarity(a, b);
    EXPECT_EQ(sim, 0.0f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityMismatchedSize) {
    std::vector<float> a = {1.0f, 0.0f};
    std::vector<float> b = {1.0f, 0.0f, 0.0f};
    float sim = CosineSimilarity(a, b);
    EXPECT_EQ(sim, 0.0f);
}

TEST_F(BucketEmbedderTest, CosineSimilarityMatrixSmall) {
    // 3 vectors of dimension 2
    std::vector<float> embeddings = {
        1.0f, 0.0f,   // vec 0: unit x
        0.0f, 1.0f,   // vec 1: unit y
        1.0f, 1.0f,   // vec 2: diagonal (not normalized)
    };

    std::vector<float> matrix(9);
    ComputeCosineSimilarityMatrix(embeddings, 3, 2, matrix);

    // Diagonal should be 1
    EXPECT_NEAR(matrix[0], 1.0f, 1e-6f);  // (0,0)
    EXPECT_NEAR(matrix[4], 1.0f, 1e-6f);  // (1,1)
    EXPECT_NEAR(matrix[8], 1.0f, 1e-6f);  // (2,2)

    // (0,1) and (1,0) should be 0 (orthogonal)
    EXPECT_NEAR(matrix[1], 0.0f, 1e-6f);
    EXPECT_NEAR(matrix[3], 0.0f, 1e-6f);

    // (0,2) and (2,0): (1,0) dot (1,1) / (1 * sqrt(2)) = 1/sqrt(2) ~ 0.707
    float expected = 1.0f / std::sqrt(2.0f);
    EXPECT_NEAR(matrix[2], expected, 1e-5f);
    EXPECT_NEAR(matrix[6], expected, 1e-5f);
}

// Config path tests

TEST_F(BucketEmbedderTest, ModelPathSupportsFilesystemPath) {
    std::filesystem::path model_path = "/opt/models/evo_distilled_v1.onnx";
    config_.model_path = model_path;
    EXPECT_EQ(config_.model_path.filename().string(), "evo_distilled_v1.onnx");
    EXPECT_EQ(config_.model_path.extension().string(), ".onnx");
}

TEST_F(BucketEmbedderTest, ChunkConfigValidation) {
    config_.chunk_size = 2048;
    config_.chunk_overlap = 256;
    EXPECT_GT(config_.chunk_size, config_.chunk_overlap);
}

// Tests that require ONNX Runtime to be available

class RealModelBucketEmbedderTest : public ::testing::Test {
protected:
    void SetUp() override {
        model_path_ = GetTestModelPath();
        onnx_available_ = IsOnnxRuntimeAvailable() &&
                          std::filesystem::exists(model_path_);
    }

    std::filesystem::path model_path_;
    bool onnx_available_ = false;
};

TEST_F(RealModelBucketEmbedderTest, CreateWithValidModelReturnsEmbedder) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;  // Use smaller for test model

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);
    EXPECT_TRUE(embedder->IsReady());
}

TEST_F(RealModelBucketEmbedderTest, EmbedBucketProducesCorrectDimension) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::string bucket_seq = GenerateRandomSequence(500, 42);
    auto result = embedder->EmbedBucket(bucket_seq);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->size(), 256);
}

TEST_F(RealModelBucketEmbedderTest, EmbeddingIsNormalized) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;
    config.normalize_output = true;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::string bucket_seq = GenerateRandomSequence(400, 123);
    auto result = embedder->EmbedBucket(bucket_seq);
    ASSERT_TRUE(result.has_value());

    float norm_sq = 0.0f;
    for (float v : *result) {
        norm_sq += v * v;
    }
    float norm = std::sqrt(norm_sq);
    EXPECT_NEAR(norm, 1.0f, 0.01f);
}

TEST_F(RealModelBucketEmbedderTest, DifferentBucketsProduceDifferentEmbeddings) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result1 = embedder->EmbedBucket(GenerateRandomSequence(300, 1));
    auto result2 = embedder->EmbedBucket(GenerateRandomSequence(300, 2));

    ASSERT_TRUE(result1.has_value());
    ASSERT_TRUE(result2.has_value());

    float diff_sum = 0.0f;
    for (size_t i = 0; i < result1->size(); ++i) {
        diff_sum += std::abs((*result1)[i] - (*result2)[i]);
    }
    EXPECT_GT(diff_sum, 0.1f);
}

TEST_F(RealModelBucketEmbedderTest, BatchEmbeddingWorks) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::vector<std::string> sequences_storage = {
        GenerateRandomSequence(200, 1),
        GenerateRandomSequence(300, 2),
        GenerateRandomSequence(250, 3),
        GenerateRandomSequence(150, 4),
    };

    std::vector<std::string_view> sequences;
    for (const auto& s : sequences_storage) {
        sequences.push_back(s);
    }

    auto result = embedder->EmbedBuckets(sequences);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->size(), 4);

    for (const auto& emb : *result) {
        EXPECT_EQ(emb.size(), 256);
    }
}

TEST_F(RealModelBucketEmbedderTest, EmbedBucketsIntoWithStats) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::vector<std::string> sequences_storage = {
        GenerateRandomSequence(200, 1),
        GenerateRandomSequence(300, 2),
    };

    std::vector<std::string_view> sequences;
    for (const auto& s : sequences_storage) {
        sequences.push_back(s);
    }

    std::vector<float> output(2 * 256);
    BucketEmbeddingStats stats;

    bool success = embedder->EmbedBucketsInto(sequences, output, &stats);
    EXPECT_TRUE(success);

    EXPECT_EQ(stats.total_buckets, 2);
    EXPECT_EQ(stats.total_nucleotides, 500);  // 200 + 300
    EXPECT_GT(stats.total_time_ms, 0.0f);
    EXPECT_GT(stats.avg_time_per_bucket_ms, 0.0f);
    EXPECT_GT(stats.throughput_buckets_per_sec, 0.0f);
}

TEST_F(RealModelBucketEmbedderTest, ConsistentEmbeddingsForSameBucket) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    std::string bucket = GenerateRandomSequence(400, 999);
    auto result1 = embedder->EmbedBucket(bucket);
    auto result2 = embedder->EmbedBucket(bucket);

    ASSERT_TRUE(result1.has_value());
    ASSERT_TRUE(result2.has_value());

    for (size_t i = 0; i < result1->size(); ++i) {
        EXPECT_FLOAT_EQ((*result1)[i], (*result2)[i]);
    }
}

TEST_F(RealModelBucketEmbedderTest, InvalidSequenceReturnsNullopt) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    auto result = embedder->EmbedBucket("ACGT123XYZ");
    EXPECT_FALSE(result.has_value());
}

TEST_F(RealModelBucketEmbedderTest, ModelInfoAccessors) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    EXPECT_EQ(embedder->ModelName(), "test_dna_embedder");
    EXPECT_FALSE(embedder->ModelVersion().empty());
    EXPECT_FALSE(embedder->IsGpuEnabled());
    EXPECT_EQ(embedder->ActiveProvider(), ExecutionProvider::CPU);
}

TEST_F(RealModelBucketEmbedderTest, SimilarSequencesHaveHighCosineSimilarity) {
    if (!onnx_available_) {
        GTEST_SKIP() << "ONNX Runtime not available or test model not found";
    }

    BucketEmbedderConfig config;
    config.model_path = model_path_;
    config.provider = ExecutionProvider::CPU;
    config.embedding_dim = 256;
    config.max_sequence_length = 512;

    auto embedder = BucketEmbedder::Create(config);
    ASSERT_NE(embedder, nullptr);

    // Two sequences with only a few mutations should be similar
    std::string seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    std::string seq2 = "ACGTACGTACGTACGTACGTACGTACGTACGT";  // identical

    auto emb1 = embedder->EmbedBucket(seq1);
    auto emb2 = embedder->EmbedBucket(seq2);

    ASSERT_TRUE(emb1.has_value());
    ASSERT_TRUE(emb2.has_value());

    float sim = CosineSimilarity(*emb1, *emb2);
    EXPECT_NEAR(sim, 1.0f, 1e-5f);
}

}  // namespace
}  // namespace llmap::ai
