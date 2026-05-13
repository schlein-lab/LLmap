#include <gtest/gtest.h>

#include "self_interference/faiss_wrapper.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace llmap::self_interference {
namespace {

// Generate random normalized vectors
std::vector<float> GenerateRandomVectors(
    size_t num_vectors,
    size_t dim,
    uint32_t seed = 42)
{
    std::mt19937 gen(seed);
    std::normal_distribution<float> dist(0.0f, 1.0f);

    std::vector<float> vectors(num_vectors * dim);
    for (size_t i = 0; i < num_vectors; ++i) {
        float* vec = vectors.data() + i * dim;
        float norm = 0.0f;
        for (size_t j = 0; j < dim; ++j) {
            vec[j] = dist(gen);
            norm += vec[j] * vec[j];
        }
        norm = std::sqrt(norm);
        if (norm > 1e-8f) {
            for (size_t j = 0; j < dim; ++j) {
                vec[j] /= norm;
            }
        }
    }
    return vectors;
}

// Compute L2 distance between two vectors
float L2Distance(const float* a, const float* b, size_t dim) {
    float dist = 0.0f;
    for (size_t i = 0; i < dim; ++i) {
        float diff = a[i] - b[i];
        dist += diff * diff;
    }
    return dist;
}

// Brute-force k-NN for verification
std::vector<int64_t> BruteForceKnn(
    const std::vector<float>& database,
    const float* query,
    size_t num_vectors,
    size_t dim,
    size_t k)
{
    std::vector<std::pair<float, int64_t>> distances(num_vectors);
    for (size_t i = 0; i < num_vectors; ++i) {
        distances[i] = {L2Distance(database.data() + i * dim, query, dim), static_cast<int64_t>(i)};
    }
    std::partial_sort(distances.begin(), distances.begin() + k, distances.end());

    std::vector<int64_t> indices(k);
    for (size_t i = 0; i < k; ++i) {
        indices[i] = distances[i].second;
    }
    return indices;
}

class FaissWrapperTest : public ::testing::Test {
protected:
    void SetUp() override {
        config_.embedding_dim = 64;
        config_.nlist = 16;
        config_.nprobe = 4;
        config_.default_k = 10;
        config_.normalize_vectors = false;  // Vectors already normalized in tests
    }

    FaissIndexConfig config_;
};

// ========== Configuration Tests ==========

TEST_F(FaissWrapperTest, DefaultConfigValues) {
    FaissIndexConfig cfg;
    EXPECT_EQ(cfg.index_type, FaissIndexType::IVFFlat);
    EXPECT_EQ(cfg.provider, FaissProvider::CPU);
    EXPECT_EQ(cfg.gpu_device_id, 0);
    EXPECT_EQ(cfg.embedding_dim, 256);
    EXPECT_EQ(cfg.nlist, 1024);
    EXPECT_EQ(cfg.nprobe, 32);
    EXPECT_EQ(cfg.m_subquantizers, 8);
    EXPECT_EQ(cfg.nbits_per_idx, 8);
    EXPECT_EQ(cfg.training_sample_size, 0);
    EXPECT_EQ(cfg.default_k, 50);
    EXPECT_TRUE(cfg.normalize_vectors);
}

TEST_F(FaissWrapperTest, ConfigCanBeModified) {
    config_.index_type = FaissIndexType::FlatL2;
    config_.provider = FaissProvider::GPU;
    config_.embedding_dim = 128;
    config_.nlist = 256;
    config_.nprobe = 16;

    EXPECT_EQ(config_.index_type, FaissIndexType::FlatL2);
    EXPECT_EQ(config_.provider, FaissProvider::GPU);
    EXPECT_EQ(config_.embedding_dim, 128);
    EXPECT_EQ(config_.nlist, 256);
    EXPECT_EQ(config_.nprobe, 16);
}

// ========== Utility Function Tests ==========

TEST(FaissUtilsTest, IsFaissAvailableReturnsBoolean) {
    bool available = IsFaissAvailable();
    // Just check it returns a boolean without crashing
    EXPECT_TRUE(available || !available);
}

TEST(FaissUtilsTest, IsFaissGpuAvailableReturnsBoolean) {
    bool available = IsFaissGpuAvailable();
    EXPECT_TRUE(available || !available);
}

TEST(FaissUtilsTest, GetFaissVersionReturnsString) {
    std::string version = GetFaissVersion();
    EXPECT_FALSE(version.empty());
}

TEST(FaissUtilsTest, ComputeRecallAtK_PerfectMatch) {
    std::vector<int64_t> approx = {0, 1, 2, 3, 4};
    std::vector<int64_t> exact = {0, 1, 2, 3, 4};
    float recall = ComputeRecallAtK(approx, exact, 1, 5);
    EXPECT_FLOAT_EQ(recall, 1.0f);
}

TEST(FaissUtilsTest, ComputeRecallAtK_NoMatch) {
    std::vector<int64_t> approx = {0, 1, 2, 3, 4};
    std::vector<int64_t> exact = {5, 6, 7, 8, 9};
    float recall = ComputeRecallAtK(approx, exact, 1, 5);
    EXPECT_FLOAT_EQ(recall, 0.0f);
}

TEST(FaissUtilsTest, ComputeRecallAtK_PartialMatch) {
    std::vector<int64_t> approx = {0, 1, 2, 5, 6};
    std::vector<int64_t> exact = {0, 1, 2, 3, 4};
    float recall = ComputeRecallAtK(approx, exact, 1, 5);
    EXPECT_FLOAT_EQ(recall, 0.6f);  // 3 out of 5
}

TEST(FaissUtilsTest, ComputeRecallAtK_MultipleQueries) {
    // Two queries, k=3
    std::vector<int64_t> approx = {0, 1, 2, 3, 4, 5};
    std::vector<int64_t> exact = {0, 1, 2, 3, 4, 5};
    float recall = ComputeRecallAtK(approx, exact, 2, 3);
    EXPECT_FLOAT_EQ(recall, 1.0f);
}

TEST(FaissUtilsTest, ComputeRecallAtK_EmptyInput) {
    std::vector<int64_t> approx = {};
    std::vector<int64_t> exact = {};
    float recall = ComputeRecallAtK(approx, exact, 0, 5);
    EXPECT_FLOAT_EQ(recall, 0.0f);  // Edge case: 0/0 handled as 0
}

// ========== Index Creation Tests ==========

TEST_F(FaissWrapperTest, CreateFlatL2Index) {
    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);

    if (IsFaissAvailable()) {
        ASSERT_NE(index, nullptr);
        EXPECT_EQ(index->Config().index_type, FaissIndexType::FlatL2);
        EXPECT_EQ(index->EmbeddingDim(), config_.embedding_dim);
        EXPECT_TRUE(index->IsTrained());  // Flat is always trained
    } else {
        EXPECT_EQ(index, nullptr);
    }
}

TEST_F(FaissWrapperTest, CreateFlatIPIndex) {
    config_.index_type = FaissIndexType::FlatIP;
    auto index = FaissIndex::Create(config_);

    if (IsFaissAvailable()) {
        ASSERT_NE(index, nullptr);
        EXPECT_EQ(index->Config().index_type, FaissIndexType::FlatIP);
        EXPECT_TRUE(index->IsTrained());
    } else {
        EXPECT_EQ(index, nullptr);
    }
}

TEST_F(FaissWrapperTest, CreateIVFFlatIndex) {
    config_.index_type = FaissIndexType::IVFFlat;
    auto index = FaissIndex::Create(config_);

    if (IsFaissAvailable()) {
        ASSERT_NE(index, nullptr);
        EXPECT_EQ(index->Config().index_type, FaissIndexType::IVFFlat);
        EXPECT_FALSE(index->IsTrained());  // IVF needs training
    } else {
        EXPECT_EQ(index, nullptr);
    }
}

TEST_F(FaissWrapperTest, CreateIVFPQIndex) {
    config_.index_type = FaissIndexType::IVFPQ;
    config_.m_subquantizers = 8;
    config_.nbits_per_idx = 8;
    auto index = FaissIndex::Create(config_);

    if (IsFaissAvailable()) {
        ASSERT_NE(index, nullptr);
        EXPECT_EQ(index->Config().index_type, FaissIndexType::IVFPQ);
        EXPECT_FALSE(index->IsTrained());
    } else {
        EXPECT_EQ(index, nullptr);
    }
}

// ========== Index Operations Tests (require FAISS) ==========

TEST_F(FaissWrapperTest, FlatL2_AddAndSearch) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 100;
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    // Add vectors
    size_t added = index->Add(vectors, num_vectors);
    EXPECT_EQ(added, num_vectors);
    EXPECT_EQ(index->NumVectors(), num_vectors);

    // Search for first vector (should find itself)
    auto result = index->Search(std::span<const float>(vectors.data(), dim), 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->indices.size(), 1);
    EXPECT_EQ(result->indices[0], 0);  // First vector should match itself
    EXPECT_NEAR(result->distances[0], 0.0f, 1e-5f);  // Distance to self is ~0
}

TEST_F(FaissWrapperTest, FlatL2_BatchSearch) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 100;
    const size_t dim = config_.embedding_dim;
    const size_t k = 5;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    index->Add(vectors, num_vectors);

    // Search for first 10 vectors
    const size_t num_queries = 10;
    auto result = index->SearchBatch(
        std::span<const float>(vectors.data(), num_queries * dim),
        num_queries,
        k);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->num_queries, num_queries);
    EXPECT_EQ(result->k, k);
    EXPECT_EQ(result->indices.size(), num_queries * k);
    EXPECT_EQ(result->distances.size(), num_queries * k);

    // Each query should find itself as the nearest neighbor
    for (size_t q = 0; q < num_queries; ++q) {
        EXPECT_EQ(result->indices[q * k], static_cast<int64_t>(q));
        EXPECT_NEAR(result->distances[q * k], 0.0f, 1e-5f);
    }
}

TEST_F(FaissWrapperTest, FlatL2_SearchBatchInto) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 50;
    const size_t dim = config_.embedding_dim;
    const size_t k = 3;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    index->Add(vectors, num_vectors);

    const size_t num_queries = 5;
    std::vector<int64_t> indices(num_queries * k);
    std::vector<float> distances(num_queries * k);

    bool success = index->SearchBatchInto(
        std::span<const float>(vectors.data(), num_queries * dim),
        num_queries,
        k,
        indices,
        distances);

    EXPECT_TRUE(success);
    for (size_t q = 0; q < num_queries; ++q) {
        EXPECT_EQ(indices[q * k], static_cast<int64_t>(q));
    }
}

TEST_F(FaissWrapperTest, IVFFlat_TrainAddAndSearch) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::IVFFlat;
    config_.nlist = 8;  // Few clusters for small test dataset
    config_.nprobe = 4;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);
    EXPECT_FALSE(index->IsTrained());

    const size_t num_vectors = 200;  // Need enough for training
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    // Train
    bool trained = index->Train(vectors, num_vectors);
    EXPECT_TRUE(trained);
    EXPECT_TRUE(index->IsTrained());

    // Add
    size_t added = index->Add(vectors, num_vectors);
    EXPECT_EQ(added, num_vectors);

    // Search
    const size_t k = 5;
    auto result = index->Search(std::span<const float>(vectors.data(), dim), k);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->indices.size(), k);

    // First result should be query itself (index 0)
    EXPECT_EQ(result->indices[0], 0);
}

TEST_F(FaissWrapperTest, IVFFlat_TrainAndAdd) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::IVFFlat;
    config_.nlist = 4;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 100;
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    // TrainAndAdd in one call
    bool success = index->TrainAndAdd(vectors, num_vectors);
    EXPECT_TRUE(success);
    EXPECT_TRUE(index->IsTrained());
    EXPECT_EQ(index->NumVectors(), num_vectors);
}

TEST_F(FaissWrapperTest, IVFFlat_AddBeforeTrainFails) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::IVFFlat;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 100;
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    // Should fail because not trained
    size_t added = index->Add(vectors, num_vectors);
    EXPECT_EQ(added, 0);
    EXPECT_FALSE(index->LastError().empty());
}

TEST_F(FaissWrapperTest, SetAndGetNprobe) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::IVFFlat;
    config_.nprobe = 8;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    // Train first
    const size_t num_vectors = 100;
    auto vectors = GenerateRandomVectors(num_vectors, config_.embedding_dim);
    index->Train(vectors, num_vectors);

    // Default nprobe from config
    // Note: nprobe is typically set after training on some implementations

    // Set new nprobe
    index->SetNprobe(16);
    EXPECT_EQ(index->GetNprobe(), 16);

    index->SetNprobe(1);
    EXPECT_EQ(index->GetNprobe(), 1);
}

TEST_F(FaissWrapperTest, GetStats) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    auto stats = index->GetStats();
    EXPECT_EQ(stats.num_vectors, 0);
    EXPECT_EQ(stats.embedding_dim, config_.embedding_dim);
    EXPECT_TRUE(stats.is_trained);

    // Add vectors
    const size_t num_vectors = 50;
    auto vectors = GenerateRandomVectors(num_vectors, config_.embedding_dim);
    index->Add(vectors, num_vectors);

    stats = index->GetStats();
    EXPECT_EQ(stats.num_vectors, num_vectors);
    EXPECT_GT(stats.add_time_ms, 0.0f);
}

// ========== Accuracy Tests ==========

TEST_F(FaissWrapperTest, FlatL2_ExactKnnMatch) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    config_.normalize_vectors = false;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 200;
    const size_t dim = config_.embedding_dim;
    const size_t k = 10;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    index->Add(vectors, num_vectors);

    // Test a few queries against brute-force
    for (size_t q = 0; q < 5; ++q) {
        auto faiss_result = index->Search(
            std::span<const float>(vectors.data() + q * dim, dim), k);
        auto brute_force = BruteForceKnn(vectors, vectors.data() + q * dim, num_vectors, dim, k);

        ASSERT_TRUE(faiss_result.has_value());

        // Flat index should give exact results
        for (size_t i = 0; i < k; ++i) {
            EXPECT_EQ(faiss_result->indices[i], brute_force[i])
                << "Query " << q << ", position " << i;
        }
    }
}

TEST_F(FaissWrapperTest, IVFFlat_HighRecall) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::IVFFlat;
    config_.nlist = 16;
    config_.nprobe = 8;  // High nprobe for good recall
    config_.normalize_vectors = false;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 500;
    const size_t dim = config_.embedding_dim;
    const size_t k = 10;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    index->TrainAndAdd(vectors, num_vectors);

    // Compute recall for several queries
    const size_t num_test_queries = 20;
    size_t total_found = 0;

    for (size_t q = 0; q < num_test_queries; ++q) {
        auto faiss_result = index->Search(
            std::span<const float>(vectors.data() + q * dim, dim), k);
        auto brute_force = BruteForceKnn(vectors, vectors.data() + q * dim, num_vectors, dim, k);

        ASSERT_TRUE(faiss_result.has_value());

        // Count how many of the true k-NN were found
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                if (faiss_result->indices[i] == brute_force[j]) {
                    ++total_found;
                    break;
                }
            }
        }
    }

    float recall = static_cast<float>(total_found) / static_cast<float>(num_test_queries * k);
    EXPECT_GT(recall, 0.8f) << "IVFFlat recall should be > 80% with nprobe=8";
}

// ========== Edge Cases ==========

TEST_F(FaissWrapperTest, EmptyIndexSearch) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    // Search on empty index
    auto vectors = GenerateRandomVectors(1, config_.embedding_dim);
    auto result = index->Search(vectors, 5);

    // Should return nullopt or empty results
    EXPECT_FALSE(result.has_value());
}

TEST_F(FaissWrapperTest, SearchKLargerThanIndex) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 5;
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    index->Add(vectors, num_vectors);

    // Search for k=10 when only 5 vectors exist
    auto result = index->Search(std::span<const float>(vectors.data(), dim), 10);
    ASSERT_TRUE(result.has_value());

    // Should return 5 results (all we have), rest may be -1 or similar
    size_t valid_count = 0;
    for (auto idx : result->indices) {
        if (idx >= 0 && idx < static_cast<int64_t>(num_vectors)) {
            ++valid_count;
        }
    }
    EXPECT_GE(valid_count, num_vectors);
}

TEST_F(FaissWrapperTest, AddWithIds) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);

    const size_t num_vectors = 50;
    const size_t dim = config_.embedding_dim;
    auto vectors = GenerateRandomVectors(num_vectors, dim);

    // Custom IDs starting at 1000
    std::vector<int64_t> ids(num_vectors);
    for (size_t i = 0; i < num_vectors; ++i) {
        ids[i] = static_cast<int64_t>(1000 + i);
    }

    size_t added = index->AddWithIds(vectors, ids, num_vectors);
    EXPECT_EQ(added, num_vectors);

    // Search should return custom IDs
    auto result = index->Search(std::span<const float>(vectors.data(), dim), 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->indices[0], 1000);  // First vector has ID 1000
}

TEST_F(FaissWrapperTest, IsReadyState) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    // Flat index is ready immediately
    config_.index_type = FaissIndexType::FlatL2;
    auto flat_index = FaissIndex::Create(config_);
    ASSERT_NE(flat_index, nullptr);
    EXPECT_TRUE(flat_index->IsReady());

    // IVF index is not ready until trained
    config_.index_type = FaissIndexType::IVFFlat;
    auto ivf_index = FaissIndex::Create(config_);
    ASSERT_NE(ivf_index, nullptr);
    EXPECT_FALSE(ivf_index->IsReady());

    // After training, should be ready
    auto vectors = GenerateRandomVectors(100, config_.embedding_dim);
    ivf_index->Train(vectors, 100);
    EXPECT_TRUE(ivf_index->IsReady());
}

// ========== Provider Tests ==========

TEST_F(FaissWrapperTest, GpuProviderFallback) {
    // Request GPU but should work (fall back to CPU if GPU unavailable)
    config_.index_type = FaissIndexType::FlatL2;
    config_.provider = FaissProvider::GPU;

    auto index = FaissIndex::Create(config_);

    if (IsFaissAvailable()) {
        ASSERT_NE(index, nullptr);
        // Should either use GPU or fall back to CPU
        auto provider = index->ActiveProvider();
        EXPECT_TRUE(provider == FaissProvider::GPU || provider == FaissProvider::CPU);
    } else {
        EXPECT_EQ(index, nullptr);
    }
}

TEST_F(FaissWrapperTest, IsGpuEnabled) {
    if (!IsFaissAvailable()) {
        GTEST_SKIP() << "FAISS not available";
    }

    config_.index_type = FaissIndexType::FlatL2;
    config_.provider = FaissProvider::CPU;
    auto index = FaissIndex::Create(config_);
    ASSERT_NE(index, nullptr);
    EXPECT_FALSE(index->IsGpuEnabled());
}

}  // namespace
}  // namespace llmap::self_interference
