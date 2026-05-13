// Unit tests for ReferenceIndex (Stage 2 reference WaveCollapse)

#include "reference_collapse/reference_index.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <random>

namespace llmap {
namespace {

// Helper to create a small test reference
std::vector<ReferenceTarget> MakeTestTargets() {
    return {
        {.name = "chr1", .length = 10'000'000, .md5 = "abc123"},
        {.name = "chr2", .length = 5'000'000, .md5 = "def456"},
        {.name = "chrM", .length = 16'569, .md5 = "mito01"},
    };
}

// Helper to generate random embeddings
std::vector<float> MakeRandomEmbeddings(std::size_t num_vectors,
                                        std::size_t dim,
                                        unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    std::vector<float> result(num_vectors * dim);
    for (auto& v : result) {
        v = dist(rng);
    }
    return result;
}

class ReferenceIndexTest : public ::testing::Test {
protected:
    void SetUp() override {
        targets_ = MakeTestTargets();
        config_ = {
            .l1_granularity = 5'000'000,
            .l2_granularity = 50'000,
            .embedding_dim = 256,
            .include_embeddings = true,
            .reference_version = "TestRef1.0",
        };
    }

    void TearDown() override {
        // Clean up temp files
        if (std::filesystem::exists(temp_path_)) {
            std::filesystem::remove(temp_path_);
        }
    }

    std::vector<ReferenceTarget> targets_;
    ReferenceIndexConfig config_;
    std::filesystem::path temp_path_ = "/tmp/test_reference_index.llmap.idx";
};

// --- Builder tests ---

TEST_F(ReferenceIndexTest, BuilderCreatesValidIndex) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    EXPECT_EQ(index->NumTargets(), 3);
    EXPECT_EQ(index->ReferenceVersion(), "TestRef1.0");
}

TEST_F(ReferenceIndexTest, BuilderFailsWithNoTargets) {
    auto builder = ReferenceIndex::Builder(config_);
    auto index = builder.Build();

    EXPECT_EQ(index, nullptr);
    EXPECT_FALSE(builder.LastError().empty());
}

TEST_F(ReferenceIndexTest, BuilderAddsTargetsIncrementally) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTarget(targets_[0])
                     .AddTarget(targets_[1])
                     .AddTarget(targets_[2])
                     .Build();

    ASSERT_NE(index, nullptr);
    EXPECT_EQ(index->NumTargets(), 3);
}

TEST_F(ReferenceIndexTest, BuilderSetsEmbeddings) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    // After build, compute expected bucket counts and set embeddings
    auto l0_count = index->Pyramid().l0_count();
    auto l1_count = index->Pyramid().l1_count();
    auto l2_count = index->Pyramid().l2_count();

    auto l0_emb = MakeRandomEmbeddings(l0_count, config_.embedding_dim);
    auto l1_emb = MakeRandomEmbeddings(l1_count, config_.embedding_dim);
    auto l2_emb = MakeRandomEmbeddings(l2_count, config_.embedding_dim);

    auto index2 = ReferenceIndex::Builder(config_)
                      .AddTargets(targets_)
                      .SetL0Embeddings(l0_emb)
                      .SetL1Embeddings(l1_emb)
                      .SetL2Embeddings(l2_emb)
                      .Build();

    ASSERT_NE(index2, nullptr);
    EXPECT_TRUE(index2->HasEmbeddings());
    EXPECT_EQ(index2->L0Embeddings().size(), l0_count * config_.embedding_dim);
    EXPECT_EQ(index2->L1Embeddings().size(), l1_count * config_.embedding_dim);
    EXPECT_EQ(index2->L2Embeddings().size(), l2_count * config_.embedding_dim);
}

// --- Pyramid structure tests ---

TEST_F(ReferenceIndexTest, PyramidHasCorrectL0Count) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    // L0: one bucket per target
    EXPECT_EQ(index->Pyramid().l0_count(), 3);
}

TEST_F(ReferenceIndexTest, PyramidHasCorrectL1Count) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    // chr1: 10MB / 5MB = 2 buckets
    // chr2: 5MB / 5MB = 1 bucket
    // chrM: 16kb / 5MB = 1 bucket (ceil)
    EXPECT_EQ(index->Pyramid().l1_count(), 4);
}

TEST_F(ReferenceIndexTest, PyramidHasCorrectL2Count) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    // chr1: 10MB / 50kb = 200 buckets
    // chr2: 5MB / 50kb = 100 buckets
    // chrM: 16569 / 50kb = 1 bucket (ceil)
    EXPECT_EQ(index->Pyramid().l2_count(), 301);
}

TEST_F(ReferenceIndexTest, PyramidValidates) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    EXPECT_TRUE(index->Pyramid().validate());
}

// --- Target lookup tests ---

TEST_F(ReferenceIndexTest, FindTargetReturnsCorrectTarget) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    auto chr1 = index->FindTarget("chr1");
    ASSERT_TRUE(chr1.has_value());
    EXPECT_EQ(chr1->name, "chr1");
    EXPECT_EQ(chr1->length, 10'000'000);

    auto chrM = index->FindTarget("chrM");
    ASSERT_TRUE(chrM.has_value());
    EXPECT_EQ(chrM->length, 16'569);
}

TEST_F(ReferenceIndexTest, FindTargetReturnsNulloptForUnknown) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    auto unknown = index->FindTarget("chrX");
    EXPECT_FALSE(unknown.has_value());
}

// --- Bucket lookup tests ---

TEST_F(ReferenceIndexTest, FindL2BucketAtPosition) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    // Position 0 should be in first bucket of chr1
    auto bucket0 = index->FindL2Bucket("chr1", 0);
    ASSERT_TRUE(bucket0.has_value());

    auto buckets = index->Pyramid().l2_buckets();
    EXPECT_EQ(buckets[*bucket0].target_id, "chr1");
    EXPECT_EQ(buckets[*bucket0].start, 0);

    // Position 5,000,001 should be in a middle bucket
    auto bucket_mid = index->FindL2Bucket("chr1", 5'000'001);
    ASSERT_TRUE(bucket_mid.has_value());
    EXPECT_GE(buckets[*bucket_mid].start, 5'000'000);
}

TEST_F(ReferenceIndexTest, FindL2BucketReturnsNulloptForInvalidTarget) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    auto bucket = index->FindL2Bucket("chrX", 100);
    EXPECT_FALSE(bucket.has_value());
}

TEST_F(ReferenceIndexTest, FindL2BucketsInRange) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    // Range spanning 3 L2 buckets (50kb each)
    auto buckets = index->FindL2BucketsInRange("chr1", 0, 150'000);
    EXPECT_EQ(buckets.size(), 3);

    // Verify they're sorted by start
    auto l2_buckets = index->Pyramid().l2_buckets();
    for (std::size_t i = 1; i < buckets.size(); ++i) {
        EXPECT_LT(l2_buckets[buckets[i - 1]].start,
                  l2_buckets[buckets[i]].start);
    }
}

TEST_F(ReferenceIndexTest, FindL2BucketsInRangeEmpty) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    // Range past end of chrM (16569 bp)
    auto buckets = index->FindL2BucketsInRange("chrM", 100'000, 200'000);
    EXPECT_TRUE(buckets.empty());
}

// --- Embedding access tests ---

TEST_F(ReferenceIndexTest, GetL2EmbeddingReturnsCorrectSlice) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    auto l2_count = index->Pyramid().l2_count();
    auto embeddings = MakeRandomEmbeddings(l2_count, config_.embedding_dim, 123);

    auto index2 = ReferenceIndex::Builder(config_)
                      .AddTargets(targets_)
                      .SetL2Embeddings(embeddings)
                      .Build();

    ASSERT_NE(index2, nullptr);

    // Check bucket 5's embedding
    auto emb = index2->GetL2Embedding(5);
    EXPECT_EQ(emb.size(), config_.embedding_dim);

    // Verify it matches the original data
    for (std::size_t i = 0; i < config_.embedding_dim; ++i) {
        EXPECT_FLOAT_EQ(emb[i], embeddings[5 * config_.embedding_dim + i]);
    }
}

TEST_F(ReferenceIndexTest, GetEmbeddingReturnsEmptyForNoEmbeddings) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    EXPECT_FALSE(index->HasEmbeddings());

    auto emb = index->GetL2Embedding(0);
    EXPECT_TRUE(emb.empty());
}

TEST_F(ReferenceIndexTest, GetEmbeddingReturnsEmptyForInvalidBucket) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    auto l2_count = index->Pyramid().l2_count();
    auto embeddings = MakeRandomEmbeddings(l2_count, config_.embedding_dim);

    auto index2 = ReferenceIndex::Builder(config_)
                      .AddTargets(targets_)
                      .SetL2Embeddings(embeddings)
                      .Build();

    ASSERT_NE(index2, nullptr);

    // Invalid bucket ID
    auto emb = index2->GetL2Embedding(999999);
    EXPECT_TRUE(emb.empty());
}

// --- Serialization tests ---

TEST_F(ReferenceIndexTest, SaveAndLoadRoundTrip) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);

    ASSERT_TRUE(index->Save(temp_path_));
    EXPECT_TRUE(std::filesystem::exists(temp_path_));

    auto loaded = ReferenceIndex::Load(temp_path_);
    ASSERT_NE(loaded, nullptr);

    EXPECT_EQ(loaded->NumTargets(), index->NumTargets());
    EXPECT_EQ(loaded->ReferenceVersion(), index->ReferenceVersion());
    EXPECT_EQ(loaded->Pyramid().l0_count(), index->Pyramid().l0_count());
    EXPECT_EQ(loaded->Pyramid().l1_count(), index->Pyramid().l1_count());
    EXPECT_EQ(loaded->Pyramid().l2_count(), index->Pyramid().l2_count());
}

TEST_F(ReferenceIndexTest, SaveAndLoadWithEmbeddings) {
    auto index_temp = ReferenceIndex::Builder(config_)
                          .AddTargets(targets_)
                          .Build();
    ASSERT_NE(index_temp, nullptr);

    auto l0_count = index_temp->Pyramid().l0_count();
    auto l1_count = index_temp->Pyramid().l1_count();
    auto l2_count = index_temp->Pyramid().l2_count();

    auto l0_emb = MakeRandomEmbeddings(l0_count, config_.embedding_dim, 1);
    auto l1_emb = MakeRandomEmbeddings(l1_count, config_.embedding_dim, 2);
    auto l2_emb = MakeRandomEmbeddings(l2_count, config_.embedding_dim, 3);

    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .SetL0Embeddings(l0_emb)
                     .SetL1Embeddings(l1_emb)
                     .SetL2Embeddings(l2_emb)
                     .Build();

    ASSERT_NE(index, nullptr);
    ASSERT_TRUE(index->Save(temp_path_));

    auto loaded = ReferenceIndex::Load(temp_path_);
    ASSERT_NE(loaded, nullptr);
    EXPECT_TRUE(loaded->HasEmbeddings());

    // Verify embeddings match
    EXPECT_EQ(loaded->L0Embeddings().size(), index->L0Embeddings().size());
    EXPECT_EQ(loaded->L1Embeddings().size(), index->L1Embeddings().size());
    EXPECT_EQ(loaded->L2Embeddings().size(), index->L2Embeddings().size());

    for (std::size_t i = 0; i < l2_emb.size(); ++i) {
        EXPECT_FLOAT_EQ(loaded->L2Embeddings()[i], l2_emb[i]);
    }
}

TEST_F(ReferenceIndexTest, LoadReturnsNullptrForInvalidFile) {
    auto loaded = ReferenceIndex::Load("/nonexistent/path.idx");
    EXPECT_EQ(loaded, nullptr);
}

TEST_F(ReferenceIndexTest, LoadReturnsNullptrForCorruptedFile) {
    // Write garbage to file
    {
        std::ofstream out(temp_path_, std::ios::binary);
        out << "not a valid index file";
    }

    auto loaded = ReferenceIndex::Load(temp_path_);
    EXPECT_EQ(loaded, nullptr);
}

// --- Stats tests ---

TEST_F(ReferenceIndexTest, StatsArePopulated) {
    auto builder = ReferenceIndex::Builder(config_);
    builder.AddTargets(targets_);
    auto index = builder.Build();

    ASSERT_NE(index, nullptr);

    const auto& stats = index->Stats();
    EXPECT_EQ(stats.num_targets, 3);
    EXPECT_EQ(stats.total_length, 10'000'000 + 5'000'000 + 16'569);
    EXPECT_EQ(stats.l0_buckets, 3);
    EXPECT_EQ(stats.l1_buckets, 4);
    EXPECT_EQ(stats.l2_buckets, 301);
    EXPECT_GT(stats.build_time_seconds, 0.0f);
}

// --- Utility function tests ---

TEST_F(ReferenceIndexTest, ComputeTotalLength) {
    auto total = ComputeTotalLength(targets_);
    EXPECT_EQ(total, 10'000'000 + 5'000'000 + 16'569);
}

TEST_F(ReferenceIndexTest, EstimateBucketCounts) {
    std::uint64_t total = 3'000'000'000;  // ~3 GB like human genome
    auto [l1, l2] = EstimateBucketCounts(total, config_);

    EXPECT_EQ(l1, 600);   // 3GB / 5MB
    EXPECT_EQ(l2, 60000); // 3GB / 50kb
}

// --- Biology hint tests ---

TEST_F(ReferenceIndexTest, BiologyHintsPersist) {
    BiologyHint hint{
        .prior_weight = 1.5f,
        .annotation = "IGH-Constants",
        .paralog_partner_bucket = 42,
        .expected_coverage_multiplier = 2.0f,
    };

    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .SetBiologyHint(10, hint)
                     .Build();

    ASSERT_NE(index, nullptr);

    auto retrieved = index->Pyramid().get_biology_hint(10);
    ASSERT_TRUE(retrieved.has_value());
    EXPECT_FLOAT_EQ(retrieved->prior_weight, 1.5f);
    EXPECT_EQ(retrieved->annotation.value_or(""), "IGH-Constants");
    EXPECT_EQ(retrieved->paralog_partner_bucket.value_or(0), 42);
    EXPECT_FLOAT_EQ(retrieved->expected_coverage_multiplier, 2.0f);
}

// --- Config tests ---

TEST_F(ReferenceIndexTest, ConfigPersistsAfterLoad) {
    auto index = ReferenceIndex::Builder(config_)
                     .AddTargets(targets_)
                     .Build();

    ASSERT_NE(index, nullptr);
    ASSERT_TRUE(index->Save(temp_path_));

    auto loaded = ReferenceIndex::Load(temp_path_);
    ASSERT_NE(loaded, nullptr);

    const auto& loaded_config = loaded->Config();
    EXPECT_EQ(loaded_config.l1_granularity, config_.l1_granularity);
    EXPECT_EQ(loaded_config.l2_granularity, config_.l2_granularity);
    EXPECT_EQ(loaded_config.embedding_dim, config_.embedding_dim);
    EXPECT_EQ(loaded_config.reference_version, config_.reference_version);
}

// --- Format version test ---

TEST_F(ReferenceIndexTest, FormatVersionIsConstant) {
    EXPECT_EQ(ReferenceIndex::FormatVersion(), 1);
}

}  // namespace
}  // namespace llmap
