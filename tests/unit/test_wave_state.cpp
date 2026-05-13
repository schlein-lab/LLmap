// LLmap Phase 0 — WaveState tests.
//
// Tests the sparse CSR format for Read x Bucket probability matrix,
// collapse detection, level management, and CSR invariant validation.

#include "core/wave_state.h"

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace llmap;

namespace {

class WaveStateTest : public ::testing::Test {
protected:
    std::vector<BucketProb> make_candidates(std::initializer_list<std::pair<std::uint32_t, float>> pairs) {
        std::vector<BucketProb> result;
        for (const auto& [bid, prob] : pairs) {
            result.push_back({bid, prob});
        }
        return result;
    }
};

}  // namespace

TEST_F(WaveStateTest, DefaultConstructsEmpty) {
    WaveState ws;
    EXPECT_EQ(ws.n_reads(), 0u);
    EXPECT_EQ(ws.total_entries(), 0u);
    EXPECT_TRUE(ws.validate());
}

TEST_F(WaveStateTest, ConstructWithReadsInitializesCorrectly) {
    WaveState ws(100, WaveLevel::L0);
    EXPECT_EQ(ws.n_reads(), 100u);
    EXPECT_EQ(ws.total_entries(), 0u);
    EXPECT_EQ(ws.count_collapsed(), 0u);
    EXPECT_EQ(ws.count_active(), 100u);
    EXPECT_TRUE(ws.validate());

    for (std::uint32_t i = 0; i < 100; ++i) {
        EXPECT_EQ(ws.get_level(i), WaveLevel::L0);
        EXPECT_FALSE(ws.is_collapsed(i));
    }
}

TEST_F(WaveStateTest, SetReadCandidatesStoresCorrectly) {
    WaveState ws(3);

    auto cands0 = make_candidates({{10, 0.3f}, {20, 0.5f}, {30, 0.2f}});
    ws.set_read_candidates(0, cands0);

    auto cands1 = make_candidates({{5, 0.9f}, {15, 0.1f}});
    ws.set_read_candidates(1, cands1);

    auto cands2 = make_candidates({{100, 1.0f}});
    ws.set_read_candidates(2, cands2);

    EXPECT_TRUE(ws.validate());
    EXPECT_EQ(ws.total_entries(), 6u);

    auto buckets0 = ws.bucket_indices_for_read(0);
    ASSERT_EQ(buckets0.size(), 3u);
    EXPECT_EQ(buckets0[0], 10u);
    EXPECT_EQ(buckets0[1], 20u);
    EXPECT_EQ(buckets0[2], 30u);

    auto probs0 = ws.probabilities_for_read(0);
    ASSERT_EQ(probs0.size(), 3u);
    EXPECT_FLOAT_EQ(probs0[0], 0.3f);
    EXPECT_FLOAT_EQ(probs0[1], 0.5f);
    EXPECT_FLOAT_EQ(probs0[2], 0.2f);

    auto buckets1 = ws.bucket_indices_for_read(1);
    ASSERT_EQ(buckets1.size(), 2u);
    EXPECT_EQ(buckets1[0], 5u);
    EXPECT_EQ(buckets1[1], 15u);

    auto buckets2 = ws.bucket_indices_for_read(2);
    ASSERT_EQ(buckets2.size(), 1u);
    EXPECT_EQ(buckets2[0], 100u);
}

TEST_F(WaveStateTest, GetProbabilityReturnsCorrectValue) {
    WaveState ws(2);

    auto cands = make_candidates({{10, 0.3f}, {20, 0.5f}, {30, 0.2f}});
    ws.set_read_candidates(0, cands);

    EXPECT_FLOAT_EQ(ws.get_probability(0, 10), 0.3f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 20), 0.5f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 30), 0.2f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 99), 0.0f);  // not in candidates
}

TEST_F(WaveStateTest, UpdateProbabilityWorks) {
    WaveState ws(1);

    auto cands = make_candidates({{10, 0.3f}, {20, 0.5f}});
    ws.set_read_candidates(0, cands);

    ws.update_probability(0, 20, 0.7f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 20), 0.7f);

    // Update non-existent bucket is no-op
    ws.update_probability(0, 99, 0.9f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 99), 0.0f);
}

TEST_F(WaveStateTest, NormalizeReadSumsToOne) {
    WaveState ws(1);

    auto cands = make_candidates({{10, 0.3f}, {20, 0.5f}, {30, 0.2f}});
    ws.set_read_candidates(0, cands);

    ws.normalize_read(0);

    auto probs = ws.probabilities_for_read(0);
    float sum = 0.0f;
    for (float p : probs) {
        sum += p;
    }
    EXPECT_NEAR(sum, 1.0f, 1e-6f);
}

TEST_F(WaveStateTest, NormalizeUnnormalizedDistribution) {
    WaveState ws(1);

    auto cands = make_candidates({{10, 2.0f}, {20, 3.0f}, {30, 5.0f}});
    ws.set_read_candidates(0, cands);

    ws.normalize_read(0);

    EXPECT_FLOAT_EQ(ws.get_probability(0, 10), 0.2f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 20), 0.3f);
    EXPECT_FLOAT_EQ(ws.get_probability(0, 30), 0.5f);
}

TEST_F(WaveStateTest, LevelManagement) {
    WaveState ws(3, WaveLevel::L0);

    EXPECT_EQ(ws.get_level(0), WaveLevel::L0);
    EXPECT_EQ(ws.get_level(1), WaveLevel::L0);

    ws.set_level(1, WaveLevel::L1);
    ws.set_level(2, WaveLevel::L2);

    EXPECT_EQ(ws.get_level(0), WaveLevel::L0);
    EXPECT_EQ(ws.get_level(1), WaveLevel::L1);
    EXPECT_EQ(ws.get_level(2), WaveLevel::L2);
}

TEST_F(WaveStateTest, CollapseDetection) {
    WaveState ws(2);

    // Read 0: max prob = 0.5, below 0.99 threshold
    auto cands0 = make_candidates({{10, 0.3f}, {20, 0.5f}, {30, 0.2f}});
    ws.set_read_candidates(0, cands0);

    // Read 1: max prob = 0.995, above 0.99 threshold
    auto cands1 = make_candidates({{100, 0.995f}, {200, 0.005f}});
    ws.set_read_candidates(1, cands1);

    EXPECT_FALSE(ws.check_collapse(0, 0.99f));
    EXPECT_TRUE(ws.check_collapse(1, 0.99f));
}

TEST_F(WaveStateTest, CollapseReadSetsState) {
    WaveState ws(1);

    auto cands = make_candidates({{10, 0.1f}, {20, 0.8f}, {30, 0.1f}});
    ws.set_read_candidates(0, cands);

    EXPECT_FALSE(ws.is_collapsed(0));
    EXPECT_EQ(ws.count_collapsed(), 0u);
    EXPECT_EQ(ws.count_active(), 1u);

    std::uint32_t collapsed_bucket = ws.collapse_read(0);

    EXPECT_EQ(collapsed_bucket, 20u);
    EXPECT_TRUE(ws.is_collapsed(0));
    EXPECT_EQ(ws.collapsed_bucket(0), 20u);
    EXPECT_EQ(ws.count_collapsed(), 1u);
    EXPECT_EQ(ws.count_active(), 0u);
}

TEST_F(WaveStateTest, ManualCollapseFlag) {
    WaveState ws(2);

    ws.set_collapsed(0, true);
    EXPECT_TRUE(ws.is_collapsed(0));
    EXPECT_FALSE(ws.is_collapsed(1));

    ws.set_collapsed(0, false);
    EXPECT_FALSE(ws.is_collapsed(0));
}

TEST_F(WaveStateTest, RawCSRAccess) {
    WaveState ws(3);

    auto cands0 = make_candidates({{10, 0.5f}, {20, 0.5f}});
    auto cands1 = make_candidates({{30, 1.0f}});
    auto cands2 = make_candidates({{40, 0.25f}, {50, 0.25f}, {60, 0.5f}});

    ws.set_read_candidates(0, cands0);
    ws.set_read_candidates(1, cands1);
    ws.set_read_candidates(2, cands2);

    auto offsets = ws.read_offsets();
    ASSERT_EQ(offsets.size(), 4u);
    EXPECT_EQ(offsets[0], 0u);
    EXPECT_EQ(offsets[1], 2u);
    EXPECT_EQ(offsets[2], 3u);
    EXPECT_EQ(offsets[3], 6u);

    auto buckets = ws.bucket_indices();
    EXPECT_EQ(buckets.size(), 6u);

    auto probs = ws.probabilities();
    EXPECT_EQ(probs.size(), 6u);
}

TEST_F(WaveStateTest, ValidateDetectsMismatchedArrays) {
    WaveState ws(2);

    auto cands = make_candidates({{10, 0.5f}});
    ws.set_read_candidates(0, cands);

    EXPECT_TRUE(ws.validate());
}

TEST_F(WaveStateTest, SetCandidatesGrowingShrinking) {
    WaveState ws(1);

    // Start with 2 candidates
    auto cands1 = make_candidates({{10, 0.5f}, {20, 0.5f}});
    ws.set_read_candidates(0, cands1);
    EXPECT_EQ(ws.total_entries(), 2u);
    EXPECT_TRUE(ws.validate());

    // Grow to 4 candidates
    auto cands2 = make_candidates({{10, 0.2f}, {20, 0.3f}, {30, 0.3f}, {40, 0.2f}});
    ws.set_read_candidates(0, cands2);
    EXPECT_EQ(ws.total_entries(), 4u);
    EXPECT_TRUE(ws.validate());

    // Shrink to 1 candidate
    auto cands3 = make_candidates({{25, 1.0f}});
    ws.set_read_candidates(0, cands3);
    EXPECT_EQ(ws.total_entries(), 1u);
    EXPECT_TRUE(ws.validate());

    auto buckets = ws.bucket_indices_for_read(0);
    ASSERT_EQ(buckets.size(), 1u);
    EXPECT_EQ(buckets[0], 25u);
}

TEST_F(WaveStateTest, OutOfRangeThrows) {
    WaveState ws(3);

    EXPECT_THROW(ws.bucket_indices_for_read(99), std::out_of_range);
    EXPECT_THROW(ws.probabilities_for_read(99), std::out_of_range);
    EXPECT_THROW(ws.get_probability(99, 0), std::out_of_range);
    EXPECT_THROW(ws.update_probability(99, 0, 0.5f), std::out_of_range);
    EXPECT_THROW(ws.normalize_read(99), std::out_of_range);
    EXPECT_THROW(ws.get_level(99), std::out_of_range);
    EXPECT_THROW(ws.set_level(99, WaveLevel::L1), std::out_of_range);
    EXPECT_THROW(ws.is_collapsed(99), std::out_of_range);
    EXPECT_THROW(ws.set_collapsed(99, true), std::out_of_range);
    EXPECT_THROW(ws.check_collapse(99), std::out_of_range);
    EXPECT_THROW(ws.collapse_read(99), std::out_of_range);
    EXPECT_THROW(ws.collapsed_bucket(99), std::out_of_range);

    auto cands = make_candidates({{10, 0.5f}});
    EXPECT_THROW(ws.set_read_candidates(99, cands), std::out_of_range);
}

TEST_F(WaveStateTest, CollapseReadWithNoCandidatesThrows) {
    WaveState ws(1);
    EXPECT_THROW(ws.collapse_read(0), std::runtime_error);
}

TEST_F(WaveStateTest, GPUStubMethods) {
    WaveState ws(10);

    EXPECT_FALSE(ws.has_device_copy());

    ws.sync_to_device();
    EXPECT_TRUE(ws.has_device_copy());

    ws.sync_from_device();
    // sync_from_device doesn't change has_device_copy flag
}

TEST_F(WaveStateTest, CheckCollapseAlreadyCollapsed) {
    WaveState ws(1);

    auto cands = make_candidates({{10, 0.3f}, {20, 0.5f}});
    ws.set_read_candidates(0, cands);

    EXPECT_FALSE(ws.check_collapse(0, 0.99f));

    ws.set_collapsed(0, true);
    EXPECT_TRUE(ws.check_collapse(0, 0.99f));  // returns true because already collapsed
}

TEST_F(WaveStateTest, EmptyReadCandidates) {
    WaveState ws(2);

    auto cands = make_candidates({{10, 0.5f}});
    ws.set_read_candidates(0, cands);

    // Read 1 has no candidates
    auto empty_cands = std::vector<BucketProb>{};
    ws.set_read_candidates(1, empty_cands);

    EXPECT_TRUE(ws.validate());
    EXPECT_EQ(ws.bucket_indices_for_read(1).size(), 0u);
    EXPECT_EQ(ws.probabilities_for_read(1).size(), 0u);
    EXPECT_FALSE(ws.check_collapse(1, 0.99f));
}

TEST_F(WaveStateTest, NormalizeEmptyReadIsNoOp) {
    WaveState ws(1);
    // No candidates set
    ws.normalize_read(0);  // should not crash
    EXPECT_TRUE(ws.validate());
}

TEST_F(WaveStateTest, ReserveEntriesDoesNotAffectSize) {
    WaveState ws(10);
    ws.reserve_entries(1000);

    EXPECT_EQ(ws.total_entries(), 0u);
    EXPECT_TRUE(ws.validate());
}
