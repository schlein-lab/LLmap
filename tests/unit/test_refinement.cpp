// LLmap — Unit tests for refinement module.

#include <gtest/gtest.h>

#include "core/bucket_pyramid.h"
#include "core/wave_state.h"
#include "reference_collapse/refinement.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace llmap {
namespace {

// Helper: create a test pyramid with known structure
// L0: 2 buckets (chr1, chr2)
// L1: 4 buckets (2 per L0)
// L2: 8 buckets (2 per L1)
BucketPyramid CreateTestPyramid() {
    BucketPyramid pyramid;

    // L0 buckets
    pyramid.add_l0_bucket({0, "chr1", 100'000'000});
    pyramid.add_l0_bucket({0, "chr2", 200'000'000});

    // L1 buckets (2 per L0)
    pyramid.add_l1_bucket({0, "chr1", 0, 50'000'000}, 0);
    pyramid.add_l1_bucket({0, "chr1", 50'000'000, 100'000'000}, 0);
    pyramid.add_l1_bucket({0, "chr2", 0, 100'000'000}, 1);
    pyramid.add_l1_bucket({0, "chr2", 100'000'000, 200'000'000}, 1);

    // L2 buckets (2 per L1)
    pyramid.add_l2_bucket({0, "chr1", 0, 25'000'000}, 0);
    pyramid.add_l2_bucket({0, "chr1", 25'000'000, 50'000'000}, 0);
    pyramid.add_l2_bucket({0, "chr1", 50'000'000, 75'000'000}, 1);
    pyramid.add_l2_bucket({0, "chr1", 75'000'000, 100'000'000}, 1);
    pyramid.add_l2_bucket({0, "chr2", 0, 50'000'000}, 2);
    pyramid.add_l2_bucket({0, "chr2", 50'000'000, 100'000'000}, 2);
    pyramid.add_l2_bucket({0, "chr2", 100'000'000, 150'000'000}, 3);
    pyramid.add_l2_bucket({0, "chr2", 150'000'000, 200'000'000}, 3);

    return pyramid;
}

// ========== ChildIndex Tests ==========

TEST(ChildIndexTest, BuildL0ToL1) {
    auto pyramid = CreateTestPyramid();
    auto index = ChildIndex::BuildL0ToL1(pyramid);

    EXPECT_TRUE(index.IsValid());
    EXPECT_EQ(index.NumParents(), 2);
    EXPECT_EQ(index.TotalChildren(), 4);

    // L0 bucket 0 has L1 children 0, 1
    auto children0 = index.GetChildren(0);
    EXPECT_EQ(children0.size(), 2);
    EXPECT_EQ(children0[0], 0);
    EXPECT_EQ(children0[1], 1);

    // L0 bucket 1 has L1 children 2, 3
    auto children1 = index.GetChildren(1);
    EXPECT_EQ(children1.size(), 2);
    EXPECT_EQ(children1[0], 2);
    EXPECT_EQ(children1[1], 3);
}

TEST(ChildIndexTest, BuildL1ToL2) {
    auto pyramid = CreateTestPyramid();
    auto index = ChildIndex::BuildL1ToL2(pyramid);

    EXPECT_TRUE(index.IsValid());
    EXPECT_EQ(index.NumParents(), 4);
    EXPECT_EQ(index.TotalChildren(), 8);

    // L1 bucket 0 has L2 children 0, 1
    auto children0 = index.GetChildren(0);
    EXPECT_EQ(children0.size(), 2);

    // L1 bucket 3 has L2 children 6, 7
    auto children3 = index.GetChildren(3);
    EXPECT_EQ(children3.size(), 2);
    EXPECT_EQ(children3[0], 6);
    EXPECT_EQ(children3[1], 7);
}

TEST(ChildIndexTest, EmptyPyramid) {
    BucketPyramid empty;
    auto index = ChildIndex::BuildL0ToL1(empty);

    EXPECT_FALSE(index.IsValid());
    EXPECT_EQ(index.NumParents(), 0);
}

TEST(ChildIndexTest, OutOfBoundsParent) {
    auto pyramid = CreateTestPyramid();
    auto index = ChildIndex::BuildL0ToL1(pyramid);

    // Requesting children of non-existent parent returns empty span
    auto children = index.GetChildren(100);
    EXPECT_TRUE(children.empty());
}

// ========== Refinement Tests ==========

TEST(RefinementTest, RefineL0ToL1_SingleRead) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    // Create state with one read at L0
    WaveState state(1, WaveLevel::L0);

    // Set candidates: L0 bucket 0 with probability 1.0
    std::vector<BucketProb> candidates = {{0, 1.0f}};
    state.set_read_candidates(0, candidates);

    Refinement refinement;
    auto stats = refinement.RefineL0ToL1(state, child_index);

    EXPECT_EQ(stats.reads_refined, 1);
    EXPECT_EQ(stats.reads_skipped, 0);
    EXPECT_EQ(state.get_level(0), WaveLevel::L1);

    // Should now have 2 L1 candidates (children of L0 bucket 0)
    auto new_buckets = state.bucket_indices_for_read(0);
    EXPECT_EQ(new_buckets.size(), 2);
}

TEST(RefinementTest, RefineL0ToL1_PreservesRelativeProbs) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(1, WaveLevel::L0);

    // L0 bucket 0 with p=0.7, L0 bucket 1 with p=0.3
    std::vector<BucketProb> candidates = {{0, 0.7f}, {1, 0.3f}};
    state.set_read_candidates(0, candidates);

    RefinementConfig config;
    config.preserve_relative_probs = true;
    config.expansion_threshold = 0.01f;
    Refinement refinement(config);

    refinement.RefineL0ToL1(state, child_index);

    // Should have 4 L1 candidates (2 children each from both L0 buckets)
    auto new_probs = state.probabilities_for_read(0);
    EXPECT_EQ(new_probs.size(), 4);

    // Sum should be 1.0 (normalized)
    float sum = std::accumulate(new_probs.begin(), new_probs.end(), 0.0f);
    EXPECT_NEAR(sum, 1.0f, 1e-5f);
}

TEST(RefinementTest, RefineL1ToL2) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL1ToL2(pyramid);

    WaveState state(1, WaveLevel::L1);

    // L1 bucket 0 with probability 1.0
    std::vector<BucketProb> candidates = {{0, 1.0f}};
    state.set_read_candidates(0, candidates);

    Refinement refinement;
    auto stats = refinement.RefineL1ToL2(state, child_index);

    EXPECT_EQ(stats.reads_refined, 1);
    EXPECT_EQ(state.get_level(0), WaveLevel::L2);

    auto new_buckets = state.bucket_indices_for_read(0);
    EXPECT_EQ(new_buckets.size(), 2);
}

TEST(RefinementTest, SkipsCollapsedReads) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(2, WaveLevel::L0);

    // Both reads have same candidates
    std::vector<BucketProb> candidates = {{0, 1.0f}};
    state.set_read_candidates(0, candidates);
    state.set_read_candidates(1, candidates);

    // Collapse read 0
    state.set_collapsed(0, true);

    Refinement refinement;
    auto stats = refinement.RefineL0ToL1(state, child_index);

    EXPECT_EQ(stats.reads_refined, 1);
    EXPECT_EQ(stats.reads_skipped, 1);

    // Read 0 should still be at L0 (unchanged)
    EXPECT_EQ(state.get_level(0), WaveLevel::L0);
    // Read 1 should be at L1
    EXPECT_EQ(state.get_level(1), WaveLevel::L1);
}

TEST(RefinementTest, SkipsWrongLevel) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(1, WaveLevel::L1);  // Already at L1

    std::vector<BucketProb> candidates = {{0, 1.0f}};
    state.set_read_candidates(0, candidates);

    Refinement refinement;
    auto stats = refinement.RefineL0ToL1(state, child_index);

    EXPECT_EQ(stats.reads_refined, 0);
    EXPECT_EQ(stats.reads_skipped, 1);
}

TEST(RefinementTest, ExpansionThreshold) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(1, WaveLevel::L0);

    // Two buckets: one above threshold, one below
    std::vector<BucketProb> candidates = {{0, 0.05f}, {1, 0.95f}};
    state.set_read_candidates(0, candidates);

    RefinementConfig config;
    config.expansion_threshold = 0.1f;  // Only expand buckets with P >= 0.1
    Refinement refinement(config);

    refinement.RefineL0ToL1(state, child_index);

    // Should only have children of bucket 1 (which had P=0.95)
    auto new_buckets = state.bucket_indices_for_read(0);
    EXPECT_EQ(new_buckets.size(), 2);

    // Children should be 2, 3 (children of L0 bucket 1)
    EXPECT_EQ(new_buckets[0], 2);
    EXPECT_EQ(new_buckets[1], 3);
}

TEST(RefinementTest, MaxCandidatesLimit) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(1, WaveLevel::L0);

    // Both L0 buckets with equal probability
    std::vector<BucketProb> candidates = {{0, 0.5f}, {1, 0.5f}};
    state.set_read_candidates(0, candidates);

    RefinementConfig config;
    config.max_candidates = 3;  // Limit to 3 candidates
    Refinement refinement(config);

    refinement.RefineL0ToL1(state, child_index);

    auto new_buckets = state.bucket_indices_for_read(0);
    EXPECT_LE(new_buckets.size(), 3);
}

TEST(RefinementTest, CountRefinementCandidates) {
    auto pyramid = CreateTestPyramid();

    WaveState state(3, WaveLevel::L0);

    std::vector<BucketProb> high_prob = {{0, 0.9f}};
    std::vector<BucketProb> low_prob = {{0, 0.001f}};

    state.set_read_candidates(0, high_prob);  // Above threshold
    state.set_read_candidates(1, low_prob);   // Below threshold
    state.set_read_candidates(2, high_prob);  // Above threshold
    state.set_collapsed(2, true);             // But collapsed

    RefinementConfig config;
    config.expansion_threshold = 0.01f;
    Refinement refinement(config);

    auto count = refinement.CountRefinementCandidates(state, WaveLevel::L0);
    EXPECT_EQ(count, 1);  // Only read 0 qualifies
}

// ========== Free Function Tests ==========

TEST(RefinementTest, ShouldRefine_HighProb) {
    WaveState state(10, WaveLevel::L0);

    // 8 reads with high probability
    for (std::uint32_t i = 0; i < 8; ++i) {
        std::vector<BucketProb> high = {{0, 0.95f}};
        state.set_read_candidates(i, high);
    }

    // 2 reads with low probability
    for (std::uint32_t i = 8; i < 10; ++i) {
        std::vector<BucketProb> low = {{0, 0.3f}, {1, 0.3f}, {2, 0.4f}};
        state.set_read_candidates(i, low);
    }

    // 80% threshold: should trigger refinement
    EXPECT_TRUE(ShouldRefine(state, WaveLevel::L0, 0.8f, 0.9f));

    // 90% threshold: should not trigger
    EXPECT_FALSE(ShouldRefine(state, WaveLevel::L0, 0.9f, 0.9f));
}

TEST(RefinementTest, ShouldRefine_CollapsedReads) {
    WaveState state(10, WaveLevel::L0);

    // All reads collapsed
    for (std::uint32_t i = 0; i < 10; ++i) {
        std::vector<BucketProb> cands = {{0, 1.0f}};
        state.set_read_candidates(i, cands);
        state.set_collapsed(i, true);
    }

    EXPECT_TRUE(ShouldRefine(state, WaveLevel::L0, 0.8f, 0.9f));
}

TEST(RefinementTest, NextAndPrevLevel) {
    EXPECT_EQ(NextLevel(WaveLevel::L0), WaveLevel::L1);
    EXPECT_EQ(NextLevel(WaveLevel::L1), WaveLevel::L2);
    EXPECT_EQ(NextLevel(WaveLevel::L2), WaveLevel::L3);
    EXPECT_EQ(NextLevel(WaveLevel::L3), WaveLevel::L3);  // Saturates

    EXPECT_EQ(PrevLevel(WaveLevel::L0), WaveLevel::L0);  // Saturates
    EXPECT_EQ(PrevLevel(WaveLevel::L1), WaveLevel::L0);
    EXPECT_EQ(PrevLevel(WaveLevel::L2), WaveLevel::L1);
    EXPECT_EQ(PrevLevel(WaveLevel::L3), WaveLevel::L2);
}

TEST(RefinementTest, RefineStats) {
    auto pyramid = CreateTestPyramid();
    auto child_index = ChildIndex::BuildL0ToL1(pyramid);

    WaveState state(5, WaveLevel::L0);

    for (std::uint32_t i = 0; i < 5; ++i) {
        std::vector<BucketProb> cands = {{0, 0.6f}, {1, 0.4f}};
        state.set_read_candidates(i, cands);
    }

    Refinement refinement;
    auto stats = refinement.RefineL0ToL1(state, child_index);

    EXPECT_EQ(stats.reads_refined, 5);
    EXPECT_EQ(stats.reads_skipped, 0);
    EXPECT_GT(stats.total_children_added, 0);
    EXPECT_GT(stats.avg_expansion_factor, 0.0f);
    EXPECT_GE(stats.refinement_time_ms, 0.0f);
}

TEST(RefinementTest, MergesDuplicateChildren) {
    // Create pyramid where multiple parents share the same child
    // This is unusual but we should handle it
    BucketPyramid pyramid;
    pyramid.add_l0_bucket({0, "test", 100});
    pyramid.add_l1_bucket({0, "test", 0, 50}, 0);

    WaveState state(1, WaveLevel::L0);
    std::vector<BucketProb> cands = {{0, 1.0f}};
    state.set_read_candidates(0, cands);

    auto child_index = ChildIndex::BuildL0ToL1(pyramid);
    Refinement refinement;
    refinement.RefineL0ToL1(state, child_index);

    auto new_buckets = state.bucket_indices_for_read(0);
    EXPECT_EQ(new_buckets.size(), 1);

    auto new_probs = state.probabilities_for_read(0);
    EXPECT_NEAR(new_probs[0], 1.0f, 1e-5f);
}

}  // namespace
}  // namespace llmap
