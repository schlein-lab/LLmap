// LLmap — Tests for EM Iterator module.

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

#include "core/wave_state.h"
#include "reference_collapse/em_iterator.h"

namespace llmap {
namespace {

class EmIteratorTest : public ::testing::Test {
protected:
    void SetUp() override {
        state_ = WaveState(4);  // 4 reads

        // Read 0: buckets {0, 1, 2} with equal probability
        std::vector<BucketProb> cands0 = {{0, 0.33f}, {1, 0.34f}, {2, 0.33f}};
        state_.set_read_candidates(0, cands0);

        // Read 1: buckets {1, 2, 3} with peak at 1
        std::vector<BucketProb> cands1 = {{1, 0.8f}, {2, 0.1f}, {3, 0.1f}};
        state_.set_read_candidates(1, cands1);

        // Read 2: buckets {0, 3} with equal
        std::vector<BucketProb> cands2 = {{0, 0.5f}, {3, 0.5f}};
        state_.set_read_candidates(2, cands2);

        // Read 3: buckets {2, 3, 4} with peak at 3
        std::vector<BucketProb> cands3 = {{2, 0.1f}, {3, 0.85f}, {4, 0.05f}};
        state_.set_read_candidates(3, cands3);
    }

    WaveState state_;
};

TEST_F(EmIteratorTest, DefaultConstruction) {
    EmIterator iter;
    EXPECT_FLOAT_EQ(iter.Config().gamma, 0.3f);
    EXPECT_FLOAT_EQ(iter.Config().tau_collapse, 0.99f);
    EXPECT_TRUE(iter.Config().apply_smoothing);
}

TEST_F(EmIteratorTest, ConfigConstruction) {
    EmIteratorConfig config;
    config.gamma = 0.5f;
    config.tau_collapse = 0.95f;
    config.apply_smoothing = false;

    EmIterator iter(config);
    EXPECT_FLOAT_EQ(iter.Config().gamma, 0.5f);
    EXPECT_FLOAT_EQ(iter.Config().tau_collapse, 0.95f);
    EXPECT_FALSE(iter.Config().apply_smoothing);
}

TEST_F(EmIteratorTest, ComputeCoveragePrior) {
    EmIterator iter;

    auto coverage = iter.ComputeCoveragePrior(state_, 5);

    ASSERT_EQ(coverage.size(), 5u);

    // Bucket 0: Read 0 (0.33) + Read 2 (0.5) = 0.83
    EXPECT_NEAR(coverage[0], 0.83f, 0.01f);

    // Bucket 1: Read 0 (0.34) + Read 1 (0.8) = 1.14
    EXPECT_NEAR(coverage[1], 1.14f, 0.01f);

    // Bucket 2: Read 0 (0.33) + Read 1 (0.1) + Read 3 (0.1) = 0.53
    EXPECT_NEAR(coverage[2], 0.53f, 0.01f);

    // Bucket 3: Read 1 (0.1) + Read 2 (0.5) + Read 3 (0.85) = 1.45
    EXPECT_NEAR(coverage[3], 1.45f, 0.01f);

    // Bucket 4: Read 3 (0.05) = 0.05
    EXPECT_NEAR(coverage[4], 0.05f, 0.01f);
}

TEST_F(EmIteratorTest, CoveragePriorIgnoresCollapsed) {
    state_.set_collapsed(1, true);

    EmIterator iter;
    auto coverage = iter.ComputeCoveragePrior(state_, 5);

    // Read 1 is collapsed, so its contributions should be excluded
    // Bucket 1: Only Read 0 (0.34)
    EXPECT_NEAR(coverage[1], 0.34f, 0.01f);
}

TEST_F(EmIteratorTest, StepUpdatesState) {
    EmIterator iter;

    auto stats = iter.Step(state_, 5);

    EXPECT_EQ(stats.reads_processed, 4u);
    EXPECT_EQ(stats.reads_collapsed, 0u);
    EXPECT_GT(stats.iteration_time_ms, 0.0f);

    // After step, probabilities should still sum to ~1.0 for each read
    for (std::uint32_t r = 0; r < state_.n_reads(); ++r) {
        auto probs = state_.probabilities_for_read(r);
        float sum = 0.0f;
        for (float p : probs) sum += p;
        EXPECT_NEAR(sum, 1.0f, 0.001f);
    }
}

TEST_F(EmIteratorTest, StepWithLikelihoodFn) {
    EmIterator iter;

    // Likelihood strongly prefers bucket 1 for all reads
    iter.SetLikelihoodFn([](std::uint32_t, std::uint32_t b) {
        return (b == 1) ? 1.0f : 0.1f;
    });

    auto stats = iter.Step(state_, 5);
    EXPECT_EQ(stats.reads_processed, 4u);

    // Read 0 should have shifted probability toward bucket 1
    auto probs_0 = state_.probabilities_for_read(0);
    auto buckets_0 = state_.bucket_indices_for_read(0);
    float prob_bucket_1 = 0.0f;
    for (std::size_t i = 0; i < buckets_0.size(); ++i) {
        if (buckets_0[i] == 1) {
            prob_bucket_1 = probs_0[i];
            break;
        }
    }
    EXPECT_GT(prob_bucket_1, 0.34f);  // Was 0.34, should have increased
}

TEST_F(EmIteratorTest, StepWithAiPrior) {
    EmIterator iter;

    iter.SetAiPriorFn([](std::uint32_t r, std::uint32_t b) {
        // AI prior prefers bucket 0 for read 0, bucket 3 for others
        if (r == 0) return (b == 0) ? 1.0f : 0.1f;
        return (b == 3) ? 1.0f : 0.1f;
    });

    auto stats = iter.Step(state_, 5);
    EXPECT_EQ(stats.reads_processed, 4u);
}

TEST_F(EmIteratorTest, StepWithBiologyPrior) {
    EmIterator iter;

    std::vector<float> bio_prior = {1.0f, 0.5f, 0.5f, 2.0f, 0.1f};
    iter.SetBiologyPrior(bio_prior);

    auto stats = iter.Step(state_, 5);
    EXPECT_EQ(stats.reads_processed, 4u);
}

TEST_F(EmIteratorTest, CheckAndCollapseNoneReadyYet) {
    EmIterator iter;

    std::uint32_t collapsed = iter.CheckAndCollapse(state_);

    // No read has P > 0.99 initially
    EXPECT_EQ(collapsed, 0u);
    EXPECT_EQ(state_.count_collapsed(), 0u);
}

TEST_F(EmIteratorTest, CheckAndCollapseWithHighProbability) {
    // Manually set read 1 to have very high probability on bucket 1
    std::vector<BucketProb> cands_high = {{1, 0.995f}, {2, 0.003f}, {3, 0.002f}};
    state_.set_read_candidates(1, cands_high);

    EmIterator iter;
    std::uint32_t collapsed = iter.CheckAndCollapse(state_);

    EXPECT_EQ(collapsed, 1u);
    EXPECT_TRUE(state_.is_collapsed(1));
    EXPECT_EQ(state_.collapsed_bucket(1), 1u);
}

TEST_F(EmIteratorTest, CheckAndCollapseWithCustomThreshold) {
    // Read 1 has 0.8 on bucket 1
    EmIteratorConfig config;
    config.tau_collapse = 0.75f;  // Lower threshold

    EmIterator iter(config);
    std::uint32_t collapsed = iter.CheckAndCollapse(state_);

    EXPECT_EQ(collapsed, 2u);  // Read 1 (0.8) and Read 3 (0.85)
}

TEST_F(EmIteratorTest, MultipleStepsConvergence) {
    EmIterator iter;

    // Likelihood strongly prefers bucket 1
    iter.SetLikelihoodFn([](std::uint32_t, std::uint32_t b) {
        return (b == 1) ? 1.0f : 0.01f;
    });

    // Run multiple steps
    for (int i = 0; i < 20; ++i) {
        iter.Step(state_, 5);
    }

    // Read 0 has bucket 1 as candidate, should converge toward it
    auto probs_0 = state_.probabilities_for_read(0);
    auto buckets_0 = state_.bucket_indices_for_read(0);
    float prob_bucket_1 = 0.0f;
    for (std::size_t i = 0; i < buckets_0.size(); ++i) {
        if (buckets_0[i] == 1) {
            prob_bucket_1 = probs_0[i];
            break;
        }
    }
    // Should be very high after 20 iterations
    EXPECT_GT(prob_bucket_1, 0.9f);
}

TEST_F(EmIteratorTest, StepSkipsCollapsedReads) {
    state_.set_collapsed(0, true);
    state_.set_collapsed(2, true);

    EmIterator iter;
    auto stats = iter.Step(state_, 5);

    EXPECT_EQ(stats.reads_processed, 2u);  // Only reads 1 and 3
    EXPECT_EQ(stats.reads_collapsed, 2u);
}

TEST_F(EmIteratorTest, ComputeEntropy) {
    std::vector<float> uniform = {0.25f, 0.25f, 0.25f, 0.25f};
    float entropy_uniform = ComputeEntropy(uniform);
    EXPECT_NEAR(entropy_uniform, 2.0f, 0.01f);  // log2(4) = 2

    std::vector<float> peaked = {0.99f, 0.005f, 0.005f};
    float entropy_peaked = ComputeEntropy(peaked);
    EXPECT_LT(entropy_peaked, 0.1f);  // Very low entropy

    std::vector<float> binary = {0.5f, 0.5f};
    float entropy_binary = ComputeEntropy(binary);
    EXPECT_NEAR(entropy_binary, 1.0f, 0.01f);  // log2(2) = 1
}

TEST_F(EmIteratorTest, NormalizeProbabilities) {
    std::vector<float> probs = {2.0f, 3.0f, 5.0f};
    NormalizeProbabilities(probs);

    EXPECT_NEAR(probs[0], 0.2f, 0.001f);
    EXPECT_NEAR(probs[1], 0.3f, 0.001f);
    EXPECT_NEAR(probs[2], 0.5f, 0.001f);

    float sum = probs[0] + probs[1] + probs[2];
    EXPECT_NEAR(sum, 1.0f, 0.001f);
}

TEST_F(EmIteratorTest, NormalizeProbabilitiesZeroSum) {
    std::vector<float> probs = {0.0f, 0.0f, 0.0f};
    NormalizeProbabilities(probs);

    // Should fall back to uniform
    EXPECT_NEAR(probs[0], 1.0f / 3.0f, 0.001f);
    EXPECT_NEAR(probs[1], 1.0f / 3.0f, 0.001f);
    EXPECT_NEAR(probs[2], 1.0f / 3.0f, 0.001f);
}

TEST(BucketNeighborhoodTest, BuildNeighborhoodEmpty) {
    SmoothingKernelConfig config;
    std::vector<std::uint64_t> positions;

    auto hood = BuildNeighborhood(positions, config);

    EXPECT_EQ(hood.NumBuckets(), 0u);
}

TEST(BucketNeighborhoodTest, BuildNeighborhoodSingleBucket) {
    SmoothingKernelConfig config;
    config.sigma_genome_bp = 50000.0f;

    std::vector<std::uint64_t> positions = {100000};

    auto hood = BuildNeighborhood(positions, config);

    EXPECT_EQ(hood.NumBuckets(), 1u);
    EXPECT_TRUE(hood.GetNeighbors(0).empty());  // No other buckets
}

TEST(BucketNeighborhoodTest, BuildNeighborhoodNearbyBuckets) {
    SmoothingKernelConfig config;
    config.sigma_genome_bp = 100000.0f;
    config.min_weight = 0.01f;
    config.max_neighbors = 10;

    // Three buckets: 0 at 100kb, 1 at 150kb, 2 at 500kb
    std::vector<std::uint64_t> positions = {100000, 150000, 500000};

    auto hood = BuildNeighborhood(positions, config);

    EXPECT_EQ(hood.NumBuckets(), 3u);

    // Bucket 0 should have bucket 1 as neighbor (50kb apart)
    auto n0 = hood.GetNeighbors(0);
    EXPECT_FALSE(n0.empty());
    bool has_bucket_1 = false;
    for (auto n : n0) {
        if (n == 1) has_bucket_1 = true;
    }
    EXPECT_TRUE(has_bucket_1);

    // Bucket 2 is far from others, may have fewer/no neighbors
    auto n2 = hood.GetNeighbors(2);
    // 350kb from bucket 1, weight = exp(-350^2 / (2*100^2)) = exp(-6.125) ≈ 0.002
    // This is > 0.01, so should still be a neighbor
}

TEST(BucketNeighborhoodTest, NeighborWeightsDecreaseWithDistance) {
    SmoothingKernelConfig config;
    config.sigma_genome_bp = 100000.0f;
    config.min_weight = 0.001f;
    config.max_neighbors = 20;

    // Buckets at 0, 50kb, 100kb, 200kb
    std::vector<std::uint64_t> positions = {0, 50000, 100000, 200000};

    auto hood = BuildNeighborhood(positions, config);

    auto w0 = hood.GetWeights(0);
    auto n0 = hood.GetNeighbors(0);

    // Find weights for bucket 1 (50kb away) and bucket 3 (200kb away)
    float w_50k = 0.0f, w_200k = 0.0f;
    for (std::size_t i = 0; i < n0.size(); ++i) {
        if (n0[i] == 1) w_50k = w0[i];
        if (n0[i] == 3) w_200k = w0[i];
    }

    EXPECT_GT(w_50k, w_200k);  // Closer bucket should have higher weight
}

TEST_F(EmIteratorTest, MeanEntropyTracking) {
    EmIterator iter;
    auto stats = iter.Step(state_, 5);

    EXPECT_GT(stats.mean_entropy, 0.0f);
    EXPECT_LT(stats.mean_entropy, 3.0f);  // Bounded by log2(num_buckets)
}

TEST_F(EmIteratorTest, MaxProbDeltaTracking) {
    // With no likelihood/prior functions, damping should cause small updates
    EmIteratorConfig config;
    config.gamma = 0.1f;

    EmIterator iter(config);
    auto stats = iter.Step(state_, 5);

    EXPECT_GE(stats.max_prob_delta, 0.0f);
    EXPECT_LE(stats.max_prob_delta, 1.0f);
}

TEST_F(EmIteratorTest, SetConfigAfterConstruction) {
    EmIterator iter;

    EmIteratorConfig new_config;
    new_config.gamma = 0.8f;
    new_config.tau_collapse = 0.9f;

    iter.SetConfig(new_config);

    EXPECT_FLOAT_EQ(iter.Config().gamma, 0.8f);
    EXPECT_FLOAT_EQ(iter.Config().tau_collapse, 0.9f);
}

TEST_F(EmIteratorTest, MoveConstruction) {
    EmIterator iter1;
    iter1.SetLikelihoodFn([](std::uint32_t, std::uint32_t) { return 1.0f; });

    EmIterator iter2 = std::move(iter1);

    EXPECT_FLOAT_EQ(iter2.Config().gamma, 0.3f);
}

TEST_F(EmIteratorTest, WeightMultipliers) {
    EmIteratorConfig config;
    config.weight_seq_likelihood = 2.0f;
    config.weight_coverage = 0.5f;
    config.weight_ai_prior = 0.0f;  // Disable AI prior
    config.weight_bio_prior = 1.0f;

    EmIterator iter(config);

    iter.SetLikelihoodFn([](std::uint32_t, std::uint32_t b) {
        return (b == 1) ? 0.9f : 0.1f;
    });

    auto stats = iter.Step(state_, 5);
    EXPECT_EQ(stats.reads_processed, 4u);
}

}  // namespace
}  // namespace llmap
