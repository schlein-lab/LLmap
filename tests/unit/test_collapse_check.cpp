// LLmap — Unit tests for collapse_check module.

#include "reference_collapse/collapse_check.h"
#include "core/wave_state.h"

#include <gtest/gtest.h>
#include <cmath>

namespace llmap {
namespace {

// Test fixture with helper methods
class CollapseCheckTest : public ::testing::Test {
protected:
    // Create a WaveState with specified probability distributions
    WaveState CreateStateWithProbs(
        const std::vector<std::vector<std::pair<uint32_t, float>>>& read_probs)
    {
        WaveState state(static_cast<uint32_t>(read_probs.size()));

        for (uint32_t r = 0; r < read_probs.size(); ++r) {
            std::vector<BucketProb> candidates;
            for (const auto& [bucket, prob] : read_probs[r]) {
                candidates.push_back({bucket, prob});
            }
            // Sort by bucket ID as required
            std::sort(candidates.begin(), candidates.end(),
                [](const auto& a, const auto& b) {
                    return a.bucket_id < b.bucket_id;
                });
            state.set_read_candidates(r, candidates);
        }

        return state;
    }
};

// --- Basic functionality tests ---

TEST_F(CollapseCheckTest, EmptyStateReturnsZeroStats) {
    WaveState state(0);
    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 0u);
    EXPECT_EQ(stats.newly_collapsed, 0u);
    EXPECT_EQ(stats.already_collapsed, 0u);
}

TEST_F(CollapseCheckTest, SingleReadHighProbCollapses) {
    // Single read with max probability 0.99
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.005f}, {2, 0.005f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 1u);
    EXPECT_EQ(stats.newly_collapsed, 1u);
    EXPECT_EQ(stats.active_after, 0u);
    EXPECT_TRUE(state.is_collapsed(0));
    EXPECT_EQ(state.collapsed_bucket(0), 0u);
}

TEST_F(CollapseCheckTest, SingleReadBelowThresholdRemainActive) {
    // Single read with max probability 0.90 (below default 0.99)
    auto state = CreateStateWithProbs({
        {{0, 0.90f}, {1, 0.05f}, {2, 0.05f}}
    });

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 1u);
    EXPECT_EQ(stats.newly_collapsed, 0u);
    EXPECT_EQ(stats.active_after, 1u);
    EXPECT_FALSE(state.is_collapsed(0));
}

TEST_F(CollapseCheckTest, CustomThresholdWorks) {
    auto state = CreateStateWithProbs({
        {{0, 0.80f}, {1, 0.15f}, {2, 0.05f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.80f;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.newly_collapsed, 1u);
    EXPECT_TRUE(state.is_collapsed(0));
}

TEST_F(CollapseCheckTest, MultipleReadsWithMixedConvergence) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.005f}, {2, 0.005f}},  // Should collapse
        {{0, 0.50f}, {1, 0.30f}, {2, 0.20f}},    // Should remain active
        {{0, 0.995f}, {1, 0.005f}},               // Should collapse
        {{0, 0.70f}, {1, 0.30f}},                 // Should remain active
    });

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 4u);
    EXPECT_EQ(stats.newly_collapsed, 2u);
    EXPECT_EQ(stats.active_after, 2u);
    EXPECT_TRUE(state.is_collapsed(0));
    EXPECT_FALSE(state.is_collapsed(1));
    EXPECT_TRUE(state.is_collapsed(2));
    EXPECT_FALSE(state.is_collapsed(3));
}

// --- Dropout classification tests ---

TEST_F(CollapseCheckTest, AmbiguousReadDetected) {
    // Two candidates with similar probability
    auto state = CreateStateWithProbs({
        {{0, 0.51f}, {1, 0.49f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    config.tau_ambiguous = 0.05f;
    CollapseChecker checker(config);

    auto info = checker.AnalyzeRead(state, 0);

    EXPECT_EQ(info.dropout_type, DropoutType::Ambiguous);
    EXPECT_FLOAT_EQ(info.max_probability, 0.51f);
    EXPECT_FLOAT_EQ(info.second_probability, 0.49f);
}

TEST_F(CollapseCheckTest, LowConfidenceReadDetected) {
    // Max probability below low-confidence threshold
    auto state = CreateStateWithProbs({
        {{0, 0.40f}, {1, 0.30f}, {2, 0.30f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    config.tau_low_confidence = 0.50f;
    CollapseChecker checker(config);

    auto info = checker.AnalyzeRead(state, 0);

    EXPECT_EQ(info.dropout_type, DropoutType::LowConfidence);
}

TEST_F(CollapseCheckTest, OrphanedReadDetected) {
    // Empty probability distribution
    WaveState state(1);  // Read with no candidates

    CollapseChecker checker;
    auto info = checker.AnalyzeRead(state, 0);

    EXPECT_EQ(info.dropout_type, DropoutType::Orphaned);
}

TEST_F(CollapseCheckTest, DropoutStatsAccumulate) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},                 // Converged
        {{0, 0.51f}, {1, 0.49f}},                 // Ambiguous
        {{0, 0.40f}, {1, 0.30f}, {2, 0.30f}},    // Low confidence
        {{0, 0.80f}, {1, 0.20f}},                 // Active
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    config.tau_ambiguous = 0.05f;
    config.tau_low_confidence = 0.50f;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 4u);
    EXPECT_EQ(stats.newly_collapsed, 1u);
    EXPECT_EQ(stats.ambiguous, 1u);
    EXPECT_EQ(stats.low_confidence, 1u);
    EXPECT_EQ(stats.active_after, 1u);
}

// --- Already-collapsed reads ---

TEST_F(CollapseCheckTest, AlreadyCollapsedReadsAreSkipped) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},
        {{0, 0.99f}, {1, 0.01f}},
    });

    // Pre-collapse first read
    state.collapse_read(0);

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_EQ(stats.total_reads, 2u);
    EXPECT_EQ(stats.already_collapsed, 1u);
    EXPECT_EQ(stats.newly_collapsed, 1u);
}

// --- Auto-collapse control ---

TEST_F(CollapseCheckTest, AutoCollapseDisabledDoesNotModifyState) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}}
    });

    CollapseCheckConfig config;
    config.auto_collapse = false;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.newly_collapsed, 1u);  // Detected as converged
    EXPECT_FALSE(state.is_collapsed(0));   // But not actually collapsed
}

// --- Detailed check tests ---

TEST_F(CollapseCheckTest, DetailedCheckReturnsPerReadInfo) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},
        {{0, 0.60f}, {1, 0.40f}},
    });

    CollapseChecker checker;
    std::vector<ReadCollapseInfo> read_info;

    auto stats = checker.CheckDetailed(state, read_info);

    ASSERT_EQ(read_info.size(), 2u);

    EXPECT_EQ(read_info[0].read_idx, 0u);
    EXPECT_EQ(read_info[0].dropout_type, DropoutType::Converged);
    EXPECT_FLOAT_EQ(read_info[0].max_probability, 0.99f);
    EXPECT_EQ(read_info[0].top_bucket, 0u);

    EXPECT_EQ(read_info[1].read_idx, 1u);
    EXPECT_EQ(read_info[1].dropout_type, DropoutType::None);
}

// --- Entropy calculation tests ---

TEST_F(CollapseCheckTest, EntropyIsZeroForDeterministic) {
    auto state = CreateStateWithProbs({
        {{0, 1.0f}}
    });

    CollapseChecker checker;
    auto info = checker.AnalyzeRead(state, 0);

    EXPECT_NEAR(info.entropy, 0.0f, 1e-6f);
}

TEST_F(CollapseCheckTest, EntropyIsMaxForUniform) {
    auto state = CreateStateWithProbs({
        {{0, 0.25f}, {1, 0.25f}, {2, 0.25f}, {3, 0.25f}}
    });

    CollapseChecker checker;
    auto info = checker.AnalyzeRead(state, 0);

    // Max entropy for 4 outcomes is log2(4) = 2.0
    EXPECT_NEAR(info.entropy, 2.0f, 1e-5f);
}

// --- Statistics tests ---

TEST_F(CollapseCheckTest, ProbabilityStatsAreCorrect) {
    auto state = CreateStateWithProbs({
        {{0, 0.80f}, {1, 0.20f}},
        {{0, 0.90f}, {1, 0.10f}},
        {{0, 0.70f}, {1, 0.30f}},
    });

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_NEAR(stats.mean_max_prob, 0.80f, 1e-5f);  // (0.8 + 0.9 + 0.7) / 3
    EXPECT_FLOAT_EQ(stats.min_max_prob, 0.70f);
    EXPECT_FLOAT_EQ(stats.max_max_prob, 0.90f);
}

TEST_F(CollapseCheckTest, CollapseRateCalculation) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},  // Collapses
        {{0, 0.99f}, {1, 0.01f}},  // Collapses
        {{0, 0.50f}, {1, 0.50f}},  // Active
        {{0, 0.50f}, {1, 0.50f}},  // Active
    });

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_FLOAT_EQ(stats.collapse_rate(), 0.5f);  // 2/4
}

// --- Convenience function tests ---

TEST_F(CollapseCheckTest, CheckAndCollapseConvenienceWorks) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},
        {{0, 0.50f}, {1, 0.50f}},
    });

    auto stats = CheckAndCollapse(state, 0.99f);

    EXPECT_EQ(stats.newly_collapsed, 1u);
    EXPECT_TRUE(state.is_collapsed(0));
}

TEST_F(CollapseCheckTest, GetDropoutTypesConvenienceWorks) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},  // Converged
        {{0, 0.51f}, {1, 0.49f}},  // Ambiguous
        {{0, 0.80f}, {1, 0.20f}},  // None (active)
    });

    auto types = GetDropoutTypes(state, 0.99f, 0.05f);

    ASSERT_EQ(types.size(), 3u);
    EXPECT_EQ(types[0], DropoutType::Converged);
    EXPECT_EQ(types[1], DropoutType::Ambiguous);
    EXPECT_EQ(types[2], DropoutType::None);
}

TEST_F(CollapseCheckTest, GetReadsByDropoutTypeWorks) {
    std::vector<DropoutType> types = {
        DropoutType::Converged,
        DropoutType::Ambiguous,
        DropoutType::Converged,
        DropoutType::None,
        DropoutType::Ambiguous,
    };

    auto converged = GetReadsByDropoutType(types, DropoutType::Converged);
    auto ambiguous = GetReadsByDropoutType(types, DropoutType::Ambiguous);

    EXPECT_EQ(converged.size(), 2u);
    EXPECT_EQ(converged[0], 0u);
    EXPECT_EQ(converged[1], 2u);

    EXPECT_EQ(ambiguous.size(), 2u);
    EXPECT_EQ(ambiguous[0], 1u);
    EXPECT_EQ(ambiguous[1], 4u);
}

// --- Edge case tests ---

TEST_F(CollapseCheckTest, SingleCandidateAlwaysCollapses) {
    auto state = CreateStateWithProbs({
        {{5, 1.0f}}  // Single candidate
    });

    CollapseChecker checker;
    auto stats = checker.Check(state);

    EXPECT_EQ(stats.newly_collapsed, 1u);
    EXPECT_TRUE(state.is_collapsed(0));
    EXPECT_EQ(state.collapsed_bucket(0), 5u);
}

TEST_F(CollapseCheckTest, VerySmallProbabilitiesHandled) {
    auto state = CreateStateWithProbs({
        {{0, 1e-10f}, {1, 1e-11f}}
    });

    CollapseChecker checker;
    auto info = checker.AnalyzeRead(state, 0);

    // Should be treated as orphaned (near-zero max prob)
    EXPECT_EQ(info.dropout_type, DropoutType::Orphaned);
}

TEST_F(CollapseCheckTest, ExactThresholdIsInclusive) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;  // Exact match
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.newly_collapsed, 1u);
}

TEST_F(CollapseCheckTest, JustBelowThresholdDoesNotCollapse) {
    auto state = CreateStateWithProbs({
        {{0, 0.9899f}, {1, 0.0101f}}
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    EXPECT_EQ(stats.newly_collapsed, 0u);
    EXPECT_EQ(stats.active_after, 1u);
}

TEST_F(CollapseCheckTest, DropoutRateCalculation) {
    auto state = CreateStateWithProbs({
        {{0, 0.99f}, {1, 0.01f}},                 // Converged
        {{0, 0.51f}, {1, 0.49f}},                 // Ambiguous
        {{0, 0.40f}, {1, 0.30f}, {2, 0.30f}},    // Low confidence
        {{0, 0.80f}, {1, 0.20f}},                 // Active
    });

    CollapseCheckConfig config;
    config.tau_collapse = 0.99f;
    config.tau_ambiguous = 0.05f;
    config.tau_low_confidence = 0.50f;
    CollapseChecker checker(config);

    auto stats = checker.Check(state);

    // Dropouts = ambiguous + low_confidence + orphaned = 1 + 1 + 0 = 2
    // Dropout rate = 2/4 = 0.5
    EXPECT_FLOAT_EQ(stats.dropout_rate(), 0.5f);
}

// --- ComputeReadEntropy standalone tests ---

TEST_F(CollapseCheckTest, ComputeReadEntropyEmptyInput) {
    std::vector<float> empty;
    float entropy = ComputeReadEntropy(empty);
    EXPECT_FLOAT_EQ(entropy, 0.0f);
}

TEST_F(CollapseCheckTest, ComputeReadEntropyBinaryEven) {
    std::vector<float> probs = {0.5f, 0.5f};
    float entropy = ComputeReadEntropy(probs);
    EXPECT_NEAR(entropy, 1.0f, 1e-5f);  // log2(2) = 1
}

}  // namespace
}  // namespace llmap
