#include <gtest/gtest.h>

#include "classical/chain.h"
#include "classical/minimizer_index.h"

using namespace llmap::classical;

class ChainTest : public ::testing::Test {
protected:
    ChainConfig config;

    void SetUp() override {
        config = ChainConfig{};
    }
};

TEST_F(ChainTest, EmptyAnchorsProducesEmptyResult) {
    std::vector<Anchor> anchors;
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    EXPECT_TRUE(result.chains.empty());
    EXPECT_EQ(result.best_score, 0);
    EXPECT_EQ(result.total_anchors, 0);
}

TEST_F(ChainTest, SingleAnchorTooFewForChain) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true}
    };

    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    EXPECT_TRUE(result.chains.empty());
    EXPECT_EQ(result.total_anchors, 1);
}

TEST_F(ChainTest, TwoAnchorsStillTooFew) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true}
    };

    config.min_chain_anchors = 3;  // Need at least 3
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    EXPECT_TRUE(result.chains.empty());
}

TEST_F(ChainTest, ThreeColinearAnchorsFormChain) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    ASSERT_EQ(result.chains.size(), 1);
    EXPECT_EQ(result.chains[0].NumAnchors(), 3);
    EXPECT_EQ(result.chains[0].ref_id, 0);
    EXPECT_TRUE(result.chains[0].is_forward);
    EXPECT_GT(result.chains[0].score, 0);
}

TEST_F(ChainTest, ForwardStrandChainCorrectOrder) {
    std::vector<Anchor> anchors = {
        {0, 100, 50, true},
        {0, 200, 150, true},
        {0, 300, 250, true},
        {0, 400, 350, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    ASSERT_EQ(result.chains.size(), 1);
    const auto& chain = result.chains[0];

    EXPECT_EQ(chain.ref_start, 100);
    EXPECT_EQ(chain.ref_end, 400);
    EXPECT_EQ(chain.query_start, 50);
    EXPECT_EQ(chain.query_end, 350);
    EXPECT_TRUE(chain.is_forward);
}

TEST_F(ChainTest, ReverseStrandChain) {
    // Reverse strand: ref increases, query decreases
    std::vector<Anchor> anchors = {
        {0, 100, 350, false},
        {0, 200, 250, false},
        {0, 300, 150, false},
        {0, 400, 50, false}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    ASSERT_EQ(result.chains.size(), 1);
    EXPECT_FALSE(result.chains[0].is_forward);
    EXPECT_EQ(result.chains[0].NumAnchors(), 4);
}

TEST_F(ChainTest, DifferentRefsProduceSeparateChains) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true},
        {1, 100, 100, true},
        {1, 200, 200, true},
        {1, 300, 300, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    ASSERT_EQ(result.chains.size(), 2);
    // Check both refs are represented
    std::set<uint32_t> refs;
    for (const auto& chain : result.chains) {
        refs.insert(chain.ref_id);
    }
    EXPECT_EQ(refs.size(), 2);
}

TEST_F(ChainTest, DifferentStrandsProduceSeparateChains) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true},
        {0, 400, 350, false},
        {0, 500, 250, false},
        {0, 600, 150, false}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    ASSERT_EQ(result.chains.size(), 2);
    bool has_forward = false, has_reverse = false;
    for (const auto& chain : result.chains) {
        if (chain.is_forward) has_forward = true;
        else has_reverse = true;
    }
    EXPECT_TRUE(has_forward);
    EXPECT_TRUE(has_reverse);
}

TEST_F(ChainTest, NonColinearAnchorsNotChained) {
    // Query positions not increasing with ref positions
    std::vector<Anchor> anchors = {
        {0, 100, 300, true},  // Out of order in query
        {0, 200, 200, true},
        {0, 300, 100, true}
    };

    config.min_chain_anchors = 2;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    // Should not form a 3-anchor chain because they're not colinear
    for (const auto& chain : result.chains) {
        EXPECT_LT(chain.NumAnchors(), 3);
    }
}

TEST_F(ChainTest, GapTooLargeBreaksChain) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true},
        {0, 10000, 10000, true},  // Gap > max_gap_ref
        {0, 10100, 10100, true},
        {0, 10200, 10200, true}
    };

    config.max_gap_ref = 5000;
    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 20000, config);

    // Should produce two separate chains
    EXPECT_EQ(result.chains.size(), 2);
}

TEST_F(ChainTest, ScoreFilteringKeepsTopChains) {
    // Create two chains with different scores
    std::vector<Anchor> anchors = {
        // Strong chain with many anchors
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true},
        {0, 400, 400, true},
        {0, 500, 500, true},
        // Weaker chain with fewer anchors
        {1, 100, 100, true},
        {1, 200, 200, true},
        {1, 300, 300, true}
    };

    config.min_chain_anchors = 3;
    config.min_score_fraction = 0.9f;
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    // First chain should have higher score
    ASSERT_GE(result.chains.size(), 1);
    EXPECT_EQ(result.best_score, result.chains[0].score);
}

TEST_F(ChainTest, AnchorPairScoreBasicForward) {
    Anchor a{0, 100, 100, true};
    Anchor b{0, 200, 200, true};

    int32_t score = AnchorPairScore(a, b, config);
    EXPECT_GT(score, 0);
}

TEST_F(ChainTest, AnchorPairScoreNegativeForDifferentRefs) {
    Anchor a{0, 100, 100, true};
    Anchor b{1, 200, 200, true};

    int32_t score = AnchorPairScore(a, b, config);
    EXPECT_EQ(score, std::numeric_limits<int32_t>::min());
}

TEST_F(ChainTest, AnchorPairScoreNegativeForDifferentStrands) {
    Anchor a{0, 100, 100, true};
    Anchor b{0, 200, 200, false};

    int32_t score = AnchorPairScore(a, b, config);
    EXPECT_EQ(score, std::numeric_limits<int32_t>::min());
}

TEST_F(ChainTest, AnchorPairScoreNegativeForBackwards) {
    Anchor a{0, 200, 200, true};
    Anchor b{0, 100, 100, true};  // Before a

    int32_t score = AnchorPairScore(a, b, config);
    EXPECT_EQ(score, std::numeric_limits<int32_t>::min());
}

TEST_F(ChainTest, IsColinearMatchesScore) {
    Anchor a{0, 100, 100, true};
    Anchor b{0, 200, 200, true};
    Anchor c{0, 200, 100, true};  // Not colinear (query doesn't increase)

    EXPECT_TRUE(IsColinear(a, b, config));
    EXPECT_FALSE(IsColinear(a, c, config));
}

TEST_F(ChainTest, ExtractChainsFromHits) {
    std::vector<MinimizerHit> hits = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true},
        {0, 400, 400, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChains(hits, 500, config);

    ASSERT_EQ(result.chains.size(), 1);
    EXPECT_EQ(result.chains[0].NumAnchors(), 4);
}

TEST_F(ChainTest, ChainResultHasCorrectMetadata) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 200, 200, true},
        {0, 300, 300, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    EXPECT_EQ(result.total_anchors, 3);
    EXPECT_GE(result.chain_time_ms, 0.0f);
}

TEST_F(ChainTest, MinChainScoreFilter) {
    std::vector<Anchor> anchors = {
        {0, 100, 100, true},
        {0, 110, 110, true},
        {0, 120, 120, true}
    };

    // With high min_chain_score, weak chains are filtered
    config.min_chain_anchors = 3;
    config.min_chain_score = 1000;  // Very high threshold
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    EXPECT_TRUE(result.chains.empty());
}

TEST_F(ChainTest, ChainSpanCalculation) {
    std::vector<Anchor> anchors = {
        {0, 100, 50, true},
        {0, 200, 150, true},
        {0, 400, 350, true}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 500, config);

    ASSERT_EQ(result.chains.size(), 1);
    const auto& chain = result.chains[0];

    EXPECT_EQ(chain.RefSpan(), 300);   // 400 - 100
    EXPECT_EQ(chain.QuerySpan(), 300); // 350 - 50
}

TEST_F(ChainTest, LargeGapDifferenceReducesScore) {
    // Two chains: one with matching gaps, one with large gap difference
    Anchor a1{0, 100, 100, true};
    Anchor b1{0, 200, 200, true};  // Gap: ref=100, query=100

    Anchor a2{0, 100, 100, true};
    Anchor b2{0, 200, 150, true};  // Gap: ref=100, query=50

    int32_t score1 = AnchorPairScore(a1, b1, config);
    int32_t score2 = AnchorPairScore(a2, b2, config);

    // Equal gaps should score better than unequal gaps
    EXPECT_GT(score1, score2);
}

TEST_F(ChainTest, MultipleChainsPreserveOrder) {
    std::vector<Anchor> anchors;
    // Create 3 chains with varying strengths
    for (int ref = 0; ref < 3; ++ref) {
        size_t num_anchors = 5 - ref;  // 5, 4, 3 anchors
        for (size_t i = 0; i < num_anchors; ++i) {
            anchors.push_back({
                static_cast<uint32_t>(ref),
                static_cast<uint32_t>(100 + i * 100),
                static_cast<uint32_t>(100 + i * 100),
                true
            });
        }
    }

    config.min_chain_anchors = 3;
    config.min_score_fraction = 0.5f;
    auto result = ExtractChainsFromAnchors(anchors, 1000, config);

    // All 3 chains should pass with 0.5 threshold
    EXPECT_EQ(result.chains.size(), 3);

    // Chains should be sorted by score descending
    for (size_t i = 1; i < result.chains.size(); ++i) {
        EXPECT_LE(result.chains[i].score, result.chains[i-1].score);
    }
}

TEST_F(ChainTest, ReverseStrandQueryCoordinates) {
    std::vector<Anchor> anchors = {
        {0, 100, 300, false},
        {0, 200, 200, false},
        {0, 300, 100, false}
    };

    config.min_chain_anchors = 3;
    auto result = ExtractChainsFromAnchors(anchors, 400, config);

    ASSERT_EQ(result.chains.size(), 1);
    const auto& chain = result.chains[0];

    // For reverse strand, query_start < query_end (sorted)
    EXPECT_LE(chain.query_start, chain.query_end);
}
