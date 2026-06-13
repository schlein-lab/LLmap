// LLmap — merge_ambiguity (WaveCollapse read-vs-read merge guard) tests.

#include "mapper/merge_ambiguity.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

namespace {

using llmap::mapper::AmbiguityConfig;
using llmap::mapper::AssessDiagonal;

AmbiguityConfig Cfg() { return AmbiguityConfig{}; }  // tol 20, frac 0.60, min 3

TEST(MergeAmbiguity, Empty) {
    const auto d = AssessDiagonal({}, Cfg());
    EXPECT_EQ(d.shared, 0u);
    EXPECT_FALSE(d.clear);
}

TEST(MergeAmbiguity, ClearOverlapCollapsesToOneDiagonal) {
    // A real overlap: every anchor agrees on offset ~100.
    std::vector<std::int64_t> off = {98, 100, 101, 99, 102, 100};
    const auto d = AssessDiagonal(off, Cfg());
    EXPECT_EQ(d.shared, 6u);
    EXPECT_DOUBLE_EQ(d.dominant_fraction, 1.0);
    EXPECT_LT(d.offset_entropy, 0.25);  // mass collapsed → low entropy
    EXPECT_TRUE(d.clear);
}

TEST(MergeAmbiguity, BlurryMatchAcrossManyDiagonalsRejected) {
    // Spurious repeat/paralog match: anchors scattered over many offsets.
    std::vector<std::int64_t> off = {10, 500, 1000, 50, 800, 300, 1500};
    const auto d = AssessDiagonal(off, Cfg());
    EXPECT_LT(d.dominant_fraction, 0.60);
    EXPECT_GT(d.offset_entropy, 0.75);  // mass spread → high entropy
    EXPECT_FALSE(d.clear);              // do NOT merge
}

TEST(MergeAmbiguity, MostlyOneDiagonalWithFewOutliersStillClear) {
    // Real overlap with a couple of spurious anchors — still collapses.
    std::vector<std::int64_t> off = {100, 101, 99, 100, 102, 700};
    const auto d = AssessDiagonal(off, Cfg());
    EXPECT_GT(d.dominant_fraction, 0.60);  // 5/6 on the diagonal
    EXPECT_TRUE(d.clear);
}

TEST(MergeAmbiguity, TooFewAnchorsCannotJudge) {
    // Concentrated but below min_anchors → not declared clear (conservative).
    std::vector<std::int64_t> off = {100, 101};
    const auto d = AssessDiagonal(off, Cfg());
    EXPECT_DOUBLE_EQ(d.dominant_fraction, 1.0);
    EXPECT_FALSE(d.clear);  // gated by min_anchors=3
}

}  // namespace
