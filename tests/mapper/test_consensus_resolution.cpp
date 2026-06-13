// LLmap — consensus_resolution (multi-level agreement) tests.

#include "mapper/consensus_resolution.h"

#include <gtest/gtest.h>

#include <vector>

namespace {

using llmap::mapper::ConsensusConfig;
using llmap::mapper::IsConsensusResolved;
using llmap::mapper::LevelPlacement;

LevelPlacement P(std::string ref, std::uint64_t pos) {
    return LevelPlacement{true, std::move(ref), pos};
}
LevelPlacement Unmapped() { return LevelPlacement{false, "", 0}; }

TEST(ConsensusResolution, LastTwoAgreeResolves) {
    std::vector<LevelPlacement> levels = {
        P("chr1", 5000),     // coarse — different (still settling)
        P("chr1", 1000),     // fine
        P("chr1", 1003),     // finer — agrees with prev within tol
    };
    const auto r = IsConsensusResolved(levels, {/*agree_k=*/2, /*tol=*/10});
    EXPECT_TRUE(r.resolved);
    EXPECT_EQ(r.agreeing_levels, 2u);
}

TEST(ConsensusResolution, LastTwoDisagreeEscalates) {
    std::vector<LevelPlacement> levels = {
        P("chr1", 1000),
        P("chr1", 1002),
        P("chr2", 9000),  // jumps → no trailing agreement
    };
    EXPECT_FALSE(IsConsensusResolved(levels, {2, 10}).resolved);
}

TEST(ConsensusResolution, AllThreeAgree) {
    std::vector<LevelPlacement> levels = {P("chr1", 1000), P("chr1", 1001),
                                          P("chr1", 1000)};
    const auto r = IsConsensusResolved(levels, {3, 10});
    EXPECT_TRUE(r.resolved);
    EXPECT_EQ(r.agreeing_levels, 3u);
}

TEST(ConsensusResolution, NotEnoughLevels) {
    std::vector<LevelPlacement> levels = {P("chr1", 1000)};
    EXPECT_FALSE(IsConsensusResolved(levels, {2, 10}).resolved);
}

TEST(ConsensusResolution, UnmappedLastBlocks) {
    std::vector<LevelPlacement> levels = {P("chr1", 1000), Unmapped()};
    EXPECT_FALSE(IsConsensusResolved(levels, {2, 10}).resolved);
}

TEST(ConsensusResolution, OutsideToleranceDoesNotAgree) {
    std::vector<LevelPlacement> levels = {P("chr1", 1000), P("chr1", 1100)};
    EXPECT_FALSE(IsConsensusResolved(levels, {2, 10}).resolved);  // 100 > tol 10
}

}  // namespace
