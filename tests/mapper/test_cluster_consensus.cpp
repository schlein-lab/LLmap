// LLmap — cluster_consensus (assemble-then-map V1) tests.

#include "mapper/cluster_consensus.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace {

using llmap::mapper::AssembleConsensus;
using llmap::mapper::ConsensusConfig;

// A non-repetitive 50 bp backbone (so 15-mers anchor uniquely).
const std::string kS = "ACGTAGCTAGGCATTCGATCCGTAGCATGCATTAGCCGATTACGGATCGA";

ConsensusConfig Cfg() { return ConsensusConfig{/*k=*/15, /*min_shared=*/3}; }

TEST(ClusterConsensus, SingleReadIsItself) {
    std::vector<std::string> reads = {kS};
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.sequence, kS);
    ASSERT_EQ(c.members.size(), 1u);
    EXPECT_TRUE(c.members[0].anchored);
    EXPECT_EQ(c.n_anchored, 1u);
}

TEST(ClusterConsensus, MajorityCorrectsAnError) {
    std::string mut = kS;
    mut[25] = (mut[25] == 'A') ? 'T' : 'A';  // single minority error at pos 25
    std::vector<std::string> reads = {kS, kS, mut};  // 2 correct, 1 mutated
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.sequence, kS);          // majority restores the true base
    EXPECT_EQ(c.n_anchored, 3u);
}

TEST(ClusterConsensus, SubstringAnchoredAtOffset) {
    const std::string sub = kS.substr(10, 35);  // offset 10
    std::vector<std::string> reads = {kS, sub};
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.sequence.size(), kS.size());  // backbone spans the contig
    ASSERT_EQ(c.members.size(), 2u);
    // The substring placed at offset 10 (backbone is at 0).
    EXPECT_TRUE(c.members[1].anchored);
    EXPECT_EQ(c.members[1].offset, 10);
}

TEST(ClusterConsensus, UnrelatedReadNotAnchored) {
    const std::string other = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGG";  // no shared 15-mer
    std::vector<std::string> reads = {kS, other};
    const auto c = AssembleConsensus(reads, Cfg());
    ASSERT_EQ(c.members.size(), 2u);
    EXPECT_TRUE(c.members[0].anchored);
    EXPECT_FALSE(c.members[1].anchored);  // kept as singleton → caller maps it alone
}

TEST(ClusterConsensus, EmptyCluster) {
    const auto c = AssembleConsensus({}, Cfg());
    EXPECT_TRUE(c.sequence.empty());
    EXPECT_TRUE(c.members.empty());
}

}  // namespace
