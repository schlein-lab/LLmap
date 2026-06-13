// LLmap — cluster_consensus (assemble-then-map, Line A) tests.
//
// probe = transient search probe (smoothing allowed, region-finding only).
// members[]/offset/strand + depth[] = lossless layout (the real output).

#include "mapper/cluster_consensus.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>

namespace {

using llmap::mapper::AssembleConsensus;
using llmap::mapper::ConsensusConfig;

// A non-repetitive 50 bp backbone (so 15-mers anchor uniquely).
const std::string kS = "ACGTAGCTAGGCATTCGATCCGTAGCATGCATTAGCCGATTACGGATCGA";

ConsensusConfig Cfg() { return ConsensusConfig{/*k=*/15, /*min_shared=*/3}; }

std::string RevComp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
        }
    }
    return r;
}

TEST(ClusterConsensus, SingleReadIsItself) {
    std::vector<std::string> reads = {kS};
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.probe, kS);
    ASSERT_EQ(c.members.size(), 1u);
    EXPECT_TRUE(c.members[0].anchored);
    EXPECT_FALSE(c.members[0].reverse);
    EXPECT_EQ(c.n_anchored, 1u);
    ASSERT_EQ(c.depth.size(), kS.size());
    EXPECT_EQ(c.depth[0], 1u);  // single read covers every column once
}

TEST(ClusterConsensus, ProbeSmoothsAnErrorButLayoutKeepsAllReads) {
    std::string mut = kS;
    mut[25] = (mut[25] == 'A') ? 'T' : 'A';  // single minority error at pos 25
    std::vector<std::string> reads = {kS, kS, mut};  // 2 correct, 1 mutated
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.probe, kS);     // probe (transient) takes the majority base
    EXPECT_EQ(c.n_anchored, 3u);
    EXPECT_EQ(c.members.size(), 3u);  // ALL three reads kept in the layout
    ASSERT_EQ(c.depth.size(), kS.size());
    EXPECT_EQ(c.depth[25], 3u);  // depth preserved: 3 reads cover the variant column
}

TEST(ClusterConsensus, SubstringAnchoredAtOffset) {
    const std::string sub = kS.substr(10, 35);  // offset 10
    std::vector<std::string> reads = {kS, sub};
    const auto c = AssembleConsensus(reads, Cfg());
    EXPECT_EQ(c.probe.size(), kS.size());  // backbone spans the contig
    ASSERT_EQ(c.members.size(), 2u);
    EXPECT_TRUE(c.members[1].anchored);
    EXPECT_EQ(c.members[1].offset, 10);
}

TEST(ClusterConsensus, ReverseStrandMemberOrientedAndTracked) {
    const std::string rc = RevComp(kS);
    std::vector<std::string> reads = {kS, rc};
    const std::vector<std::uint8_t> rev = {0, 1};  // read 1 is reverse strand
    const auto c = AssembleConsensus(reads, Cfg(), rev);
    ASSERT_EQ(c.members.size(), 2u);
    EXPECT_FALSE(c.members[0].reverse);
    EXPECT_TRUE(c.members[1].reverse);    // strand recorded, never dropped
    EXPECT_TRUE(c.members[1].anchored);   // revcomp'd into frame → anchors to kS
    EXPECT_EQ(c.members[1].offset, 0);
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
    EXPECT_TRUE(c.probe.empty());
    EXPECT_TRUE(c.members.empty());
    EXPECT_TRUE(c.depth.empty());
}

}  // namespace
