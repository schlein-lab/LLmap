// Unit tests for consensus-contig placement propagation (assemble-then-map).

#include "mapper/consensus_propagation.h"
#include "mapper/cluster_consensus.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace llmap::mapper {
namespace {

ConsensusMember Mem(std::size_t idx, std::int64_t offset, bool anchored, bool rev) {
    ConsensusMember m;
    m.read_idx = idx;
    m.offset = offset;
    m.anchored = anchored;
    m.reverse = rev;
    return m;
}

TEST(ConsensusPropagation, ForwardProbePropagatesOffsetAndStrand) {
    // Two members of a length-100 contig; probe maps forward at ref 1000.
    std::vector<std::string> reads = {std::string(100, 'A'), std::string(40, 'C')};
    ConsensusContig c;
    c.length = 100;
    c.members = {Mem(0, 0, true, false), Mem(1, 30, true, true)};

    const auto pl = PropagatePlacement(c, reads, {"chr1", 1000, false});
    ASSERT_EQ(pl.size(), 2u);
    // member 0: ref_start = 1000 + 0; forward
    EXPECT_EQ(pl[0].ref_start, 1000);
    EXPECT_FALSE(pl[0].reverse);
    EXPECT_TRUE(pl[0].anchored);
    EXPECT_EQ(pl[0].ref_name, "chr1");
    // member 1 at offset 30, reverse-in-contig
    EXPECT_EQ(pl[1].ref_start, 1030);
    EXPECT_TRUE(pl[1].reverse);     // forward probe ⇒ effective strand = m.reverse
}

TEST(ConsensusPropagation, ReverseProbeFlipsCoordinatesAndStrand) {
    // Contig length 100; probe maps REVERSE at ref 1000. A forward member of
    // length 40 at offset 30 spans contig [30,70) → ref_start = 1000 + 100 - 70 = 1030.
    std::vector<std::string> reads = {std::string(40, 'C')};
    ConsensusContig c;
    c.length = 100;
    c.members = {Mem(0, 30, true, false)};

    const auto pl = PropagatePlacement(c, reads, {"chr1", 1000, true});
    ASSERT_EQ(pl.size(), 1u);
    EXPECT_EQ(pl[0].ref_start, 1000 + 100 - (30 + 40));   // 1030
    EXPECT_TRUE(pl[0].reverse);     // reverse probe ⇒ effective = !m.reverse = !false
}

TEST(ConsensusPropagation, ReverseProbeWithReverseMemberIsForwardOnRef) {
    std::vector<std::string> reads = {std::string(40, 'C')};
    ConsensusContig c;
    c.length = 100;
    c.members = {Mem(0, 30, true, true)};
    const auto pl = PropagatePlacement(c, reads, {"chr1", 1000, true});
    ASSERT_EQ(pl.size(), 1u);
    EXPECT_FALSE(pl[0].reverse);    // !m.reverse = !true = false
}

TEST(ConsensusPropagation, UnanchoredMemberRoutedSolo) {
    std::vector<std::string> reads = {std::string(100, 'A'), std::string(50, 'G')};
    ConsensusContig c;
    c.length = 100;
    c.members = {Mem(0, 0, true, false), Mem(1, 0, false, false)};
    const auto pl = PropagatePlacement(c, reads, {"chr1", 1000, false});
    ASSERT_EQ(pl.size(), 2u);
    EXPECT_TRUE(pl[0].anchored);
    EXPECT_FALSE(pl[1].anchored);   // caller maps read 1 solo (lossless fallback)
}

TEST(ConsensusPropagation, EmptyContigEmptyResult) {
    std::vector<std::string> reads;
    ConsensusContig c;
    c.length = 0;
    const auto pl = PropagatePlacement(c, reads, {"chr1", 1000, false});
    EXPECT_TRUE(pl.empty());
}

}  // namespace
}  // namespace llmap::mapper
