// LLmap — assemble_pipeline (cluster→assemble orchestrator) tests.
//
// Focus: clusters formed correctly AND member.read_idx remapped local→GLOBAL,
// strand threaded through, nothing dropped (lossless, Line A).

#include "mapper/assemble_pipeline.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

namespace {

using llmap::mapper::AssembleClusters;
using llmap::mapper::ConsensusConfig;
using llmap::mapper::MinimizerClusterConfig;

std::string Backbone(std::size_t len, std::uint64_t seed) {
    static const char kB[] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(len);
    std::uint64_t x = seed;
    for (std::size_t i = 0; i < len; ++i) {
        x = x * 6364136223846793005ULL + 1442695040888963407ULL;
        s.push_back(kB[(x >> 33) & 3]);
    }
    return s;
}

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

MinimizerClusterConfig CCfg() {
    return MinimizerClusterConfig{/*k=*/15, /*w=*/5, /*min_shared=*/3,
                                  /*min_overlap_bp=*/30, /*max_freq_frac=*/0.6};
}
ConsensusConfig ACfg() { return ConsensusConfig{/*k=*/15, /*min_shared=*/3}; }

bool HasReadIdx(const llmap::mapper::ConsensusContig& c, std::size_t idx) {
    return std::any_of(c.members.begin(), c.members.end(),
                       [&](const auto& m) { return m.read_idx == idx; });
}

TEST(AssemblePipeline, Empty) {
    EXPECT_TRUE(AssembleClusters({}, CCfg(), ACfg()).empty());
}

TEST(AssemblePipeline, OverlapClusteredUnrelatedSeparate_GlobalIdx) {
    const std::string bb = Backbone(400, 7);
    // read 0 and read 2 overlap; read 1 is unrelated (sits between them).
    std::vector<std::string> reads = {bb.substr(0, 200), Backbone(200, 999999),
                                      bb.substr(100, 200)};
    const auto contigs = AssembleClusters(reads, CCfg(), ACfg());
    ASSERT_EQ(contigs.size(), 2u);

    // Find the 2-member contig (the overlap cluster) and the singleton.
    const auto* pair = contigs[0].members.size() == 2 ? &contigs[0] : &contigs[1];
    const auto* solo = contigs[0].members.size() == 2 ? &contigs[1] : &contigs[0];
    ASSERT_EQ(pair->members.size(), 2u);
    ASSERT_EQ(solo->members.size(), 1u);

    // GLOBAL indices: overlap cluster carries reads 0 and 2; singleton carries 1.
    EXPECT_TRUE(HasReadIdx(*pair, 0));
    EXPECT_TRUE(HasReadIdx(*pair, 2));
    EXPECT_EQ(solo->members[0].read_idx, 1u);

    // Lossless: every input read appears in exactly one contig member.
    std::size_t total = 0;
    for (const auto& c : contigs) total += c.members.size();
    EXPECT_EQ(total, reads.size());
}

TEST(AssemblePipeline, ReverseStrandMemberThreadedAndGlobalIdx) {
    const std::string bb = Backbone(400, 11);
    std::vector<std::string> reads = {bb.substr(0, 200),
                                      RevComp(bb.substr(50, 200))};  // reverse strand
    const auto contigs = AssembleClusters(reads, CCfg(), ACfg());
    ASSERT_EQ(contigs.size(), 1u);  // overlap → one cluster
    const auto& c = contigs[0];
    ASSERT_EQ(c.members.size(), 2u);
    EXPECT_TRUE(HasReadIdx(c, 0));
    EXPECT_TRUE(HasReadIdx(c, 1));
    // One member forward, one reverse — strand tracked, not dropped.
    EXPECT_NE(c.members[0].reverse, c.members[1].reverse);
}

}  // namespace
