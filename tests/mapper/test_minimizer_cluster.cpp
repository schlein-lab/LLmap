// LLmap — minimizer_cluster (CPU overlap clustering) tests.

#include "mapper/minimizer_cluster.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace {

using llmap::mapper::ClusterByMinimizers;
using llmap::mapper::GroupClusters;
using llmap::mapper::MinimizerClusterConfig;

// Deterministic non-repetitive backbone (LCG over ACGT) — long enough for many
// minimizers so overlaps clear the shared/span floors.
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

// Small thresholds so ~100 bp synthetic overlaps qualify.
MinimizerClusterConfig Cfg() {
    return MinimizerClusterConfig{/*k=*/15, /*w=*/5, /*min_shared=*/3,
                                  /*min_overlap_bp=*/30, /*max_freq_frac=*/0.6};
}

TEST(MinimizerCluster, Empty) {
    const auto a = ClusterByMinimizers({}, Cfg());
    EXPECT_TRUE(a.empty());
}

TEST(MinimizerCluster, SingletonOwnCluster) {
    std::vector<std::string> reads = {Backbone(200, 1)};
    const auto a = ClusterByMinimizers(reads, Cfg());
    ASSERT_EQ(a.size(), 1u);
    EXPECT_FALSE(a[0].reverse);
}

TEST(MinimizerCluster, OverlappingReadsSameCluster) {
    const std::string bb = Backbone(400, 7);
    std::vector<std::string> reads = {bb.substr(0, 200), bb.substr(100, 200)};
    const auto a = ClusterByMinimizers(reads, Cfg());
    ASSERT_EQ(a.size(), 2u);
    EXPECT_EQ(a[0].cluster_id, a[1].cluster_id);  // 100 bp overlap → merged
    EXPECT_EQ(a[0].reverse, a[1].reverse);        // same strand
}

TEST(MinimizerCluster, RevCompOverlapMergedWithStrandTracked) {
    const std::string bb = Backbone(400, 11);
    std::vector<std::string> reads = {
        bb.substr(0, 200),            // forward
        bb.substr(100, 200),          // forward, overlaps r0
        RevComp(bb.substr(150, 200))  // reverse strand of the same locus
    };
    const auto a = ClusterByMinimizers(reads, Cfg());
    ASSERT_EQ(a.size(), 3u);
    EXPECT_EQ(a[0].cluster_id, a[1].cluster_id);
    EXPECT_EQ(a[0].cluster_id, a[2].cluster_id);  // merged across strand
    EXPECT_EQ(a[0].reverse, a[1].reverse);        // r0,r1 same orientation
    EXPECT_NE(a[0].reverse, a[2].reverse);        // r2 flipped, but tracked
}

TEST(MinimizerCluster, UnrelatedReadsSeparate) {
    std::vector<std::string> reads = {Backbone(200, 3), Backbone(200, 999999)};
    const auto a = ClusterByMinimizers(reads, Cfg());
    ASSERT_EQ(a.size(), 2u);
    EXPECT_NE(a[0].cluster_id, a[1].cluster_id);
}

TEST(MinimizerCluster, GroupClustersPartitions) {
    const std::string bb = Backbone(400, 5);
    std::vector<std::string> reads = {bb.substr(0, 200), bb.substr(100, 200),
                                      Backbone(200, 424242)};
    const auto a = ClusterByMinimizers(reads, Cfg());
    const auto g = GroupClusters(a);
    std::size_t total = 0;
    for (const auto& c : g) total += c.size();
    EXPECT_EQ(total, 3u);  // every read assigned exactly once
}

}  // namespace
