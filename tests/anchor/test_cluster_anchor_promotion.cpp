// LLmap — ClusterAnchorPromotion tests.
//
// Synthetic Stage-1 SelfWaveCollapseResult mirror; verify policy gates,
// rank ordering, tag content, and incremental AnchorStore update.

#include "anchor/cluster_anchor_promotion.h"

#include "anchor/anchor_store.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

using namespace llmap::anchor;

namespace {

ReadClusterAssignment MakeAssignment(std::uint32_t read_idx,
                                       std::uint32_t cluster_id,
                                       float confidence,
                                       bool collapsed = true,
                                       std::uint32_t anchor_read = 0) {
    ReadClusterAssignment a;
    a.read_idx = read_idx;
    a.cluster_id = cluster_id;
    a.confidence = confidence;
    a.collapsed = collapsed;
    a.anchor_read = anchor_read;
    return a;
}

}  // namespace

TEST(ClusterAnchorPromotion, PromotesQualifyingCluster) {
    AnchorStore store;
    std::vector<std::string> reads = {
        "ACGTACGT",  // 0
        "ACGTACGT",  // 1 — same seq, cluster member
        "ACGTACGT",  // 2
        "GGGGGGGG",  // 3 — singleton, won't promote
    };
    std::vector<ReadClusterAssignment> assignments = {
        MakeAssignment(0, /*cid*/ 7, 0.9f, true, /*anchor_read*/ 0),
        MakeAssignment(1, 7, 0.85f),
        MakeAssignment(2, 7, 0.95f),
        MakeAssignment(3, /*cid*/ 8, 0.95f, true),  // singleton
    };

    ClusterPromotionConfig cfg;
    cfg.min_cluster_size = 3;
    auto stats = PromoteClusterAnchors(assignments, reads, cfg, store);

    EXPECT_EQ(stats.clusters_considered, 2u);
    EXPECT_EQ(stats.clusters_promoted,   1u);
    EXPECT_EQ(stats.clusters_skipped_too_small, 1u);  // cluster 8
    EXPECT_EQ(stats.anchors_emitted, 3u);             // capped by member count
    EXPECT_EQ(store.size(), 3u);

    // anchor_ids follow CLUSTER:7:0, CLUSTER:7:1, CLUSTER:7:2 — rank-0 is
    // the highest-confidence member (rid 2).
    auto* rank0 = store.ById("CLUSTER:7:0");
    ASSERT_NE(rank0, nullptr);
    EXPECT_EQ(rank0->source, AnchorSource::Computed_Cluster);
    EXPECT_EQ(rank0->sequence, "ACGTACGT");
    // tags include computed_cluster, cluster_size:3, cluster_conf:0.90,
    // cluster_representative on rank-0.
    bool has_repr_tag = false, has_size_tag = false;
    for (const auto& t : rank0->tags) {
        if (t == "cluster_representative")  has_repr_tag = true;
        if (t == "cluster_size:3")          has_size_tag = true;
    }
    EXPECT_TRUE(has_repr_tag);
    EXPECT_TRUE(has_size_tag);
}

TEST(ClusterAnchorPromotion, SkipsLowConfidenceClusters) {
    AnchorStore store;
    std::vector<std::string> reads(5, "ACGTACGT");
    std::vector<ReadClusterAssignment> assignments;
    for (std::uint32_t i = 0; i < 5; ++i) {
        assignments.push_back(MakeAssignment(i, /*cid*/ 1, 0.3f));
    }
    ClusterPromotionConfig cfg;
    cfg.min_cluster_size = 3;
    cfg.min_cluster_confidence = 0.7f;

    auto stats = PromoteClusterAnchors(assignments, reads, cfg, store);
    EXPECT_EQ(stats.clusters_promoted, 0u);
    EXPECT_EQ(stats.clusters_skipped_low_confidence, 1u);
    EXPECT_EQ(store.size(), 0u);
}

TEST(ClusterAnchorPromotion, RespectMaxAnchorsPerCluster) {
    AnchorStore store;
    std::vector<std::string> reads(20, "ACGTACGTACGT");
    std::vector<ReadClusterAssignment> assignments;
    for (std::uint32_t i = 0; i < 20; ++i) {
        assignments.push_back(
            MakeAssignment(i, /*cid*/ 1,
                            // Decreasing confidence so rank order is stable
                            1.0f - 0.01f * i));
    }
    ClusterPromotionConfig cfg;
    cfg.min_cluster_size = 3;
    cfg.max_anchors_per_cluster = 5;
    cfg.subsample_large_clusters = true;
    auto stats = PromoteClusterAnchors(assignments, reads, cfg, store);
    EXPECT_EQ(stats.anchors_emitted, 5u);
    EXPECT_EQ(store.size(), 5u);
}

TEST(ClusterAnchorPromotion, CustomNamespaceUsedInAnchorId) {
    AnchorStore store;
    std::vector<std::string> reads(4, "AAAAAAAA");
    std::vector<ReadClusterAssignment> assignments = {
        MakeAssignment(0, 11, 0.95f),
        MakeAssignment(1, 11, 0.90f),
        MakeAssignment(2, 11, 0.85f),
    };
    ClusterPromotionConfig cfg;
    cfg.min_cluster_size = 3;
    cfg.anchor_id_namespace = "STAGE1_CLUSTER";
    PromoteClusterAnchors(assignments, reads, cfg, store);
    EXPECT_NE(store.ById("STAGE1_CLUSTER:11:0"), nullptr);
}

TEST(ClusterAnchorPromotion, EmptyAssignmentsYieldsEmptyStats) {
    AnchorStore store;
    std::vector<std::string> reads;
    std::vector<ReadClusterAssignment> assignments;
    auto stats = PromoteClusterAnchors(assignments, reads, {}, store);
    EXPECT_EQ(stats.clusters_considered, 0u);
    EXPECT_EQ(stats.clusters_promoted, 0u);
    EXPECT_EQ(store.size(), 0u);
}

TEST(ClusterAnchorPromotion, AdditiveToExistingStore) {
    AnchorStore store;
    // Pre-populate one GENCODE-style anchor so we test the additive
    // behaviour rather than replace.
    AnchorRecord pre;
    pre.anchor_id = "GENCODE:ENST.1:exon1";
    pre.source = AnchorSource::Gencode;
    store.AddAnchor(std::move(pre));
    EXPECT_EQ(store.size(), 1u);

    std::vector<std::string> reads(3, "AAAA");
    std::vector<ReadClusterAssignment> assignments = {
        MakeAssignment(0, 1, 0.95f),
        MakeAssignment(1, 1, 0.90f),
        MakeAssignment(2, 1, 0.85f),
    };
    ClusterPromotionConfig cfg;
    cfg.min_cluster_size = 3;
    auto stats = PromoteClusterAnchors(assignments, reads, cfg, store);
    EXPECT_EQ(stats.anchors_emitted, 3u);
    EXPECT_EQ(store.size(), 4u);  // 1 pre + 3 cluster
    EXPECT_NE(store.ById("GENCODE:ENST.1:exon1"), nullptr);
    EXPECT_NE(store.ById("CLUSTER:1:0"), nullptr);
}
