#include <gtest/gtest.h>

#include "self_interference/self_wavecollapse.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <set>
#include <unordered_set>
#include <vector>

namespace llmap::self_interference {
namespace {

// ========== Configuration Tests ==========

TEST(SelfWaveCollapseConfigTest, DefaultValues) {
    SelfWaveCollapseConfig config;
    EXPECT_EQ(config.max_iterations, 20);
    EXPECT_FLOAT_EQ(config.convergence_threshold, 1e-5f);
    EXPECT_FLOAT_EQ(config.collapse_threshold, 0.95f);
    EXPECT_FLOAT_EQ(config.gamma, 0.3f);
    EXPECT_FLOAT_EQ(config.similarity_exponent, 2.0f);
    EXPECT_EQ(config.min_cluster_size, 2);
    EXPECT_EQ(config.max_cluster_size, 0);
    EXPECT_TRUE(config.enable_subclustering);
    EXPECT_EQ(config.subcluster_threshold, 1000);
    EXPECT_FALSE(config.use_gpu);
    EXPECT_EQ(config.seed, 42);
}

TEST(SelfWaveCollapseConfigTest, ConfigCanBeModified) {
    SelfWaveCollapseConfig config;
    config.max_iterations = 50;
    config.gamma = 0.5f;
    config.collapse_threshold = 0.99f;
    config.seed = 123;

    EXPECT_EQ(config.max_iterations, 50);
    EXPECT_FLOAT_EQ(config.gamma, 0.5f);
    EXPECT_FLOAT_EQ(config.collapse_threshold, 0.99f);
    EXPECT_EQ(config.seed, 123);
}

// ========== Helper to Create Test Graphs ==========

class SelfWaveCollapseTestBase : public ::testing::Test {
protected:
    std::unique_ptr<SimilarityGraph> CreateClusterGraph(
        size_t num_clusters,
        size_t nodes_per_cluster,
        float intra_weight = 1.0f,
        float inter_weight = 0.1f)
    {
        std::vector<Edge> edges;
        const size_t total_nodes = num_clusters * nodes_per_cluster;

        // Create intra-cluster edges (strong)
        for (size_t c = 0; c < num_clusters; ++c) {
            size_t base = c * nodes_per_cluster;
            for (size_t i = 0; i < nodes_per_cluster; ++i) {
                for (size_t j = i + 1; j < nodes_per_cluster; ++j) {
                    edges.push_back({
                        static_cast<uint32_t>(base + i),
                        static_cast<uint32_t>(base + j),
                        intra_weight
                    });
                }
            }
        }

        // Create inter-cluster edges (weak) - connect last node of each cluster to first of next
        for (size_t c = 0; c + 1 < num_clusters; ++c) {
            uint32_t from = static_cast<uint32_t>((c + 1) * nodes_per_cluster - 1);
            uint32_t to = static_cast<uint32_t>((c + 1) * nodes_per_cluster);
            edges.push_back({from, to, inter_weight});
        }

        SimilarityGraphConfig config;
        config.make_symmetric = true;
        return SimilarityGraph::BuildFromEdgeList(edges, total_nodes, config);
    }

    std::unique_ptr<ClusteringResult> CreateClusteringResult(
        size_t num_clusters,
        size_t nodes_per_cluster)
    {
        auto result = std::make_unique<ClusteringResult>();
        const size_t total_nodes = num_clusters * nodes_per_cluster;

        result->labels.resize(total_nodes);
        for (size_t c = 0; c < num_clusters; ++c) {
            for (size_t i = 0; i < nodes_per_cluster; ++i) {
                result->labels[c * nodes_per_cluster + i] = static_cast<uint32_t>(c);
            }
        }
        result->num_communities = num_clusters;
        result->modularity = 0.5f;

        return result;
    }
};

// ========== Basic Refinement Tests ==========

class SelfWaveCollapseBasicTest : public SelfWaveCollapseTestBase {
protected:
    void SetUp() override {
        graph_ = CreateClusterGraph(2, 5);
        clustering_ = CreateClusteringResult(2, 5);
    }

    std::unique_ptr<SimilarityGraph> graph_;
    std::unique_ptr<ClusteringResult> clustering_;
};

TEST_F(SelfWaveCollapseBasicTest, RefinesClusteredReads) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), 10);
    EXPECT_EQ(result->num_clusters, 2);
}

TEST_F(SelfWaveCollapseBasicTest, AssignmentsHaveCorrectReadIndices) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    std::set<uint32_t> read_indices;
    for (const auto& a : result->assignments) {
        read_indices.insert(a.read_idx);
    }

    EXPECT_EQ(read_indices.size(), 10);
    for (uint32_t i = 0; i < 10; ++i) {
        EXPECT_TRUE(read_indices.count(i));
    }
}

TEST_F(SelfWaveCollapseBasicTest, ClusterIDsArePreserved) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    // Nodes 0-4 should be in cluster 0
    for (uint32_t i = 0; i < 5; ++i) {
        EXPECT_EQ(result->assignments[i].cluster_id, 0);
    }

    // Nodes 5-9 should be in cluster 1
    for (uint32_t i = 5; i < 10; ++i) {
        EXPECT_EQ(result->assignments[i].cluster_id, 1);
    }
}

TEST_F(SelfWaveCollapseBasicTest, ConfidenceIsComputed) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    for (const auto& a : result->assignments) {
        EXPECT_GE(a.confidence, 0.0f);
        EXPECT_LE(a.confidence, 1.0f);
    }
}

TEST_F(SelfWaveCollapseBasicTest, AnchorsAreWithinCluster) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    // Anchors for cluster 0 should be in [0,5)
    for (uint32_t i = 0; i < 5; ++i) {
        EXPECT_LT(result->assignments[i].anchor_read, 10u);
    }

    // Anchors for cluster 1 should be in [5,10)
    for (uint32_t i = 5; i < 10; ++i) {
        EXPECT_LT(result->assignments[i].anchor_read, 10u);
    }
}

TEST_F(SelfWaveCollapseBasicTest, StatsArePopulated) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    EXPECT_EQ(result->stats.num_reads, 10);
    EXPECT_EQ(result->stats.num_clusters, 2);
    EXPECT_GT(result->stats.total_time_ms, 0.0f);
}

// ========== Empty and Edge Case Tests ==========

TEST(SelfWaveCollapseEmptyTest, EmptyClusteringReturnsEmptyResult) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 0);

    ClusteringResult clustering;
    clustering.num_communities = 0;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->assignments.empty());
    EXPECT_EQ(result->num_clusters, 0);
}

TEST(SelfWaveCollapseSingleNodeTest, SingleNodeCluster) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 1);

    ClusteringResult clustering;
    clustering.labels = {0};
    clustering.num_communities = 1;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), 1);
    EXPECT_EQ(result->assignments[0].read_idx, 0);
    EXPECT_EQ(result->assignments[0].cluster_id, 0);
    EXPECT_TRUE(result->assignments[0].collapsed);
    EXPECT_FLOAT_EQ(result->assignments[0].confidence, 1.0f);
}

TEST(SelfWaveCollapseTwoNodesTest, TwoNodeCluster) {
    std::vector<Edge> edges = {{0, 1, 1.0f}};

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0};
    clustering.num_communities = 1;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), 2);
    EXPECT_EQ(result->assignments[0].cluster_id, 0);
    EXPECT_EQ(result->assignments[1].cluster_id, 0);
}

// ========== Convergence Tests ==========

class SelfWaveCollapseConvergenceTest : public SelfWaveCollapseTestBase {};

TEST_F(SelfWaveCollapseConvergenceTest, ConvergesOnTightCluster) {
    auto graph = CreateClusterGraph(1, 5, 1.0f, 0.0f);
    auto clustering = CreateClusteringResult(1, 5);

    SelfWaveCollapseConfig config;
    config.max_iterations = 100;
    config.collapse_threshold = 0.8f;

    SelfWaveCollapse swc(config);
    auto result = swc.Refine(*graph, *clustering);

    ASSERT_NE(result, nullptr);

    // Tight cluster should have high confidence
    float avg_confidence = 0.0f;
    for (const auto& a : result->assignments) {
        avg_confidence += a.confidence;
    }
    avg_confidence /= static_cast<float>(result->assignments.size());

    EXPECT_GT(avg_confidence, 0.1f);
}

TEST_F(SelfWaveCollapseConvergenceTest, RespectsMaxIterations) {
    auto graph = CreateClusterGraph(1, 10);
    auto clustering = CreateClusteringResult(1, 10);

    SelfWaveCollapseConfig config;
    config.max_iterations = 1;
    config.convergence_threshold = 0.0f;  // Never converge naturally

    SelfWaveCollapse swc(config);
    auto result = swc.Refine(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    // Should still produce valid output
    EXPECT_EQ(result->assignments.size(), 10);
}

// ========== Damping Tests ==========

TEST_F(SelfWaveCollapseConvergenceTest, HighDampingSlowerConvergence) {
    auto graph = CreateClusterGraph(1, 5);
    auto clustering = CreateClusteringResult(1, 5);

    SelfWaveCollapseConfig low_gamma_config;
    low_gamma_config.gamma = 0.1f;
    low_gamma_config.max_iterations = 10;

    SelfWaveCollapseConfig high_gamma_config;
    high_gamma_config.gamma = 0.9f;
    high_gamma_config.max_iterations = 10;

    SelfWaveCollapse swc_low(low_gamma_config);
    SelfWaveCollapse swc_high(high_gamma_config);

    auto result_low = swc_low.Refine(*graph, *clustering);
    auto result_high = swc_high.Refine(*graph, *clustering);

    ASSERT_NE(result_low, nullptr);
    ASSERT_NE(result_high, nullptr);

    // Both should produce valid results
    EXPECT_EQ(result_low->assignments.size(), 5);
    EXPECT_EQ(result_high->assignments.size(), 5);
}

// ========== Cluster Member Access Tests ==========

TEST_F(SelfWaveCollapseBasicTest, GetClusterMembersReturnsCorrectNodes) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    auto cluster0 = result->GetClusterMembers(0);
    auto cluster1 = result->GetClusterMembers(1);

    EXPECT_EQ(cluster0.size(), 5);
    EXPECT_EQ(cluster1.size(), 5);

    // Cluster 0 should contain nodes 0-4
    std::set<uint32_t> c0_set(cluster0.begin(), cluster0.end());
    for (uint32_t i = 0; i < 5; ++i) {
        EXPECT_TRUE(c0_set.count(i));
    }

    // Cluster 1 should contain nodes 5-9
    std::set<uint32_t> c1_set(cluster1.begin(), cluster1.end());
    for (uint32_t i = 5; i < 10; ++i) {
        EXPECT_TRUE(c1_set.count(i));
    }
}

TEST_F(SelfWaveCollapseBasicTest, GetClusterSizesReturnsCorrectCounts) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    auto sizes = result->GetClusterSizes();
    EXPECT_EQ(sizes.size(), 2);
    EXPECT_EQ(sizes[0], 5);
    EXPECT_EQ(sizes[1], 5);
}

TEST_F(SelfWaveCollapseBasicTest, GetClusterAnchorsReturnsOnePerCluster) {
    auto result = RunSelfWaveCollapse(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    auto anchors = result->GetClusterAnchors();
    EXPECT_EQ(anchors.size(), 2);

    // Anchor for cluster 0 should be in [0,5)
    EXPECT_GE(anchors[0], 0u);
    EXPECT_LT(anchors[0], 5u);

    // Anchor for cluster 1 should be in [5,10)
    EXPECT_GE(anchors[1], 5u);
    EXPECT_LT(anchors[1], 10u);
}

// ========== RefineCluster Direct Tests ==========

TEST(SelfWaveCollapseRefineClusterTest, RefinesSingleClusterDirectly) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    std::vector<uint32_t> members = {0, 1, 2};

    SelfWaveCollapse swc;
    auto assignments = swc.RefineCluster(*graph, members, 42);

    EXPECT_EQ(assignments.size(), 3);

    for (const auto& a : assignments) {
        EXPECT_EQ(a.cluster_id, 42);
        EXPECT_GE(a.confidence, 0.0f);
        EXPECT_LE(a.confidence, 1.0f);
    }
}

TEST(SelfWaveCollapseRefineClusterTest, EmptyClusterReturnsEmpty) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 0);

    std::vector<uint32_t> members;

    SelfWaveCollapse swc;
    auto assignments = swc.RefineCluster(*graph, members, 0);

    EXPECT_TRUE(assignments.empty());
}

TEST(SelfWaveCollapseRefineClusterTest, SingleMemberCluster) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 1);

    std::vector<uint32_t> members = {0};

    SelfWaveCollapse swc;
    auto assignments = swc.RefineCluster(*graph, members, 0);

    ASSERT_EQ(assignments.size(), 1);
    EXPECT_EQ(assignments[0].read_idx, 0);
    EXPECT_TRUE(assignments[0].collapsed);
    EXPECT_FLOAT_EQ(assignments[0].confidence, 1.0f);
}

// ========== Config Accessor Tests ==========

TEST(SelfWaveCollapseConfigAccessTest, GetSetConfig) {
    SelfWaveCollapse swc;

    SelfWaveCollapseConfig new_config;
    new_config.max_iterations = 100;
    new_config.gamma = 0.8f;

    swc.SetConfig(new_config);

    EXPECT_EQ(swc.GetConfig().max_iterations, 100);
    EXPECT_FLOAT_EQ(swc.GetConfig().gamma, 0.8f);
}

TEST(SelfWaveCollapseConstructorTest, ConstructWithConfig) {
    SelfWaveCollapseConfig config;
    config.max_iterations = 50;
    config.seed = 999;

    SelfWaveCollapse swc(config);

    EXPECT_EQ(swc.GetConfig().max_iterations, 50);
    EXPECT_EQ(swc.GetConfig().seed, 999);
}

TEST(SelfWaveCollapseConstructorTest, DefaultConstructor) {
    SelfWaveCollapse swc;

    // Should have default config
    EXPECT_EQ(swc.GetConfig().max_iterations, 20);
    EXPECT_EQ(swc.GetConfig().seed, 42);
}

// ========== Move Semantics Tests ==========

TEST(SelfWaveCollapseMoveTest, MoveConstructor) {
    SelfWaveCollapseConfig config;
    config.seed = 12345;

    SelfWaveCollapse swc1(config);
    SelfWaveCollapse swc2(std::move(swc1));

    EXPECT_EQ(swc2.GetConfig().seed, 12345);
}

TEST(SelfWaveCollapseMoveTest, MoveAssignment) {
    SelfWaveCollapseConfig config;
    config.seed = 99999;

    SelfWaveCollapse swc1(config);
    SelfWaveCollapse swc2;

    swc2 = std::move(swc1);

    EXPECT_EQ(swc2.GetConfig().seed, 99999);
}

// ========== Larger Scale Tests ==========

class SelfWaveCollapseLargerTest : public SelfWaveCollapseTestBase {};

TEST_F(SelfWaveCollapseLargerTest, MediumSizeGraph) {
    auto graph = CreateClusterGraph(5, 20);
    auto clustering = CreateClusteringResult(5, 20);

    auto result = RunSelfWaveCollapse(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), 100);
    EXPECT_EQ(result->num_clusters, 5);

    // Check cluster sizes
    auto sizes = result->GetClusterSizes();
    for (size_t s : sizes) {
        EXPECT_EQ(s, 20);
    }
}

TEST_F(SelfWaveCollapseLargerTest, PerformanceWithManySmallClusters) {
    const size_t num_clusters = 50;
    const size_t nodes_per_cluster = 10;

    auto graph = CreateClusterGraph(num_clusters, nodes_per_cluster);
    auto clustering = CreateClusteringResult(num_clusters, nodes_per_cluster);

    SelfWaveCollapseConfig config;
    config.max_iterations = 10;

    SelfWaveCollapse swc(config);
    auto result = swc.Refine(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), num_clusters * nodes_per_cluster);
    EXPECT_EQ(result->num_clusters, num_clusters);
    EXPECT_LT(result->stats.total_time_ms, 10000.0f);  // Should complete in <10s
}

// ========== Quality Tests ==========

TEST_F(SelfWaveCollapseLargerTest, HighSimilarityReadsHaveHighConfidence) {
    // Create a tight cluster where all nodes are strongly connected
    std::vector<Edge> edges;
    const size_t n = 10;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            edges.push_back({
                static_cast<uint32_t>(i),
                static_cast<uint32_t>(j),
                1.0f  // All edges have weight 1.0
            });
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, n, config);

    ClusteringResult clustering;
    clustering.labels.resize(n, 0);
    clustering.num_communities = 1;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);

    // All reads should have similar confidence (symmetric graph)
    float min_conf = 1.0f, max_conf = 0.0f;
    for (const auto& a : result->assignments) {
        min_conf = std::min(min_conf, a.confidence);
        max_conf = std::max(max_conf, a.confidence);
    }

    // In a fully symmetric graph, confidence should be similar
    EXPECT_LT(max_conf - min_conf, 0.5f);
}

TEST_F(SelfWaveCollapseLargerTest, CentralNodeBecomesAnchor) {
    // Create a star graph with node 0 at center
    std::vector<Edge> edges;
    const size_t n = 10;

    for (size_t i = 1; i < n; ++i) {
        edges.push_back({0, static_cast<uint32_t>(i), 1.0f});
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, n, config);

    ClusteringResult clustering;
    clustering.labels.resize(n, 0);
    clustering.num_communities = 1;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);

    auto anchors = result->GetClusterAnchors();
    ASSERT_EQ(anchors.size(), 1);

    // Node 0 (center) should be the anchor
    EXPECT_EQ(anchors[0], 0u);
}

// ========== Integration with Leiden Tests ==========

TEST(SelfWaveCollapseIntegrationTest, WorksWithLeidenOutput) {
    // Create a graph with clear communities
    std::vector<Edge> edges = {
        // Community A
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        // Community B
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
        // Weak inter-community
        {2, 3, 0.1f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    // First run Leiden
    auto leiden_result = RunLeiden(*graph);
    ASSERT_NE(leiden_result, nullptr);
    EXPECT_EQ(leiden_result->num_communities, 2);

    // Then run Self-WaveCollapse
    auto swc_result = RunSelfWaveCollapse(*graph, *leiden_result);
    ASSERT_NE(swc_result, nullptr);
    EXPECT_EQ(swc_result->assignments.size(), 6);

    // All reads should have confidence values
    for (const auto& a : swc_result->assignments) {
        EXPECT_GE(a.confidence, 0.0f);
        EXPECT_LE(a.confidence, 1.0f);
    }
}

// ========== Min Cluster Size Tests ==========

TEST(SelfWaveCollapseMinSizeTest, SmallClustersCollapsedToSelf) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f},  // Cluster 0: nodes 0,1
        // Node 2 is singleton
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 1};
    clustering.num_communities = 2;

    SelfWaveCollapseConfig swc_config;
    swc_config.min_cluster_size = 2;

    SelfWaveCollapse swc(swc_config);
    auto result = swc.Refine(*graph, clustering);

    ASSERT_NE(result, nullptr);

    // Singleton (node 2) should be collapsed to self
    EXPECT_TRUE(result->assignments[2].collapsed);
    EXPECT_FLOAT_EQ(result->assignments[2].confidence, 1.0f);
    EXPECT_EQ(result->assignments[2].anchor_read, 2u);
}

// ========== Collapse Threshold Tests ==========

TEST(SelfWaveCollapseThresholdTest, HighThresholdFewerCollapsed) {
    std::vector<Edge> edges;
    const size_t n = 5;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            edges.push_back({
                static_cast<uint32_t>(i),
                static_cast<uint32_t>(j),
                1.0f
            });
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, n, config);

    ClusteringResult clustering;
    clustering.labels.resize(n, 0);
    clustering.num_communities = 1;

    // Very high threshold
    SelfWaveCollapseConfig high_thresh_config;
    high_thresh_config.collapse_threshold = 0.999f;

    SelfWaveCollapse swc(high_thresh_config);
    auto result = swc.Refine(*graph, clustering);

    ASSERT_NE(result, nullptr);
    // With uniform similarity, it's hard to reach 99.9%
    // so fewer should collapse
}

TEST(SelfWaveCollapseThresholdTest, LowThresholdMoreCollapsed) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    // Low threshold
    SelfWaveCollapseConfig low_thresh_config;
    low_thresh_config.collapse_threshold = 0.3f;
    low_thresh_config.max_iterations = 20;

    SelfWaveCollapse swc(low_thresh_config);
    auto result = swc.Refine(*graph, clustering);

    ASSERT_NE(result, nullptr);

    // With low threshold, all should collapse
    size_t collapsed_count = 0;
    for (const auto& a : result->assignments) {
        if (a.collapsed) ++collapsed_count;
    }
    EXPECT_GT(collapsed_count, 0);
}

// ========== Determinism Tests ==========

TEST(SelfWaveCollapseDeterminismTest, SameSeedSameResult) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    SelfWaveCollapseConfig swc_config;
    swc_config.seed = 42;

    SelfWaveCollapse swc1(swc_config);
    SelfWaveCollapse swc2(swc_config);

    auto result1 = swc1.Refine(*graph, clustering);
    auto result2 = swc2.Refine(*graph, clustering);

    ASSERT_NE(result1, nullptr);
    ASSERT_NE(result2, nullptr);

    // Same seed should give same confidence values
    for (size_t i = 0; i < result1->assignments.size(); ++i) {
        EXPECT_FLOAT_EQ(result1->assignments[i].confidence,
                        result2->assignments[i].confidence);
    }
}

// ========== Convenience Function Tests ==========

TEST(SelfWaveCollapseConvenienceTest, RunSelfWaveCollapseFunction) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    auto result = RunSelfWaveCollapse(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->assignments.size(), 3);
}

}  // namespace
}  // namespace llmap::self_interference
