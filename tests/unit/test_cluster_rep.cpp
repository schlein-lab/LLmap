#include <gtest/gtest.h>

#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
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

TEST(ClusterRepConfigTest, DefaultValues) {
    ClusterRepConfig config;
    EXPECT_EQ(config.method, ClusterRepConfig::Method::Medoid);
    EXPECT_EQ(config.max_medoid_candidates, 0);
    EXPECT_FALSE(config.use_approximate_medoid);
    EXPECT_EQ(config.approx_sample_size, 100);
    EXPECT_EQ(config.seed, 42);
    EXPECT_TRUE(config.use_confidence_tiebreaker);
    EXPECT_EQ(config.min_cluster_size, 1);
    EXPECT_TRUE(config.parallel);
    EXPECT_EQ(config.parallel_threshold, 1000);
}

TEST(ClusterRepConfigTest, ConfigCanBeModified) {
    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxDegree;
    config.max_medoid_candidates = 50;
    config.use_approximate_medoid = true;
    config.seed = 123;

    EXPECT_EQ(config.method, ClusterRepConfig::Method::MaxDegree);
    EXPECT_EQ(config.max_medoid_candidates, 50);
    EXPECT_TRUE(config.use_approximate_medoid);
    EXPECT_EQ(config.seed, 123);
}

// ========== Helper to Create Test Graphs ==========

class ClusterRepTestBase : public ::testing::Test {
protected:
    std::unique_ptr<SimilarityGraph> CreateClusterGraph(
        size_t num_clusters,
        size_t nodes_per_cluster,
        float intra_weight = 1.0f,
        float inter_weight = 0.1f)
    {
        std::vector<Edge> edges;
        const size_t total_nodes = num_clusters * nodes_per_cluster;

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

    std::unique_ptr<SimilarityGraph> CreateStarGraph(size_t n) {
        std::vector<Edge> edges;
        for (size_t i = 1; i < n; ++i) {
            edges.push_back({0, static_cast<uint32_t>(i), 1.0f});
        }

        SimilarityGraphConfig config;
        config.make_symmetric = true;
        return SimilarityGraph::BuildFromEdgeList(edges, n, config);
    }
};

// ========== Basic Selection Tests ==========

class ClusterRepBasicTest : public ClusterRepTestBase {
protected:
    void SetUp() override {
        graph_ = CreateClusterGraph(2, 5);
        clustering_ = CreateClusteringResult(2, 5);
    }

    std::unique_ptr<SimilarityGraph> graph_;
    std::unique_ptr<ClusteringResult> clustering_;
};

TEST_F(ClusterRepBasicTest, SelectsRepresentativesForEachCluster) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 2);
    EXPECT_EQ(result->representative_reads.size(), 2);
}

TEST_F(ClusterRepBasicTest, RepresentativesAreWithinClusters) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    // Cluster 0 rep should be in [0, 5)
    auto rep0 = result->GetRepresentative(0);
    ASSERT_TRUE(rep0.has_value());
    EXPECT_GE(rep0->read_idx, 0u);
    EXPECT_LT(rep0->read_idx, 5u);

    // Cluster 1 rep should be in [5, 10)
    auto rep1 = result->GetRepresentative(1);
    ASSERT_TRUE(rep1.has_value());
    EXPECT_GE(rep1->read_idx, 5u);
    EXPECT_LT(rep1->read_idx, 10u);
}

TEST_F(ClusterRepBasicTest, RepresentativeInfoHasCorrectClusterSize) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    for (const auto& rep : result->representatives) {
        EXPECT_EQ(rep.cluster_size, 5);
    }
}

TEST_F(ClusterRepBasicTest, CentralityScoreIsComputed) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    for (const auto& rep : result->representatives) {
        EXPECT_GE(rep.centrality_score, 0.0f);
    }
}

TEST_F(ClusterRepBasicTest, StatsArePopulated) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    EXPECT_EQ(result->stats.num_clusters, 2);
    EXPECT_EQ(result->stats.num_representatives, 2);
    EXPECT_EQ(result->stats.clusters_skipped, 0);
    EXPECT_GT(result->stats.total_time_ms, 0.0f);
}

// ========== Empty and Edge Case Tests ==========

TEST(ClusterRepEmptyTest, EmptyClusteringReturnsEmptyResult) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 0);

    ClusteringResult clustering;
    clustering.num_communities = 0;

    auto result = SelectClusterRepresentatives(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_TRUE(result->representatives.empty());
    EXPECT_TRUE(result->representative_reads.empty());
}

TEST(ClusterRepSingleNodeTest, SingleNodeCluster) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 1);

    ClusteringResult clustering;
    clustering.labels = {0};
    clustering.num_communities = 1;

    auto result = SelectClusterRepresentatives(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    EXPECT_EQ(result->representatives[0].read_idx, 0);
    EXPECT_EQ(result->representatives[0].cluster_size, 1);
}

TEST(ClusterRepTwoNodesTest, TwoNodeCluster) {
    std::vector<Edge> edges = {{0, 1, 1.0f}};

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0};
    clustering.num_communities = 1;

    auto result = SelectClusterRepresentatives(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    // Either node could be the representative (symmetric graph)
    EXPECT_TRUE(result->representatives[0].read_idx == 0 ||
                result->representatives[0].read_idx == 1);
}

// ========== Star Graph Tests ==========

class ClusterRepStarGraphTest : public ClusterRepTestBase {};

TEST_F(ClusterRepStarGraphTest, CenterNodeIsSelectedAsMedoid) {
    auto graph = CreateStarGraph(10);

    ClusteringResult clustering;
    clustering.labels.resize(10, 0);
    clustering.num_communities = 1;

    auto result = SelectClusterRepresentatives(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    // Node 0 is center, should be the medoid
    EXPECT_EQ(result->representatives[0].read_idx, 0u);
}

TEST_F(ClusterRepStarGraphTest, CenterHasHighestCentrality) {
    auto graph = CreateStarGraph(10);

    ClusteringResult clustering;
    clustering.labels.resize(10, 0);
    clustering.num_communities = 1;

    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxDegree;

    ClusterRepSelector selector(config);
    auto result = selector.Select(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives[0].read_idx, 0u);
}

// ========== Method Comparison Tests ==========

class ClusterRepMethodTest : public ClusterRepTestBase {};

TEST_F(ClusterRepMethodTest, MedoidMethod) {
    auto graph = CreateClusterGraph(1, 5);
    auto clustering = CreateClusteringResult(1, 5);

    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::Medoid;

    ClusterRepSelector selector(config);
    auto result = selector.Select(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
}

TEST_F(ClusterRepMethodTest, MaxDegreeMethod) {
    auto graph = CreateClusterGraph(1, 5);
    auto clustering = CreateClusteringResult(1, 5);

    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxDegree;

    ClusterRepSelector selector(config);
    auto result = selector.Select(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
}

TEST_F(ClusterRepMethodTest, MaxConfidenceMethodFallsBackWithoutConfidence) {
    auto graph = CreateClusterGraph(1, 5);
    auto clustering = CreateClusteringResult(1, 5);

    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxConfidence;

    ClusterRepSelector selector(config);
    auto result = selector.Select(*graph, *clustering);

    // Should fall back to MaxDegree
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
}

// ========== SelectForCluster Direct Tests ==========

TEST(ClusterRepSelectForClusterTest, SelectsSingleClusterDirectly) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    std::vector<uint32_t> members = {0, 1, 2};

    ClusterRepSelector selector;
    auto rep_opt = selector.SelectForCluster(*graph, members);

    ASSERT_TRUE(rep_opt.has_value());
    EXPECT_LT(*rep_opt, 3u);
}

TEST(ClusterRepSelectForClusterTest, EmptyClusterReturnsNullopt) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 0);

    std::vector<uint32_t> members;

    ClusterRepSelector selector;
    auto rep_opt = selector.SelectForCluster(*graph, members);

    EXPECT_FALSE(rep_opt.has_value());
}

TEST(ClusterRepSelectForClusterTest, SingleMemberReturnsItself) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 1);

    std::vector<uint32_t> members = {0};

    ClusterRepSelector selector;
    auto rep_opt = selector.SelectForCluster(*graph, members);

    ASSERT_TRUE(rep_opt.has_value());
    EXPECT_EQ(*rep_opt, 0u);
}

// ========== Confidence-Based Selection Tests ==========

class ClusterRepWithConfidenceTest : public ClusterRepTestBase {};

TEST_F(ClusterRepWithConfidenceTest, UsesConfidenceFromSWCResult) {
    auto graph = CreateClusterGraph(1, 5);
    auto clustering = CreateClusteringResult(1, 5);

    // Create SWC result with node 2 having highest confidence
    SelfWaveCollapseResult swc_result;
    swc_result.num_clusters = 1;
    swc_result.assignments.resize(5);
    for (uint32_t i = 0; i < 5; ++i) {
        swc_result.assignments[i].read_idx = i;
        swc_result.assignments[i].cluster_id = 0;
        swc_result.assignments[i].confidence = 0.5f;  // Base confidence
    }
    swc_result.assignments[2].confidence = 0.99f;  // Node 2 has highest

    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxConfidence;

    ClusterRepSelector selector(config);
    auto result = selector.SelectWithConfidence(*graph, *clustering, swc_result);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    EXPECT_EQ(result->representatives[0].read_idx, 2u);
}

TEST_F(ClusterRepWithConfidenceTest, ConfidenceIsTiebreaker) {
    // Create a fully connected cluster (all edges equal weight)
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

    // All nodes have equal distances; confidence should break tie
    SelfWaveCollapseResult swc_result;
    swc_result.num_clusters = 1;
    swc_result.assignments.resize(n);
    for (uint32_t i = 0; i < n; ++i) {
        swc_result.assignments[i].read_idx = i;
        swc_result.assignments[i].cluster_id = 0;
        swc_result.assignments[i].confidence = 0.5f;
    }
    swc_result.assignments[3].confidence = 0.95f;  // Node 3 highest

    ClusterRepConfig rep_config;
    rep_config.method = ClusterRepConfig::Method::Medoid;
    rep_config.use_confidence_tiebreaker = true;

    ClusterRepSelector selector(rep_config);
    auto result = selector.SelectWithConfidence(*graph, clustering, swc_result);

    ASSERT_NE(result, nullptr);
    // In a symmetric graph, confidence tiebreaker might select node 3
    EXPECT_EQ(result->representatives.size(), 1);
}

// ========== Result Query Tests ==========

TEST_F(ClusterRepBasicTest, GetRepresentativeReturnsCorrectInfo) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    auto rep0 = result->GetRepresentative(0);
    ASSERT_TRUE(rep0.has_value());
    EXPECT_EQ(rep0->cluster_id, 0);

    auto rep1 = result->GetRepresentative(1);
    ASSERT_TRUE(rep1.has_value());
    EXPECT_EQ(rep1->cluster_id, 1);

    // Non-existent cluster
    auto rep99 = result->GetRepresentative(99);
    EXPECT_FALSE(rep99.has_value());
}

TEST_F(ClusterRepBasicTest, IsRepresentativeWorksCorrectly) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    // Representatives should return true
    for (uint32_t rep_idx : result->representative_reads) {
        EXPECT_TRUE(result->IsRepresentative(rep_idx));
    }

    // At least some non-representatives should exist
    size_t non_rep_count = 0;
    for (uint32_t i = 0; i < 10; ++i) {
        if (!result->IsRepresentative(i)) {
            ++non_rep_count;
        }
    }
    EXPECT_EQ(non_rep_count, 8);  // 10 total - 2 representatives
}

TEST_F(ClusterRepBasicTest, GetRepresentativeReadsReturnsSorted) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    auto sorted_reps = result->GetRepresentativeReads();
    EXPECT_EQ(sorted_reps.size(), 2);

    // Should be sorted
    EXPECT_LE(sorted_reps[0], sorted_reps[1]);
}

TEST_F(ClusterRepBasicTest, GetClusterForRepresentativeWorks) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    for (const auto& rep : result->representatives) {
        auto cluster_opt = result->GetClusterForRepresentative(rep.read_idx);
        ASSERT_TRUE(cluster_opt.has_value());
        EXPECT_EQ(*cluster_opt, rep.cluster_id);
    }

    // Non-representative should return nullopt
    uint32_t non_rep = 0;
    for (uint32_t i = 0; i < 10; ++i) {
        if (!result->IsRepresentative(i)) {
            non_rep = i;
            break;
        }
    }
    auto cluster_opt = result->GetClusterForRepresentative(non_rep);
    EXPECT_FALSE(cluster_opt.has_value());
}

// ========== Config Accessor Tests ==========

TEST(ClusterRepConfigAccessTest, GetSetConfig) {
    ClusterRepSelector selector;

    ClusterRepConfig new_config;
    new_config.method = ClusterRepConfig::Method::MaxDegree;
    new_config.seed = 999;

    selector.SetConfig(new_config);

    EXPECT_EQ(selector.GetConfig().method, ClusterRepConfig::Method::MaxDegree);
    EXPECT_EQ(selector.GetConfig().seed, 999);
}

TEST(ClusterRepConstructorTest, ConstructWithConfig) {
    ClusterRepConfig config;
    config.method = ClusterRepConfig::Method::MaxConfidence;
    config.use_approximate_medoid = true;

    ClusterRepSelector selector(config);

    EXPECT_EQ(selector.GetConfig().method, ClusterRepConfig::Method::MaxConfidence);
    EXPECT_TRUE(selector.GetConfig().use_approximate_medoid);
}

TEST(ClusterRepConstructorTest, DefaultConstructor) {
    ClusterRepSelector selector;

    EXPECT_EQ(selector.GetConfig().method, ClusterRepConfig::Method::Medoid);
    EXPECT_FALSE(selector.GetConfig().use_approximate_medoid);
}

// ========== Move Semantics Tests ==========

TEST(ClusterRepMoveTest, MoveConstructor) {
    ClusterRepConfig config;
    config.seed = 12345;

    ClusterRepSelector sel1(config);
    ClusterRepSelector sel2(std::move(sel1));

    EXPECT_EQ(sel2.GetConfig().seed, 12345);
}

TEST(ClusterRepMoveTest, MoveAssignment) {
    ClusterRepConfig config;
    config.seed = 99999;

    ClusterRepSelector sel1(config);
    ClusterRepSelector sel2;

    sel2 = std::move(sel1);

    EXPECT_EQ(sel2.GetConfig().seed, 99999);
}

// ========== Larger Scale Tests ==========

class ClusterRepLargerTest : public ClusterRepTestBase {};

TEST_F(ClusterRepLargerTest, MediumSizeGraph) {
    auto graph = CreateClusterGraph(10, 20);
    auto clustering = CreateClusteringResult(10, 20);

    auto result = SelectClusterRepresentatives(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 10);
    EXPECT_EQ(result->stats.num_clusters, 10);
    EXPECT_EQ(result->stats.num_representatives, 10);
}

TEST_F(ClusterRepLargerTest, PerformanceWithManyClusters) {
    const size_t num_clusters = 100;
    const size_t nodes_per_cluster = 10;

    auto graph = CreateClusterGraph(num_clusters, nodes_per_cluster);
    auto clustering = CreateClusteringResult(num_clusters, nodes_per_cluster);

    auto result = SelectClusterRepresentatives(*graph, *clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), num_clusters);
    EXPECT_LT(result->stats.total_time_ms, 10000.0f);  // Should complete in <10s
}

// ========== Min Cluster Size Tests ==========

TEST(ClusterRepMinSizeTest, SkipsSmallClusters) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f},  // Cluster 0: nodes 0,1 (size 2)
        // Node 2 is singleton cluster
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 1};
    clustering.num_communities = 2;

    ClusterRepConfig rep_config;
    rep_config.min_cluster_size = 2;

    ClusterRepSelector selector(rep_config);
    auto result = selector.Select(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    EXPECT_EQ(result->stats.clusters_skipped, 1);
}

// ========== Approximate Medoid Tests ==========

TEST_F(ClusterRepLargerTest, ApproximateMedoidForLargeClusters) {
    // Create a single large cluster
    const size_t n = 200;
    std::vector<Edge> edges;

    // Create a sparse graph
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < std::min(i + 5, n); ++j) {
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

    ClusterRepConfig rep_config;
    rep_config.use_approximate_medoid = true;
    rep_config.approx_sample_size = 50;

    ClusterRepSelector selector(rep_config);
    auto result = selector.Select(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
    EXPECT_LT(result->representatives[0].read_idx, n);
}

// ========== Determinism Tests ==========

TEST(ClusterRepDeterminismTest, SameSeedSameResult) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    ClusterRepConfig rep_config;
    rep_config.seed = 42;

    ClusterRepSelector sel1(rep_config);
    ClusterRepSelector sel2(rep_config);

    auto result1 = sel1.Select(*graph, clustering);
    auto result2 = sel2.Select(*graph, clustering);

    ASSERT_NE(result1, nullptr);
    ASSERT_NE(result2, nullptr);

    // Same seed should give same representative
    EXPECT_EQ(result1->representatives[0].read_idx,
              result2->representatives[0].read_idx);
}

// ========== Convenience Function Tests ==========

TEST(ClusterRepConvenienceTest, SelectClusterRepresentativesFunction) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    auto result = SelectClusterRepresentatives(*graph, clustering);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
}

TEST(ClusterRepConvenienceTest, SelectClusterRepresentativesWithSWCFunction) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    ClusteringResult clustering;
    clustering.labels = {0, 0, 0};
    clustering.num_communities = 1;

    SelfWaveCollapseResult swc_result;
    swc_result.num_clusters = 1;
    swc_result.assignments.resize(3);
    for (uint32_t i = 0; i < 3; ++i) {
        swc_result.assignments[i].read_idx = i;
        swc_result.assignments[i].cluster_id = 0;
        swc_result.assignments[i].confidence = 0.8f;
    }

    auto result = SelectClusterRepresentatives(*graph, clustering, swc_result);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->representatives.size(), 1);
}

// ========== Average Distance Computation Tests ==========

TEST_F(ClusterRepBasicTest, AvgDistanceIsComputed) {
    auto result = SelectClusterRepresentatives(*graph_, *clustering_);

    ASSERT_NE(result, nullptr);

    for (const auto& rep : result->representatives) {
        // With fully connected clusters (weight=1), distance = 1 - weight = 0
        EXPECT_GE(rep.avg_distance_to_members, 0.0f);
        EXPECT_LE(rep.avg_distance_to_members, 1.0f);
    }
}

// ========== Integration with Leiden Tests ==========

TEST(ClusterRepIntegrationTest, WorksWithLeidenOutput) {
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

    // Then select representatives
    auto rep_result = SelectClusterRepresentatives(*graph, *leiden_result);
    ASSERT_NE(rep_result, nullptr);
    EXPECT_EQ(rep_result->representatives.size(), 2);
}

// ========== Integration with Self-WaveCollapse Tests ==========

TEST(ClusterRepIntegrationTest, WorksWithSWCOutput) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
        {2, 3, 0.1f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    auto leiden_result = RunLeiden(*graph);
    ASSERT_NE(leiden_result, nullptr);

    auto swc_result = RunSelfWaveCollapse(*graph, *leiden_result);
    ASSERT_NE(swc_result, nullptr);

    auto rep_result = SelectClusterRepresentatives(*graph, *leiden_result, *swc_result);
    ASSERT_NE(rep_result, nullptr);
    EXPECT_EQ(rep_result->representatives.size(), 2);

    // Confidence should be from SWC
    for (const auto& rep : rep_result->representatives) {
        EXPECT_GE(rep.confidence, 0.0f);
        EXPECT_LE(rep.confidence, 1.0f);
    }
}

}  // namespace
}  // namespace llmap::self_interference
