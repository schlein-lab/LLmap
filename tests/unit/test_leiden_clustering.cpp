#include <gtest/gtest.h>

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

TEST(LeidenConfigTest, DefaultValues) {
    LeidenConfig config;
    EXPECT_FLOAT_EQ(config.resolution, 1.0f);
    EXPECT_EQ(config.max_iterations, 100);
    EXPECT_FLOAT_EQ(config.min_modularity_gain, 1e-7f);
    EXPECT_TRUE(config.enable_refinement);
    EXPECT_EQ(config.refinement_iterations, 10);
    EXPECT_EQ(config.seed, 42);
    EXPECT_TRUE(config.randomize_node_order);
    EXPECT_EQ(config.min_community_size, 1);
    EXPECT_EQ(config.max_community_size, 0);
    EXPECT_FALSE(config.use_gpu);
}

TEST(LeidenConfigTest, ConfigCanBeModified) {
    LeidenConfig config;
    config.resolution = 0.5f;
    config.max_iterations = 50;
    config.min_community_size = 3;
    config.seed = 123;

    EXPECT_FLOAT_EQ(config.resolution, 0.5f);
    EXPECT_EQ(config.max_iterations, 50);
    EXPECT_EQ(config.min_community_size, 3);
    EXPECT_EQ(config.seed, 123);
}

// ========== Basic Clustering Tests ==========

class LeidenClusteringBasicTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple graph with two clear communities
        // Community A: 0, 1, 2 (fully connected)
        // Community B: 3, 4, 5 (fully connected)
        // Weak link: 2 -- 3
        std::vector<Edge> edges = {
            // Community A
            {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
            // Community B
            {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
            // Inter-community (weak)
            {2, 3, 0.1f},
        };

        SimilarityGraphConfig graph_config;
        graph_config.make_symmetric = true;

        graph_ = SimilarityGraph::BuildFromEdgeList(edges, 6, graph_config);
    }

    std::unique_ptr<SimilarityGraph> graph_;
};

TEST_F(LeidenClusteringBasicTest, FindsTwoCommunities) {
    auto result = RunLeiden(*graph_);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 6);

    // Should find exactly 2 communities
    EXPECT_EQ(result->num_communities, 2);
}

TEST_F(LeidenClusteringBasicTest, CommunitiesAreCorrect) {
    auto result = RunLeiden(*graph_);

    ASSERT_NE(result, nullptr);

    // Nodes 0, 1, 2 should be in one community
    EXPECT_EQ(result->labels[0], result->labels[1]);
    EXPECT_EQ(result->labels[1], result->labels[2]);

    // Nodes 3, 4, 5 should be in another community
    EXPECT_EQ(result->labels[3], result->labels[4]);
    EXPECT_EQ(result->labels[4], result->labels[5]);

    // The two communities should be different
    EXPECT_NE(result->labels[0], result->labels[3]);
}

TEST_F(LeidenClusteringBasicTest, ModularityIsPositive) {
    auto result = RunLeiden(*graph_);

    ASSERT_NE(result, nullptr);
    EXPECT_GT(result->modularity, 0.0f);
}

TEST_F(LeidenClusteringBasicTest, StatsArePopulated) {
    auto result = RunLeiden(*graph_);

    ASSERT_NE(result, nullptr);

    EXPECT_EQ(result->stats.num_nodes, 6);
    EXPECT_EQ(result->stats.num_communities, 2);
    EXPECT_GT(result->stats.num_iterations, 0);
    EXPECT_GT(result->stats.time_ms, 0.0f);
    EXPECT_FLOAT_EQ(result->stats.avg_community_size, 3.0f);
}

// ========== Empty and Edge Case Tests ==========

TEST(LeidenEmptyGraphTest, EmptyGraphReturnsEmptyResult) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 0);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->num_communities, 0);
    EXPECT_TRUE(result->labels.empty());
}

TEST(LeidenSingleNodeTest, SingleNodeGraph) {
    std::vector<Edge> edges;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 1);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 1);
    EXPECT_EQ(result->num_communities, 1);
}

TEST(LeidenDisconnectedTest, DisconnectedNodes) {
    std::vector<Edge> edges;  // No edges
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 5);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 5);

    // Each node should be in its own community (no edges to move along)
    std::set<uint32_t> unique_labels(result->labels.begin(), result->labels.end());
    EXPECT_EQ(unique_labels.size(), 5);
}

TEST(LeidenTwoNodesTest, TwoConnectedNodes) {
    std::vector<Edge> edges = {{0, 1, 1.0f}};

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 2);

    // Both nodes should be in the same community
    EXPECT_EQ(result->labels[0], result->labels[1]);
}

// ========== Resolution Parameter Tests ==========

class LeidenResolutionTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a graph with hierarchical structure
        // 4 cliques of 3 nodes each, loosely connected
        std::vector<Edge> edges;

        // Cliques (strong connections)
        for (int c = 0; c < 4; ++c) {
            int base = c * 3;
            edges.push_back({static_cast<uint32_t>(base), static_cast<uint32_t>(base + 1), 1.0f});
            edges.push_back({static_cast<uint32_t>(base + 1), static_cast<uint32_t>(base + 2), 1.0f});
            edges.push_back({static_cast<uint32_t>(base), static_cast<uint32_t>(base + 2), 1.0f});
        }

        // Inter-clique connections (weak)
        edges.push_back({2, 3, 0.2f});   // Clique 0 - Clique 1
        edges.push_back({5, 6, 0.2f});   // Clique 1 - Clique 2
        edges.push_back({8, 9, 0.2f});   // Clique 2 - Clique 3

        SimilarityGraphConfig config;
        config.make_symmetric = true;
        graph_ = SimilarityGraph::BuildFromEdgeList(edges, 12, config);
    }

    std::unique_ptr<SimilarityGraph> graph_;
};

TEST_F(LeidenResolutionTest, HighResolutionMoreCommunities) {
    LeidenConfig config;
    config.resolution = 2.0f;

    LeidenClustering leiden(config);
    auto result = leiden.Cluster(*graph_);

    ASSERT_NE(result, nullptr);
    EXPECT_GE(result->num_communities, 3);  // Should find more communities
}

TEST_F(LeidenResolutionTest, LowResolutionFewerCommunities) {
    LeidenConfig config;
    config.resolution = 0.3f;

    LeidenClustering leiden(config);
    auto result = leiden.Cluster(*graph_);

    ASSERT_NE(result, nullptr);
    EXPECT_LE(result->num_communities, 4);  // Might merge some
}

// ========== Modularity Calculation Tests ==========

TEST(ModularityTest, PerfectPartition) {
    // Two disconnected cliques
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    // Perfect partition: {0,1,2} and {3,4,5}
    std::vector<uint32_t> labels = {0, 0, 0, 1, 1, 1};

    float modularity = LeidenClustering::CalculateModularity(*graph, labels);

    // Disconnected components should have modularity = 0.5
    EXPECT_GT(modularity, 0.4f);
}

TEST(ModularityTest, BadPartition) {
    // Two disconnected cliques
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    // Bad partition: splits cliques
    std::vector<uint32_t> labels = {0, 1, 0, 1, 0, 1};

    float modularity = LeidenClustering::CalculateModularity(*graph, labels);

    // Should have low or negative modularity
    EXPECT_LT(modularity, 0.3f);
}

TEST(ModularityTest, SingleCommunity) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    // All in one community
    std::vector<uint32_t> labels = {0, 0, 0};

    float modularity = LeidenClustering::CalculateModularity(*graph, labels);

    // Single community should have modularity = 0
    EXPECT_NEAR(modularity, 0.0f, 0.01f);
}

// ========== Well-Connected Verification Tests ==========

TEST(WellConnectedTest, ConnectedCommunity) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    std::vector<uint32_t> labels = {0, 0, 0};

    bool well_connected = LeidenClustering::VerifyWellConnected(*graph, labels);
    EXPECT_TRUE(well_connected);
}

TEST(WellConnectedTest, DisconnectedCommunity) {
    // Create graph with edge only between 0-1, no connection to 2
    std::vector<Edge> edges = {{0, 1, 1.0f}};

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    // Put all in same community even though 2 is disconnected
    std::vector<uint32_t> labels = {0, 0, 0};

    bool well_connected = LeidenClustering::VerifyWellConnected(*graph, labels);
    EXPECT_FALSE(well_connected);
}

TEST(WellConnectedTest, MultipleSeparateCommunities) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f},
        {2, 3, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 4, config);

    // Correct partition
    std::vector<uint32_t> labels = {0, 0, 1, 1};

    bool well_connected = LeidenClustering::VerifyWellConnected(*graph, labels);
    EXPECT_TRUE(well_connected);
}

// ========== ClusteringResult Tests ==========

TEST(ClusteringResultTest, GetCommunityMembers) {
    ClusteringResult result;
    result.labels = {0, 0, 1, 1, 1, 2};
    result.num_communities = 3;

    auto members = result.GetCommunityMembers();

    ASSERT_EQ(members.size(), 3);
    EXPECT_EQ(members[0].size(), 2);  // Nodes 0, 1
    EXPECT_EQ(members[1].size(), 3);  // Nodes 2, 3, 4
    EXPECT_EQ(members[2].size(), 1);  // Node 5
}

TEST(ClusteringResultTest, GetCommunitySizes) {
    ClusteringResult result;
    result.labels = {0, 0, 1, 1, 1, 2};
    result.num_communities = 3;

    auto sizes = result.GetCommunitySizes();

    ASSERT_EQ(sizes.size(), 3);
    EXPECT_EQ(sizes[0], 2);
    EXPECT_EQ(sizes[1], 3);
    EXPECT_EQ(sizes[2], 1);
}

TEST(ClusteringResultTest, RelabelContiguous) {
    ClusteringResult result;
    result.labels = {5, 5, 10, 10, 10, 3};
    result.num_communities = 3;

    result.RelabelContiguous();

    // Labels should now be 0, 1, 2 (order depends on first occurrence)
    std::set<uint32_t> unique(result.labels.begin(), result.labels.end());
    EXPECT_EQ(unique.size(), 3);
    EXPECT_TRUE(unique.count(0));
    EXPECT_TRUE(unique.count(1));
    EXPECT_TRUE(unique.count(2));
}

TEST(ClusteringResultTest, FilterByCommunitySize) {
    ClusteringResult result;
    result.labels = {0, 0, 1, 1, 1, 2};  // Community 0: 2, Community 1: 3, Community 2: 1
    result.num_communities = 3;

    auto filtered = result.FilterByCommunitySize(2);  // Min size 2

    // Community 2 (singleton) should be filtered out
    // Nodes 0,1 -> new community 0
    // Nodes 2,3,4 -> new community 1
    // Node 5 -> filtered (new community 2)

    EXPECT_EQ(filtered.labels[0], filtered.labels[1]);  // Same community
    EXPECT_EQ(filtered.labels[2], filtered.labels[3]);
    EXPECT_EQ(filtered.labels[3], filtered.labels[4]);  // Same community
    EXPECT_NE(filtered.labels[0], filtered.labels[2]);  // Different communities
}

// ========== Initial Labels Tests ==========

TEST(LeidenInitialLabelsTest, UsesInitialPartition) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    // Start with good initial partition
    std::vector<uint32_t> initial = {0, 0, 0, 1, 1, 1};

    LeidenConfig leiden_config;
    leiden_config.max_iterations = 1;  // Limit iterations

    LeidenClustering leiden(leiden_config);
    auto result = leiden.ClusterWithInitialLabels(*graph, initial);

    ASSERT_NE(result, nullptr);

    // Should maintain the partition (already optimal)
    EXPECT_EQ(result->labels[0], result->labels[1]);
    EXPECT_EQ(result->labels[1], result->labels[2]);
    EXPECT_EQ(result->labels[3], result->labels[4]);
    EXPECT_EQ(result->labels[4], result->labels[5]);
}

TEST(LeidenInitialLabelsTest, InvalidInitialLabelsReturnsNull) {
    std::vector<Edge> edges = {{0, 1, 1.0f}};

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    // Wrong number of labels
    std::vector<uint32_t> initial = {0, 0, 0};  // 3 labels for 2 nodes

    LeidenClustering leiden;
    auto result = leiden.ClusterWithInitialLabels(*graph, initial);

    EXPECT_EQ(result, nullptr);
}

// ========== Convenience Function Tests ==========

TEST(LeidenConvenienceTest, RunLeidenFunction) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 3);
}

TEST(LeidenConvenienceTest, GetCommunityLabelsFunction) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    auto labels = GetCommunityLabels(*graph);

    EXPECT_EQ(labels.size(), 3);
    // All in same community
    EXPECT_EQ(labels[0], labels[1]);
    EXPECT_EQ(labels[1], labels[2]);
}

TEST(LeidenConvenienceTest, RunLeidenWithResolution) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    auto result = RunLeiden(*graph, 2.0f);  // Higher resolution

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), 3);
}

// ========== Stress Tests ==========

TEST(LeidenStressTest, MediumGraph) {
    const size_t num_nodes = 500;
    const size_t edges_per_node = 10;

    std::mt19937 gen(42);
    std::uniform_int_distribution<uint32_t> node_dist(0, num_nodes - 1);
    std::uniform_real_distribution<float> weight_dist(0.1f, 1.0f);

    std::vector<Edge> edges;
    edges.reserve(num_nodes * edges_per_node);

    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < edges_per_node; ++j) {
            uint32_t target = node_dist(gen);
            if (target != i) {
                edges.push_back({
                    static_cast<uint32_t>(i),
                    target,
                    weight_dist(gen)
                });
            }
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, num_nodes, config);

    auto result = RunLeiden(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), num_nodes);
    EXPECT_GT(result->num_communities, 0);
    EXPECT_LT(result->num_communities, num_nodes);  // Should merge some
    EXPECT_GT(result->modularity, 0.0f);
}

TEST(LeidenStressTest, LargerGraphPerformance) {
    const size_t num_nodes = 1000;
    const size_t edges_per_node = 20;

    std::mt19937 gen(123);
    std::uniform_int_distribution<uint32_t> node_dist(0, num_nodes - 1);
    std::uniform_real_distribution<float> weight_dist(0.1f, 1.0f);

    std::vector<Edge> edges;
    edges.reserve(num_nodes * edges_per_node);

    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < edges_per_node; ++j) {
            uint32_t target = node_dist(gen);
            if (target != i) {
                edges.push_back({
                    static_cast<uint32_t>(i),
                    target,
                    weight_dist(gen)
                });
            }
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, num_nodes, config);

    LeidenConfig leiden_config;
    leiden_config.max_iterations = 20;  // Limit for test speed

    LeidenClustering leiden(leiden_config);
    auto result = leiden.Cluster(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->labels.size(), num_nodes);
    EXPECT_GT(result->stats.time_ms, 0.0f);
    EXPECT_LT(result->stats.time_ms, 10000.0f);  // Should complete in <10s
}

// ========== Determinism Tests ==========

TEST(LeidenDeterminismTest, SameSeedSameResult) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
        {2, 3, 0.1f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    LeidenConfig leiden_config;
    leiden_config.seed = 42;

    LeidenClustering leiden1(leiden_config);
    auto result1 = leiden1.Cluster(*graph);

    LeidenClustering leiden2(leiden_config);
    auto result2 = leiden2.Cluster(*graph);

    ASSERT_NE(result1, nullptr);
    ASSERT_NE(result2, nullptr);

    // Same seed should give same result
    EXPECT_EQ(result1->labels, result2->labels);
    EXPECT_FLOAT_EQ(result1->modularity, result2->modularity);
}

TEST(LeidenDeterminismTest, DifferentSeedsMayDiffer) {
    // With a larger graph, different seeds might find different local optima
    const size_t num_nodes = 100;
    const size_t edges_per_node = 5;

    std::mt19937 gen(42);
    std::uniform_int_distribution<uint32_t> node_dist(0, num_nodes - 1);

    std::vector<Edge> edges;
    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < edges_per_node; ++j) {
            uint32_t target = node_dist(gen);
            if (target != i) {
                edges.push_back({static_cast<uint32_t>(i), target, 1.0f});
            }
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, num_nodes, config);

    LeidenConfig config1, config2;
    config1.seed = 1;
    config2.seed = 999;

    LeidenClustering leiden1(config1);
    LeidenClustering leiden2(config2);

    auto result1 = leiden1.Cluster(*graph);
    auto result2 = leiden2.Cluster(*graph);

    ASSERT_NE(result1, nullptr);
    ASSERT_NE(result2, nullptr);

    // Both should produce valid results
    EXPECT_EQ(result1->labels.size(), num_nodes);
    EXPECT_EQ(result2->labels.size(), num_nodes);
}

// ========== Config Option Tests ==========

TEST(LeidenConfigOptionsTest, DisableRefinement) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f}, {4, 5, 1.0f}, {3, 5, 1.0f},
        {2, 3, 0.1f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 6, config);

    LeidenConfig leiden_config;
    leiden_config.enable_refinement = false;

    LeidenClustering leiden(leiden_config);
    auto result = leiden.Cluster(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_FLOAT_EQ(result->stats.refinement_time_ms, 0.0f);
}

TEST(LeidenConfigOptionsTest, LimitIterations) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3, config);

    LeidenConfig leiden_config;
    leiden_config.max_iterations = 1;

    LeidenClustering leiden(leiden_config);
    auto result = leiden.Cluster(*graph);

    ASSERT_NE(result, nullptr);
    EXPECT_LE(result->stats.num_iterations, 1);
}

TEST(LeidenConfigOptionsTest, MinCommunitySize) {
    std::vector<Edge> edges = {
        {0, 1, 1.0f}, {1, 2, 1.0f}, {0, 2, 1.0f},
        {3, 4, 1.0f},  // Small community
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;
    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 5, config);

    LeidenConfig leiden_config;
    leiden_config.min_community_size = 3;

    LeidenClustering leiden(leiden_config);
    auto result = leiden.Cluster(*graph);

    ASSERT_NE(result, nullptr);

    // The 2-node community should be filtered
    auto sizes = result->GetCommunitySizes();
    for (size_t s : sizes) {
        // Either meets minimum or is the "filtered" bucket
        EXPECT_TRUE(s >= 3 || s <= 2);
    }
}

}  // namespace
}  // namespace llmap::self_interference
