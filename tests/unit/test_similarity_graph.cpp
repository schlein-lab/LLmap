#include <gtest/gtest.h>

#include "self_interference/similarity_graph.h"
#include "self_interference/faiss_wrapper.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <vector>

namespace llmap::self_interference {
namespace {

// ========== Configuration Tests ==========

TEST(SimilarityGraphConfigTest, DefaultValues) {
    SimilarityGraphConfig config;
    EXPECT_TRUE(config.convert_distance_to_similarity);
    EXPECT_FLOAT_EQ(config.distance_scale, 1.0f);
    EXPECT_FLOAT_EQ(config.min_weight_threshold, 0.0f);
    EXPECT_FLOAT_EQ(config.max_distance_threshold, -1.0f);
    EXPECT_EQ(config.max_edges_per_node, 0);
    EXPECT_TRUE(config.make_symmetric);
    EXPECT_TRUE(config.remove_self_loops);
    EXPECT_TRUE(config.merge_duplicates);
    EXPECT_EQ(config.expected_num_nodes, 0);
}

TEST(SimilarityGraphConfigTest, ConfigCanBeModified) {
    SimilarityGraphConfig config;
    config.convert_distance_to_similarity = false;
    config.distance_scale = 0.5f;
    config.min_weight_threshold = 0.1f;
    config.max_distance_threshold = 2.0f;
    config.max_edges_per_node = 10;
    config.make_symmetric = false;
    config.remove_self_loops = false;
    config.merge_duplicates = false;
    config.expected_num_nodes = 1000;

    EXPECT_FALSE(config.convert_distance_to_similarity);
    EXPECT_FLOAT_EQ(config.distance_scale, 0.5f);
    EXPECT_FLOAT_EQ(config.min_weight_threshold, 0.1f);
    EXPECT_FLOAT_EQ(config.max_distance_threshold, 2.0f);
    EXPECT_EQ(config.max_edges_per_node, 10);
    EXPECT_FALSE(config.make_symmetric);
    EXPECT_FALSE(config.remove_self_loops);
    EXPECT_FALSE(config.merge_duplicates);
    EXPECT_EQ(config.expected_num_nodes, 1000);
}

// ========== Utility Function Tests ==========

TEST(SimilarityGraphUtilsTest, DistanceToSimilarity_ZeroDistance) {
    EXPECT_FLOAT_EQ(DistanceToSimilarity(0.0f), 1.0f);
}

TEST(SimilarityGraphUtilsTest, DistanceToSimilarity_UnitDistance) {
    EXPECT_FLOAT_EQ(DistanceToSimilarity(1.0f), 0.5f);
}

TEST(SimilarityGraphUtilsTest, DistanceToSimilarity_LargeDistance) {
    float sim = DistanceToSimilarity(9.0f);
    EXPECT_FLOAT_EQ(sim, 0.1f);
}

TEST(SimilarityGraphUtilsTest, DistanceToSimilarity_WithScale) {
    float sim = DistanceToSimilarity(2.0f, 0.5f);
    EXPECT_FLOAT_EQ(sim, 0.5f);  // 1 / (1 + 2 * 0.5) = 0.5
}

TEST(SimilarityGraphUtilsTest, InnerProductToDistance) {
    EXPECT_FLOAT_EQ(InnerProductToDistance(1.0f), 0.0f);
    EXPECT_FLOAT_EQ(InnerProductToDistance(0.0f), 1.0f);
    EXPECT_FLOAT_EQ(InnerProductToDistance(0.5f), 0.5f);
}

// ========== Edge Tests ==========

TEST(EdgeTest, Equality) {
    Edge e1{0, 1, 0.5f};
    Edge e2{0, 1, 0.8f};  // Different weight
    Edge e3{0, 2, 0.5f};  // Different target

    EXPECT_TRUE(e1 == e2);   // Only source/target matter
    EXPECT_FALSE(e1 == e3);
}

// ========== BuildFromKnn Tests ==========

class SimilarityGraphKnnTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create simple k-NN result: 5 queries, k=3
        // Query 0: neighbors 1, 2, 3
        // Query 1: neighbors 0, 2, 4
        // etc.
        num_queries_ = 5;
        k_ = 3;

        indices_ = {
            1, 2, 3,   // Query 0
            0, 2, 4,   // Query 1
            0, 1, 3,   // Query 2
            0, 2, 4,   // Query 3
            1, 3, 0,   // Query 4
        };

        distances_ = {
            0.1f, 0.2f, 0.3f,   // Query 0
            0.1f, 0.15f, 0.25f, // Query 1
            0.2f, 0.3f, 0.4f,   // Query 2
            0.1f, 0.2f, 0.5f,   // Query 3
            0.05f, 0.1f, 0.2f,  // Query 4
        };
    }

    size_t num_queries_;
    size_t k_;
    std::vector<int64_t> indices_;
    std::vector<float> distances_;
};

TEST_F(SimilarityGraphKnnTest, BasicConstruction) {
    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumNodes(), num_queries_);
    EXPECT_GT(graph->NumEdges(), 0);
    EXPECT_TRUE(graph->IsSymmetric());
    EXPECT_FALSE(graph->IsEmpty());
}

TEST_F(SimilarityGraphKnnTest, AsymmetricGraph) {
    SimilarityGraphConfig config;
    config.make_symmetric = false;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);
    EXPECT_FALSE(graph->IsSymmetric());

    // Directed edges: 5 queries * 3 neighbors = 15 edges (minus self-loops if any)
    EXPECT_LE(graph->NumEdges(), num_queries_ * k_);
}

TEST_F(SimilarityGraphKnnTest, SelfLoopsRemoved) {
    // Add a self-loop to the test data
    // indices_[0] is Query 0's first neighbor → set to 0 (itself)
    indices_[0] = 0;

    SimilarityGraphConfig config;
    config.remove_self_loops = true;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Verify no self-loops
    for (size_t i = 0; i < graph->NumNodes(); ++i) {
        EXPECT_FALSE(graph->HasEdge(i, i)) << "Self-loop found at node " << i;
    }

    auto stats = graph->GetStats();
    EXPECT_EQ(stats.num_self_loops_removed, 1);
}

TEST_F(SimilarityGraphKnnTest, DistanceToSimilarityConversion) {
    SimilarityGraphConfig config;
    config.convert_distance_to_similarity = true;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Check that weights are in (0, 1] range (similarity)
    auto stats = graph->GetStats();
    EXPECT_GT(stats.min_weight, 0.0f);
    EXPECT_LE(stats.max_weight, 1.0f);
}

TEST_F(SimilarityGraphKnnTest, RawDistanceWeights) {
    SimilarityGraphConfig config;
    config.convert_distance_to_similarity = false;
    config.make_symmetric = false;  // Directed graph to preserve exact distances
    config.merge_duplicates = false;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Weights should match original distances exactly (no merging)
    auto stats = graph->GetStats();
    float min_dist = *std::min_element(distances_.begin(), distances_.end());
    float max_dist = *std::max_element(distances_.begin(), distances_.end());

    EXPECT_FLOAT_EQ(stats.min_weight, min_dist);
    EXPECT_FLOAT_EQ(stats.max_weight, max_dist);
}

TEST_F(SimilarityGraphKnnTest, DistanceThreshold) {
    SimilarityGraphConfig config;
    config.max_distance_threshold = 0.15f;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Only edges with distance <= 0.15 should be present
    auto edges = graph->ToEdgeList();
    // With symmetric, the edge count is variable, but should be less than full
    EXPECT_LT(graph->NumEdges(), num_queries_ * k_ * 2);
}

TEST_F(SimilarityGraphKnnTest, WeightThreshold) {
    SimilarityGraphConfig config;
    config.convert_distance_to_similarity = true;
    config.min_weight_threshold = 0.8f;  // Only very similar pairs

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Few edges should remain with high threshold
    auto stats = graph->GetStats();
    EXPECT_GE(stats.min_weight, 0.8f);
}

TEST_F(SimilarityGraphKnnTest, MaxEdgesPerNode) {
    SimilarityGraphConfig config;
    config.max_edges_per_node = 1;
    config.make_symmetric = false;

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_, config);

    ASSERT_NE(graph, nullptr);

    // Each node should have at most 1 outgoing edge
    for (size_t i = 0; i < graph->NumNodes(); ++i) {
        EXPECT_LE(graph->Degree(i), 1) << "Node " << i << " has degree > 1";
    }
}

TEST_F(SimilarityGraphKnnTest, InvalidInput) {
    std::vector<int64_t> empty_indices;
    std::vector<float> empty_distances;

    auto graph = SimilarityGraph::BuildFromKnn(
        empty_indices, empty_distances, 0, 0);

    // Should handle gracefully
    ASSERT_NE(graph, nullptr);
    EXPECT_TRUE(graph->IsEmpty());
}

TEST_F(SimilarityGraphKnnTest, NegativeIndicesSkipped) {
    indices_[0] = -1;  // Invalid index
    indices_[1] = -1;  // Invalid index

    auto graph = SimilarityGraph::BuildFromKnn(indices_, distances_, num_queries_, k_);

    ASSERT_NE(graph, nullptr);
    // Should have fewer edges due to skipped invalid indices
}

// ========== BuildFromEdgeList Tests ==========

TEST(SimilarityGraphEdgeListTest, BasicConstruction) {
    std::vector<Edge> edges = {
        {0, 1, 0.5f},
        {1, 2, 0.3f},
        {2, 0, 0.7f},
    };

    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 3);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumNodes(), 3);
    EXPECT_GT(graph->NumEdges(), 0);
}

TEST(SimilarityGraphEdgeListTest, SymmetricFromDirected) {
    std::vector<Edge> edges = {
        {0, 1, 0.5f},
    };

    SimilarityGraphConfig config;
    config.make_symmetric = true;

    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    ASSERT_NE(graph, nullptr);

    // Both directions should exist
    EXPECT_TRUE(graph->HasEdge(0, 1));
    EXPECT_TRUE(graph->HasEdge(1, 0));
    EXPECT_EQ(graph->NumEdges(), 2);
}

TEST(SimilarityGraphEdgeListTest, DuplicatesMerged) {
    std::vector<Edge> edges = {
        {0, 1, 0.5f},
        {0, 1, 0.8f},  // Duplicate with different weight
        {0, 1, 0.3f},  // Another duplicate
    };

    SimilarityGraphConfig config;
    config.merge_duplicates = true;
    config.make_symmetric = false;

    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 2, config);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumEdges(), 1);  // Merged into one
    EXPECT_FLOAT_EQ(graph->GetEdgeWeight(0, 1), 0.8f);  // Max weight kept
}

TEST(SimilarityGraphEdgeListTest, EmptyEdgeList) {
    std::vector<Edge> edges;

    auto graph = SimilarityGraph::BuildFromEdgeList(edges, 5);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumNodes(), 5);
    EXPECT_EQ(graph->NumEdges(), 0);
    EXPECT_TRUE(graph->IsEmpty());
}

// ========== Graph Query Tests ==========

class SimilarityGraphQueryTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple graph: 0 -- 1 -- 2
        //                             \-- 3
        std::vector<Edge> edges = {
            {0, 1, 0.9f},
            {1, 2, 0.8f},
            {1, 3, 0.7f},
        };

        SimilarityGraphConfig config;
        config.make_symmetric = true;

        graph_ = SimilarityGraph::BuildFromEdgeList(edges, 4, config);
    }

    std::unique_ptr<SimilarityGraph> graph_;
};

TEST_F(SimilarityGraphQueryTest, Neighbors) {
    auto neighbors = graph_->Neighbors(1);

    EXPECT_EQ(neighbors.size(), 3);  // 0, 2, 3

    // Should be sorted
    std::vector<uint32_t> expected = {0, 2, 3};
    std::vector<uint32_t> actual(neighbors.begin(), neighbors.end());
    EXPECT_EQ(actual, expected);
}

TEST_F(SimilarityGraphQueryTest, NeighborWeights) {
    auto weights = graph_->NeighborWeights(1);

    EXPECT_EQ(weights.size(), 3);
}

TEST_F(SimilarityGraphQueryTest, Degree) {
    EXPECT_EQ(graph_->Degree(0), 1);  // Connected to 1
    EXPECT_EQ(graph_->Degree(1), 3);  // Connected to 0, 2, 3
    EXPECT_EQ(graph_->Degree(2), 1);  // Connected to 1
    EXPECT_EQ(graph_->Degree(3), 1);  // Connected to 1
}

TEST_F(SimilarityGraphQueryTest, HasEdge) {
    EXPECT_TRUE(graph_->HasEdge(0, 1));
    EXPECT_TRUE(graph_->HasEdge(1, 0));  // Symmetric
    EXPECT_TRUE(graph_->HasEdge(1, 2));
    EXPECT_TRUE(graph_->HasEdge(1, 3));

    EXPECT_FALSE(graph_->HasEdge(0, 2));  // Not connected
    EXPECT_FALSE(graph_->HasEdge(0, 3));
    EXPECT_FALSE(graph_->HasEdge(2, 3));
}

TEST_F(SimilarityGraphQueryTest, GetEdgeWeight) {
    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(0, 1), 0.9f);
    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(1, 0), 0.9f);  // Symmetric
    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(1, 2), 0.8f);
    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(1, 3), 0.7f);

    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(0, 2), 0.0f);  // Not connected
}

TEST_F(SimilarityGraphQueryTest, InvalidNodeQueries) {
    EXPECT_EQ(graph_->Degree(999), 0);
    EXPECT_TRUE(graph_->Neighbors(999).empty());
    EXPECT_FALSE(graph_->HasEdge(999, 0));
    EXPECT_FLOAT_EQ(graph_->GetEdgeWeight(999, 0), 0.0f);
}

TEST_F(SimilarityGraphQueryTest, Statistics) {
    auto stats = graph_->GetStats();

    EXPECT_EQ(stats.num_nodes, 4);
    EXPECT_EQ(stats.num_edges, 6);  // 3 edges * 2 (symmetric)
    EXPECT_FLOAT_EQ(stats.avg_degree, 1.5f);
    EXPECT_FLOAT_EQ(stats.max_weight, 0.9f);
    EXPECT_FLOAT_EQ(stats.min_weight, 0.7f);
}

TEST_F(SimilarityGraphQueryTest, AverageDegree) {
    EXPECT_FLOAT_EQ(graph_->AverageDegree(), 1.5f);
}

TEST_F(SimilarityGraphQueryTest, MaxDegree) {
    EXPECT_EQ(graph_->MaxDegree(), 3);  // Node 1 has degree 3
}

// ========== CSR Access Tests ==========

TEST_F(SimilarityGraphQueryTest, RawCsrAccess) {
    auto row_offsets = graph_->RowOffsets();
    auto col_indices = graph_->ColIndices();
    auto weights = graph_->Weights();

    EXPECT_EQ(row_offsets.size(), 5);  // num_nodes + 1
    EXPECT_EQ(col_indices.size(), 6);  // num_edges
    EXPECT_EQ(weights.size(), 6);

    // Row offsets should be non-decreasing
    for (size_t i = 1; i < row_offsets.size(); ++i) {
        EXPECT_GE(row_offsets[i], row_offsets[i - 1]);
    }
}

// ========== Export Tests ==========

TEST_F(SimilarityGraphQueryTest, ToEdgeList) {
    auto edges = graph_->ToEdgeList();

    EXPECT_EQ(edges.size(), 6);  // 3 * 2 (symmetric)

    // Check that all edges are present
    bool found_0_1 = false, found_1_2 = false, found_1_3 = false;
    for (const auto& e : edges) {
        if (e.source == 0 && e.target == 1) found_0_1 = true;
        if (e.source == 1 && e.target == 2) found_1_2 = true;
        if (e.source == 1 && e.target == 3) found_1_3 = true;
    }
    EXPECT_TRUE(found_0_1);
    EXPECT_TRUE(found_1_2);
    EXPECT_TRUE(found_1_3);
}

TEST_F(SimilarityGraphQueryTest, ToEdgeListString) {
    auto str = graph_->ToEdgeListString();

    EXPECT_FALSE(str.empty());
    EXPECT_NE(str.find("SimilarityGraph"), std::string::npos);
    EXPECT_NE(str.find("4 nodes"), std::string::npos);
    EXPECT_NE(str.find("6 edges"), std::string::npos);
}

// ========== Serialization Tests ==========

TEST_F(SimilarityGraphQueryTest, SaveAndLoad) {
    std::filesystem::path temp_path =
        std::filesystem::temp_directory_path() / "test_graph.llmg";

    ASSERT_TRUE(graph_->Save(temp_path.string()));

    auto loaded = SimilarityGraph::Load(temp_path.string());

    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(loaded->NumNodes(), graph_->NumNodes());
    EXPECT_EQ(loaded->NumEdges(), graph_->NumEdges());
    EXPECT_EQ(loaded->IsSymmetric(), graph_->IsSymmetric());

    // Verify edges
    EXPECT_TRUE(loaded->HasEdge(0, 1));
    EXPECT_FLOAT_EQ(loaded->GetEdgeWeight(0, 1), 0.9f);

    std::filesystem::remove(temp_path);
}

TEST(SimilarityGraphSerializationTest, LoadNonexistentFile) {
    auto graph = SimilarityGraph::Load("/nonexistent/path/file.llmg");
    EXPECT_EQ(graph, nullptr);
}

TEST(SimilarityGraphSerializationTest, LoadInvalidFile) {
    std::filesystem::path temp_path =
        std::filesystem::temp_directory_path() / "invalid.llmg";

    // Write invalid data
    std::ofstream f(temp_path, std::ios::binary);
    f << "INVALID DATA";
    f.close();

    auto graph = SimilarityGraph::Load(temp_path.string());
    EXPECT_EQ(graph, nullptr);

    std::filesystem::remove(temp_path);
}

// ========== Integration with BatchKnnResult ==========

TEST(SimilarityGraphIntegrationTest, BuildFromBatchKnnResult) {
    BatchKnnResult result;
    result.num_queries = 3;
    result.k = 2;
    result.indices = {1, 2, 0, 2, 0, 1};
    result.distances = {0.1f, 0.2f, 0.1f, 0.15f, 0.2f, 0.15f};

    auto graph = SimilarityGraph::BuildFromBatchKnnResult(result);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumNodes(), 3);
    EXPECT_GT(graph->NumEdges(), 0);
}

// ========== Stress Tests ==========

TEST(SimilarityGraphStressTest, LargeGraph) {
    const size_t num_nodes = 1000;
    const size_t k = 50;

    // Generate random k-NN results
    std::mt19937 gen(42);
    std::uniform_int_distribution<int64_t> node_dist(0, num_nodes - 1);
    std::uniform_real_distribution<float> dist_dist(0.0f, 2.0f);

    std::vector<int64_t> indices(num_nodes * k);
    std::vector<float> distances(num_nodes * k);

    for (size_t i = 0; i < num_nodes * k; ++i) {
        indices[i] = node_dist(gen);
        distances[i] = dist_dist(gen);
    }

    SimilarityGraphConfig config;
    config.expected_num_nodes = num_nodes;

    auto graph = SimilarityGraph::BuildFromKnn(indices, distances, num_nodes, k, config);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumNodes(), num_nodes);
    EXPECT_GT(graph->NumEdges(), 0);

    auto stats = graph->GetStats();
    EXPECT_GT(stats.avg_degree, 0.0f);
    EXPECT_GT(stats.build_time_ms, 0.0f);
}

TEST(SimilarityGraphStressTest, HighConnectivityGraph) {
    // Create a dense graph where every node connects to every other
    const size_t num_nodes = 100;

    std::vector<Edge> edges;
    edges.reserve(num_nodes * num_nodes);

    for (size_t i = 0; i < num_nodes; ++i) {
        for (size_t j = 0; j < num_nodes; ++j) {
            if (i != j) {
                edges.push_back({
                    static_cast<uint32_t>(i),
                    static_cast<uint32_t>(j),
                    1.0f / (1.0f + std::abs(static_cast<float>(i) - static_cast<float>(j)))
                });
            }
        }
    }

    SimilarityGraphConfig config;
    config.make_symmetric = false;  // Already have both directions

    auto graph = SimilarityGraph::BuildFromEdgeList(edges, num_nodes, config);

    ASSERT_NE(graph, nullptr);
    EXPECT_EQ(graph->NumEdges(), num_nodes * (num_nodes - 1));

    // Each node should have degree = num_nodes - 1
    for (size_t i = 0; i < num_nodes; ++i) {
        EXPECT_EQ(graph->Degree(i), num_nodes - 1);
    }
}

}  // namespace
}  // namespace llmap::self_interference
