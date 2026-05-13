// LLmap — Unit tests for member propagation module.

#include <gtest/gtest.h>

#include "core/wave_state.h"
#include "reference_collapse/member_propagation.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace llmap {
namespace {

// ========== Configuration Tests ==========

TEST(MemberPropagationConfigTest, DefaultValues) {
    MemberPropagationConfig config;
    EXPECT_FLOAT_EQ(config.base_confidence_scaling, 0.8f);
    EXPECT_FLOAT_EQ(config.similarity_weight, 0.5f);
    EXPECT_FLOAT_EQ(config.min_confidence, 0.1f);
    EXPECT_FLOAT_EQ(config.max_confidence, 0.95f);
    EXPECT_FLOAT_EQ(config.cohesion_weight, 0.3f);
    EXPECT_FLOAT_EQ(config.outlier_penalty, 0.5f);
    EXPECT_TRUE(config.propagate_top_candidates_only);
    EXPECT_EQ(config.max_propagated_candidates, 10u);
}

TEST(MemberPropagationConfigTest, ConfigCanBeModified) {
    MemberPropagationConfig config;
    config.base_confidence_scaling = 0.9f;
    config.similarity_weight = 0.7f;
    config.min_confidence = 0.2f;
    config.max_confidence = 0.99f;
    config.cohesion_weight = 0.5f;
    config.outlier_penalty = 0.3f;
    config.propagate_top_candidates_only = false;
    config.max_propagated_candidates = 20;

    EXPECT_FLOAT_EQ(config.base_confidence_scaling, 0.9f);
    EXPECT_FLOAT_EQ(config.similarity_weight, 0.7f);
    EXPECT_FLOAT_EQ(config.min_confidence, 0.2f);
    EXPECT_FLOAT_EQ(config.max_confidence, 0.99f);
    EXPECT_FLOAT_EQ(config.cohesion_weight, 0.5f);
    EXPECT_FLOAT_EQ(config.outlier_penalty, 0.3f);
    EXPECT_FALSE(config.propagate_top_candidates_only);
    EXPECT_EQ(config.max_propagated_candidates, 20u);
}

// ========== Stats Structure Tests ==========

TEST(MemberPropagationStatsTest, DefaultValues) {
    MemberPropagationStats stats;
    EXPECT_EQ(stats.num_clusters, 0u);
    EXPECT_EQ(stats.num_representatives, 0u);
    EXPECT_EQ(stats.num_members_propagated, 0u);
    EXPECT_EQ(stats.num_members_skipped, 0u);
    EXPECT_FLOAT_EQ(stats.avg_member_confidence, 0.0f);
    EXPECT_FLOAT_EQ(stats.avg_similarity_to_rep, 0.0f);
    EXPECT_FLOAT_EQ(stats.avg_cluster_cohesion, 0.0f);
    EXPECT_EQ(stats.total_buckets_propagated, 0u);
    EXPECT_FLOAT_EQ(stats.propagation_time_ms, 0.0f);
}

TEST(MemberPropagationEntryTest, DefaultValues) {
    MemberPropagationEntry entry{};
    entry.member_idx = 5;
    entry.cluster_id = 2;
    entry.representative_idx = 3;
    entry.similarity_to_rep = 0.8f;
    entry.propagated_confidence = 0.7f;
    entry.buckets_received = 5;

    EXPECT_EQ(entry.member_idx, 5u);
    EXPECT_EQ(entry.cluster_id, 2u);
    EXPECT_EQ(entry.representative_idx, 3u);
    EXPECT_FLOAT_EQ(entry.similarity_to_rep, 0.8f);
    EXPECT_FLOAT_EQ(entry.propagated_confidence, 0.7f);
    EXPECT_EQ(entry.buckets_received, 5u);
}

// ========== Helper Classes for Testing ==========

class MemberPropagationTestBase : public ::testing::Test {
protected:
    std::unique_ptr<self_interference::SimilarityGraph> CreateClusterGraph(
        std::size_t num_clusters,
        std::size_t nodes_per_cluster,
        float intra_weight = 1.0f,
        float inter_weight = 0.1f)
    {
        std::vector<self_interference::Edge> edges;
        const std::size_t total_nodes = num_clusters * nodes_per_cluster;

        for (std::size_t c = 0; c < num_clusters; ++c) {
            std::size_t base = c * nodes_per_cluster;
            for (std::size_t i = 0; i < nodes_per_cluster; ++i) {
                for (std::size_t j = i + 1; j < nodes_per_cluster; ++j) {
                    edges.push_back({
                        static_cast<std::uint32_t>(base + i),
                        static_cast<std::uint32_t>(base + j),
                        intra_weight
                    });
                }
            }
        }

        for (std::size_t c = 0; c + 1 < num_clusters; ++c) {
            std::uint32_t from = static_cast<std::uint32_t>((c + 1) * nodes_per_cluster - 1);
            std::uint32_t to = static_cast<std::uint32_t>((c + 1) * nodes_per_cluster);
            edges.push_back({from, to, inter_weight});
        }

        self_interference::SimilarityGraphConfig config;
        config.make_symmetric = true;
        return self_interference::SimilarityGraph::BuildFromEdgeList(
            edges, total_nodes, config);
    }

    std::unique_ptr<self_interference::ClusteringResult> CreateClusteringResult(
        std::size_t num_clusters,
        std::size_t nodes_per_cluster)
    {
        auto result = std::make_unique<self_interference::ClusteringResult>();
        const std::size_t total_nodes = num_clusters * nodes_per_cluster;

        result->labels.resize(total_nodes);
        for (std::size_t c = 0; c < num_clusters; ++c) {
            for (std::size_t i = 0; i < nodes_per_cluster; ++i) {
                result->labels[c * nodes_per_cluster + i] = static_cast<std::uint32_t>(c);
            }
        }
        result->num_communities = num_clusters;
        result->modularity = 0.5f;

        return result;
    }

    std::unique_ptr<self_interference::ClusterRepResult> CreateRepResult(
        std::size_t num_clusters,
        std::size_t nodes_per_cluster)
    {
        auto result = std::make_unique<self_interference::ClusterRepResult>();

        for (std::size_t c = 0; c < num_clusters; ++c) {
            self_interference::RepresentativeInfo info;
            info.cluster_id = static_cast<std::uint32_t>(c);
            info.read_idx = static_cast<std::uint32_t>(c * nodes_per_cluster);
            info.confidence = 0.9f;
            info.cluster_size = nodes_per_cluster;
            info.avg_distance_to_members = 0.2f;
            info.centrality_score = 0.8f;

            result->representatives.push_back(info);
            result->representative_reads.push_back(info.read_idx);
        }

        result->stats.num_clusters = num_clusters;
        result->stats.num_representatives = num_clusters;

        return result;
    }

    WaveState CreateWaveStateWithRepCandidates(
        std::size_t num_reads,
        const self_interference::ClusterRepResult& rep_result,
        std::size_t buckets_per_rep = 5)
    {
        WaveState state(static_cast<std::uint32_t>(num_reads), WaveLevel::L1);

        for (const auto& rep_info : rep_result.representatives) {
            std::vector<BucketProb> candidates;
            float prob_sum = 0.0f;

            for (std::size_t i = 0; i < buckets_per_rep; ++i) {
                BucketProb bp;
                bp.bucket_id = static_cast<std::uint32_t>(rep_info.cluster_id * 10 + i);
                bp.probability = 1.0f / static_cast<float>(i + 1);
                prob_sum += bp.probability;
                candidates.push_back(bp);
            }

            for (auto& bp : candidates) {
                bp.probability /= prob_sum;
            }

            std::sort(candidates.begin(), candidates.end(),
                [](const BucketProb& a, const BucketProb& b) {
                    return a.bucket_id < b.bucket_id;
                });

            state.set_read_candidates(rep_info.read_idx, candidates);
        }

        return state;
    }
};

// ========== Basic Propagation Tests ==========

class MemberPropagationBasicTest : public MemberPropagationTestBase {
protected:
    void SetUp() override {
        graph_ = CreateClusterGraph(2, 5);
        clustering_ = CreateClusteringResult(2, 5);
        rep_result_ = CreateRepResult(2, 5);
        state_ = CreateWaveStateWithRepCandidates(10, *rep_result_);
    }

    std::unique_ptr<self_interference::SimilarityGraph> graph_;
    std::unique_ptr<self_interference::ClusteringResult> clustering_;
    std::unique_ptr<self_interference::ClusterRepResult> rep_result_;
    WaveState state_;
};

TEST_F(MemberPropagationBasicTest, PropagatePopulatesMembers) {
    MemberPropagation propagation;
    auto result = propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    EXPECT_GT(result.stats.num_members_propagated, 0u);
    EXPECT_EQ(result.stats.num_clusters, 2u);
    EXPECT_EQ(result.stats.num_representatives, 2u);
}

TEST_F(MemberPropagationBasicTest, MembersGetBucketsFromRep) {
    MemberPropagation propagation;
    auto result = propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    for (std::uint32_t i = 1; i < 5; ++i) {
        auto buckets = state_.bucket_indices_for_read(i);
        EXPECT_GT(buckets.size(), 0u) << "Member " << i << " should have buckets";
    }
}

TEST_F(MemberPropagationBasicTest, RepresentativesNotModified) {
    auto rep0_buckets_before = state_.bucket_indices_for_read(0);
    auto rep5_buckets_before = state_.bucket_indices_for_read(5);

    std::vector<std::uint32_t> rep0_copy(rep0_buckets_before.begin(), rep0_buckets_before.end());
    std::vector<std::uint32_t> rep5_copy(rep5_buckets_before.begin(), rep5_buckets_before.end());

    MemberPropagation propagation;
    propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    auto rep0_buckets_after = state_.bucket_indices_for_read(0);
    auto rep5_buckets_after = state_.bucket_indices_for_read(5);

    EXPECT_EQ(rep0_buckets_after.size(), rep0_copy.size());
    EXPECT_EQ(rep5_buckets_after.size(), rep5_copy.size());
}

TEST_F(MemberPropagationBasicTest, ConvenienceFunctionWorks) {
    auto result = PropagateToMembers(state_, *clustering_, *rep_result_, *graph_);

    EXPECT_GT(result.stats.num_members_propagated, 0u);
    EXPECT_GT(result.stats.propagation_time_ms, 0.0f);
}

// ========== Single Cluster Propagation Tests ==========

TEST_F(MemberPropagationBasicTest, PropagateCluster_SingleCluster) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members = {0, 1, 2, 3, 4};
    auto entries = propagation.PropagateCluster(
        state_, 0, members, *graph_, 0);

    EXPECT_EQ(entries.size(), 4u);
    for (const auto& entry : entries) {
        EXPECT_NE(entry.member_idx, 0u);
        EXPECT_EQ(entry.cluster_id, 0u);
        EXPECT_EQ(entry.representative_idx, 0u);
        EXPECT_GT(entry.similarity_to_rep, 0.0f);
        EXPECT_GT(entry.propagated_confidence, 0.0f);
        EXPECT_GT(entry.buckets_received, 0u);
    }
}

TEST_F(MemberPropagationBasicTest, PropagateCluster_SecondCluster) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members = {5, 6, 7, 8, 9};
    auto entries = propagation.PropagateCluster(
        state_, 5, members, *graph_, 1);

    EXPECT_EQ(entries.size(), 4u);
    for (const auto& entry : entries) {
        EXPECT_GE(entry.member_idx, 6u);
        EXPECT_LE(entry.member_idx, 9u);
        EXPECT_EQ(entry.cluster_id, 1u);
        EXPECT_EQ(entry.representative_idx, 5u);
    }
}

// ========== Single Member Propagation Tests ==========

TEST_F(MemberPropagationBasicTest, PropagateMember_ReturnsBucketCount) {
    MemberPropagation propagation;

    std::uint32_t buckets = propagation.PropagateMember(
        state_, 0, 1, 0.8f, 0.7f);

    EXPECT_GT(buckets, 0u);
}

TEST_F(MemberPropagationBasicTest, PropagateMember_MemberGetsBuckets) {
    MemberPropagation propagation;

    propagation.PropagateMember(state_, 0, 1, 0.9f, 0.8f);

    auto member_buckets = state_.bucket_indices_for_read(1);
    auto member_probs = state_.probabilities_for_read(1);

    EXPECT_GT(member_buckets.size(), 0u);
    EXPECT_EQ(member_buckets.size(), member_probs.size());
}

TEST_F(MemberPropagationBasicTest, PropagateMember_ProbabilitiesNormalized) {
    MemberPropagation propagation;

    propagation.PropagateMember(state_, 0, 1, 0.9f, 0.8f);

    auto member_probs = state_.probabilities_for_read(1);
    float sum = 0.0f;
    for (float p : member_probs) {
        sum += p;
    }

    EXPECT_NEAR(sum, 1.0f, 0.01f);
}

// ========== Cluster Cohesion Tests ==========

TEST_F(MemberPropagationBasicTest, ComputeClusterCohesion_HighCohesion) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members = {0, 1, 2, 3, 4};
    float cohesion = propagation.ComputeClusterCohesion(*graph_, members);

    EXPECT_GT(cohesion, 0.5f);
}

TEST_F(MemberPropagationBasicTest, ComputeClusterCohesion_SingleMember) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members = {0};
    float cohesion = propagation.ComputeClusterCohesion(*graph_, members);

    EXPECT_FLOAT_EQ(cohesion, 1.0f);
}

TEST_F(MemberPropagationBasicTest, ComputeClusterCohesion_Empty) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members;
    float cohesion = propagation.ComputeClusterCohesion(*graph_, members);

    EXPECT_FLOAT_EQ(cohesion, 1.0f);
}

// ========== Confidence Scaling Tests ==========

TEST_F(MemberPropagationBasicTest, PropagationReportsConfidence) {
    MemberPropagation propagation;

    std::vector<std::uint32_t> members = {0, 1, 2, 3, 4};
    auto entries = propagation.PropagateCluster(
        state_, 0, members, *graph_, 0);

    for (const auto& entry : entries) {
        EXPECT_GT(entry.propagated_confidence, 0.0f);
        EXPECT_LE(entry.propagated_confidence, 1.0f);
    }
}

TEST_F(MemberPropagationBasicTest, OutlierPenaltyReducesConfidence) {
    MemberPropagationConfig config;
    config.outlier_penalty = 0.1f;
    MemberPropagation propagation(config);

    std::vector<std::uint32_t> members = {0, 1, 2};
    auto entries = propagation.PropagateCluster(
        state_, 0, members, *graph_, 0);

    for (const auto& entry : entries) {
        if (entry.similarity_to_rep < 0.3f) {
            EXPECT_LT(entry.propagated_confidence, config.max_confidence);
        }
    }
}

// ========== Result Query Tests ==========

TEST_F(MemberPropagationBasicTest, ResultGetEntry) {
    MemberPropagation propagation;
    auto result = propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    auto entry = result.GetEntry(1);
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->member_idx, 1u);
}

TEST_F(MemberPropagationBasicTest, ResultGetEntry_NotFound) {
    MemberPropagation propagation;
    auto result = propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    auto entry = result.GetEntry(0);
    EXPECT_FALSE(entry.has_value());
}

TEST_F(MemberPropagationBasicTest, ResultGetClusterMembers) {
    MemberPropagation propagation;
    auto result = propagation.Propagate(state_, *clustering_, *rep_result_, *graph_);

    auto cluster0_members = result.GetClusterMembers(0);
    auto cluster1_members = result.GetClusterMembers(1);

    EXPECT_EQ(cluster0_members.size(), 4u);
    EXPECT_EQ(cluster1_members.size(), 4u);
}

// ========== Empty/Edge Case Tests ==========

TEST(MemberPropagationEdgeCaseTest, EmptyClusteringResult) {
    MemberPropagation propagation;
    WaveState state(0);
    self_interference::ClusteringResult clustering;
    self_interference::ClusterRepResult rep_result;

    auto graph = self_interference::SimilarityGraph::BuildFromEdgeList({}, 0, {});

    auto result = propagation.Propagate(state, clustering, rep_result, *graph);

    EXPECT_EQ(result.stats.num_members_propagated, 0u);
}

TEST(MemberPropagationEdgeCaseTest, NoRepresentatives) {
    WaveState state(10);

    self_interference::ClusteringResult clustering;
    clustering.labels.resize(10, 0);
    clustering.num_communities = 1;

    self_interference::ClusterRepResult rep_result;

    std::vector<self_interference::Edge> edges;
    for (std::uint32_t i = 0; i < 9; ++i) {
        edges.push_back({i, static_cast<std::uint32_t>(i + 1), 0.8f});
    }
    auto graph = self_interference::SimilarityGraph::BuildFromEdgeList(edges, 10, {});

    MemberPropagation propagation;
    auto result = propagation.Propagate(state, clustering, rep_result, *graph);

    EXPECT_EQ(result.stats.num_members_propagated, 0u);
}

TEST(MemberPropagationEdgeCaseTest, SingleMemberCluster) {
    WaveState state(1, WaveLevel::L1);

    std::vector<BucketProb> candidates = {{0, 0.5f}, {1, 0.5f}};
    state.set_read_candidates(0, candidates);

    self_interference::ClusteringResult clustering;
    clustering.labels = {0};
    clustering.num_communities = 1;

    self_interference::ClusterRepResult rep_result;
    self_interference::RepresentativeInfo rep_info;
    rep_info.cluster_id = 0;
    rep_info.read_idx = 0;
    rep_info.confidence = 0.9f;
    rep_info.cluster_size = 1;
    rep_result.representatives.push_back(rep_info);
    rep_result.representative_reads.push_back(0);
    rep_result.stats.num_representatives = 1;

    auto graph = self_interference::SimilarityGraph::BuildFromEdgeList({}, 1, {});

    MemberPropagation propagation;
    auto result = propagation.Propagate(state, clustering, rep_result, *graph);

    EXPECT_EQ(result.stats.num_members_propagated, 0u);
}

// ========== Move Semantics Tests ==========

TEST(MemberPropagationMoveTest, MoveConstruction) {
    MemberPropagationConfig config;
    config.base_confidence_scaling = 0.75f;

    MemberPropagation propagation1(config);
    MemberPropagation propagation2(std::move(propagation1));

    EXPECT_FLOAT_EQ(propagation2.Config().base_confidence_scaling, 0.75f);
}

TEST(MemberPropagationMoveTest, MoveAssignment) {
    MemberPropagationConfig config;
    config.similarity_weight = 0.6f;

    MemberPropagation propagation1(config);
    MemberPropagation propagation2;

    propagation2 = std::move(propagation1);

    EXPECT_FLOAT_EQ(propagation2.Config().similarity_weight, 0.6f);
}

}  // namespace
}  // namespace llmap
