// LLmap — AllpairPipeline processing stages.
//
// Contains: BuildSimilarityGraph, RunClustering, RunSelfWaveCollapse, SelectRepresentatives

#include "self_interference/allpair_pipeline.h"
#include "self_interference/allpair_pipeline_internal.h"

#include <algorithm>
#include <cstdio>

#include "self_interference/cluster_rep.h"
#include "self_interference/faiss_wrapper.h"
#include "self_interference/leiden_clustering.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

namespace llmap::self_interference {

bool AllpairPipeline::BuildSimilarityGraph() {
    const size_t num_reads = state_->reads.size();

    if (num_reads < 2) {
        state_->last_error = "Need at least 2 reads to build similarity graph";
        return false;
    }

    // Create FAISS index
    FaissIndexConfig faiss_config;
    faiss_config.embedding_dim = config_.embedding_dim;
    faiss_config.nlist = std::min(config_.faiss_nlist, num_reads / 10 + 1);
    faiss_config.nprobe = config_.faiss_nprobe;
    faiss_config.default_k = config_.faiss_k;
    faiss_config.normalize_vectors = true;

    if (config_.use_gpu_faiss) {
        faiss_config.provider = FaissProvider::GPU;
        faiss_config.gpu_device_id = config_.gpu_device_id;
    }

    // For small datasets, use exact search
    if (num_reads < 10000) {
        faiss_config.index_type = FaissIndexType::FlatL2;
    } else {
        faiss_config.index_type = FaissIndexType::IVFFlat;
    }

    state_->faiss_index = FaissIndex::Create(faiss_config);
    if (!state_->faiss_index) {
        state_->last_error = "Failed to create FAISS index";
        return false;
    }

    // Train and add vectors
    if (!state_->faiss_index->TrainAndAdd(state_->embeddings, num_reads)) {
        state_->last_error = "Failed to build FAISS index: " +
                             state_->faiss_index->LastError();
        return false;
    }

    // Search for k-NN
    const size_t k = std::min(config_.faiss_k, num_reads - 1);
    auto knn_result = state_->faiss_index->SearchBatch(
        state_->embeddings, num_reads, k);

    if (!knn_result) {
        state_->last_error = "k-NN search failed";
        return false;
    }

    // Build similarity graph from k-NN results
    SimilarityGraphConfig graph_config;
    graph_config.convert_distance_to_similarity = true;
    graph_config.make_symmetric = true;
    graph_config.remove_self_loops = true;

    state_->graph = SimilarityGraph::BuildFromKnn(
        knn_result->indices,
        knn_result->distances,
        num_reads,
        k,
        graph_config);

    if (!state_->graph) {
        state_->last_error = "Failed to build similarity graph";
        return false;
    }

    return true;
}

bool AllpairPipeline::RunClustering() {
    if (!state_->graph) {
        state_->last_error = "Similarity graph not built";
        return false;
    }

    LeidenConfig leiden_config;
    leiden_config.resolution = config_.leiden_resolution;
    leiden_config.max_iterations = config_.leiden_max_iterations;
    leiden_config.min_community_size = config_.min_cluster_size;

    LeidenClustering leiden(leiden_config);
    state_->clustering = leiden.Cluster(*state_->graph);

    if (!state_->clustering) {
        state_->last_error = "Leiden clustering failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "Found %zu clusters, modularity = %.4f\n",
                     state_->clustering->num_communities,
                     state_->clustering->modularity);
    }

    return true;
}

bool AllpairPipeline::RunSelfWaveCollapse() {
    if (!state_->graph || !state_->clustering) {
        state_->last_error = "Graph or clustering not available";
        return false;
    }

    SelfWaveCollapseConfig swc_config;
    swc_config.max_iterations = config_.swc_max_iterations;
    swc_config.convergence_threshold = config_.swc_convergence_threshold;
    swc_config.collapse_threshold = config_.swc_collapse_threshold;
    swc_config.gamma = config_.swc_gamma;

    SelfWaveCollapse swc(swc_config);
    state_->swc_result = swc.Refine(*state_->graph, *state_->clustering);

    if (!state_->swc_result) {
        state_->last_error = "Self-WaveCollapse refinement failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "SWC: %zu reads collapsed, %zu clusters\n",
                     state_->swc_result->stats.reads_collapsed,
                     state_->swc_result->num_clusters);
    }

    return true;
}

bool AllpairPipeline::SelectRepresentatives() {
    if (!state_->graph || !state_->clustering) {
        state_->last_error = "Graph or clustering not available";
        return false;
    }

    ClusterRepConfig rep_config;
    rep_config.method = ClusterRepConfig::Method::Medoid;
    rep_config.use_confidence_tiebreaker = true;

    ClusterRepSelector selector(rep_config);

    if (state_->swc_result) {
        state_->rep_result = selector.SelectWithConfidence(
            *state_->graph, *state_->clustering, *state_->swc_result);
    } else {
        state_->rep_result = selector.Select(*state_->graph, *state_->clustering);
    }

    if (!state_->rep_result) {
        state_->last_error = "Representative selection failed";
        return false;
    }

    if (config_.verbose) {
        std::fprintf(stderr, "Selected %zu representatives\n",
                     state_->rep_result->representatives.size());
    }

    return true;
}

}  // namespace llmap::self_interference
