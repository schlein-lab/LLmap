// LLmap — AllpairPipeline internal state (private header).
//
// Shared by the split implementation files.

#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "io/fastq_reader.h"
#include "self_interference/allpair_pipeline.h"
#include "self_interference/cluster_rep.h"
#include "self_interference/faiss_wrapper.h"
#include "self_interference/self_wavecollapse.h"
#include "self_interference/similarity_graph.h"

namespace llmap::self_interference {

// Internal state shared across pipeline stages
struct AllpairPipeline::InternalState {
    // Loaded reads
    std::vector<io::FastqRecord> reads;

    // Computed embeddings (flattened: [num_reads * embedding_dim])
    std::vector<float> embeddings;

    // FAISS index
    std::unique_ptr<FaissIndex> faiss_index;

    // Similarity graph
    std::unique_ptr<SimilarityGraph> graph;

    // Clustering result
    std::unique_ptr<ClusteringResult> clustering;

    // Self-WaveCollapse result
    std::unique_ptr<SelfWaveCollapseResult> swc_result;

    // Representative selection result
    std::unique_ptr<ClusterRepResult> rep_result;

    // Timing
    std::chrono::steady_clock::time_point start_time;

    // Error message
    std::string last_error;
};

}  // namespace llmap::self_interference
