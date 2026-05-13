#include "self_interference/similarity_graph.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>

namespace llmap::self_interference {

std::vector<Edge> SimilarityGraph::ToEdgeList() const {
    std::vector<Edge> edges;
    edges.reserve(col_indices_.size());

    for (size_t src = 0; src < num_nodes_; ++src) {
        const size_t start = row_offsets_[src];
        const size_t end = row_offsets_[src + 1];
        for (size_t i = start; i < end; ++i) {
            edges.push_back({
                static_cast<uint32_t>(src),
                col_indices_[i],
                weights_[i]
            });
        }
    }

    return edges;
}

std::string SimilarityGraph::ToEdgeListString() const {
    std::ostringstream oss;
    oss << "# SimilarityGraph: " << num_nodes_ << " nodes, "
        << col_indices_.size() << " edges\n";
    oss << "# source target weight\n";

    for (size_t src = 0; src < num_nodes_; ++src) {
        const size_t start = row_offsets_[src];
        const size_t end = row_offsets_[src + 1];
        for (size_t i = start; i < end; ++i) {
            oss << src << " " << col_indices_[i] << " " << weights_[i] << "\n";
        }
    }

    return oss.str();
}

bool SimilarityGraph::Save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    if (!file) return false;

    const char magic[] = "LLMG";
    file.write(magic, 4);

    const uint32_t version = 1;
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));

    const auto n = static_cast<uint64_t>(num_nodes_);
    const auto e = static_cast<uint64_t>(col_indices_.size());
    file.write(reinterpret_cast<const char*>(&n), sizeof(n));
    file.write(reinterpret_cast<const char*>(&e), sizeof(e));

    const uint8_t sym = is_symmetric_ ? 1 : 0;
    file.write(reinterpret_cast<const char*>(&sym), sizeof(sym));

    file.write(reinterpret_cast<const char*>(row_offsets_.data()),
               row_offsets_.size() * sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(col_indices_.data()),
               col_indices_.size() * sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(weights_.data()),
               weights_.size() * sizeof(float));

    return file.good();
}

std::unique_ptr<SimilarityGraph> SimilarityGraph::Load(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return nullptr;

    char magic[4];
    file.read(magic, 4);
    if (std::memcmp(magic, "LLMG", 4) != 0) return nullptr;

    uint32_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) return nullptr;

    uint64_t n, e;
    file.read(reinterpret_cast<char*>(&n), sizeof(n));
    file.read(reinterpret_cast<char*>(&e), sizeof(e));

    uint8_t sym;
    file.read(reinterpret_cast<char*>(&sym), sizeof(sym));

    auto graph = std::unique_ptr<SimilarityGraph>(new SimilarityGraph());
    graph->num_nodes_ = static_cast<size_t>(n);
    graph->is_symmetric_ = (sym != 0);

    graph->row_offsets_.resize(n + 1);
    graph->col_indices_.resize(e);
    graph->weights_.resize(e);

    file.read(reinterpret_cast<char*>(graph->row_offsets_.data()),
              graph->row_offsets_.size() * sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(graph->col_indices_.data()),
              graph->col_indices_.size() * sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(graph->weights_.data()),
              graph->weights_.size() * sizeof(float));

    if (!file.good()) return nullptr;

    graph->build_stats_.num_nodes = graph->num_nodes_;
    graph->build_stats_.num_edges = graph->col_indices_.size();
    if (!graph->weights_.empty()) {
        graph->build_stats_.max_weight = *std::max_element(
            graph->weights_.begin(), graph->weights_.end());
        graph->build_stats_.min_weight = *std::min_element(
            graph->weights_.begin(), graph->weights_.end());
    }
    if (graph->num_nodes_ > 0) {
        graph->build_stats_.avg_degree =
            static_cast<float>(graph->col_indices_.size()) /
            static_cast<float>(graph->num_nodes_);
    }

    return graph;
}

}  // namespace llmap::self_interference
