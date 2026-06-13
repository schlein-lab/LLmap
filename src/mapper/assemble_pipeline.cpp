// LLmap — Assemble-then-map orchestrator implementation.

#include "mapper/assemble_pipeline.h"

#include <cstdint>
#include <string>
#include <vector>

namespace llmap::mapper {

std::vector<ConsensusContig> AssembleClusters(
    std::span<const std::string> reads, const MinimizerClusterConfig& cluster_cfg,
    const ConsensusConfig& consensus_cfg) {
    std::vector<ConsensusContig> out;
    if (reads.empty()) return out;

    const std::vector<ReadCluster> assignment =
        ClusterByMinimizers(reads, cluster_cfg);
    const std::vector<std::vector<std::size_t>> groups = GroupClusters(assignment);

    out.reserve(groups.size());
    for (const auto& group : groups) {
        if (group.empty()) continue;
        // Gather this cluster's reads + strand flags (local-indexed subset).
        std::vector<std::string> sub;
        std::vector<std::uint8_t> rev;
        sub.reserve(group.size());
        rev.reserve(group.size());
        for (const std::size_t gi : group) {
            sub.push_back(reads[gi]);
            rev.push_back(assignment[gi].reverse ? 1u : 0u);
        }
        ConsensusContig contig = AssembleConsensus(sub, consensus_cfg, rev);
        // Remap member.read_idx: local subset index → global reads index.
        for (auto& m : contig.members) {
            m.read_idx = group[m.read_idx];
        }
        out.push_back(std::move(contig));
    }
    return out;
}

}  // namespace llmap::mapper
