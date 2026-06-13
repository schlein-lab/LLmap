// LLmap — Assemble-then-map orchestrator (Line A glue).
//
// Ties the CPU clustering fallback to the contig assembly: cluster reads by
// shared canonical-minimizer overlap (minimizer_cluster), group them, and
// assemble each cluster into a probe + lossless layout (cluster_consensus). The
// output is one ConsensusContig per cluster, ready for the bucket-driven
// map→propagate loop.
//
// KEY CONTRACT: AssembleConsensus reports member.read_idx LOCAL to the cluster
// subset it was handed. This orchestrator REMAPS every member.read_idx back to
// the GLOBAL index into the original `reads` span, so the downstream propagation
// addresses real reads directly. Strand (ReadCluster.reverse) is threaded into
// AssembleConsensus so reverse-strand members are oriented into the contig frame
// and recorded — never dropped (lossless, Line A).

#pragma once

#include <span>
#include <string>
#include <vector>

#include "mapper/cluster_consensus.h"
#include "mapper/minimizer_cluster.h"

namespace llmap::mapper {

// Assemble all clusters of `reads`. Returns one ConsensusContig per cluster
// (singletons included), with member.read_idx in GLOBAL `reads` coordinates.
[[nodiscard]] std::vector<ConsensusContig> AssembleClusters(
    std::span<const std::string> reads,
    const MinimizerClusterConfig& cluster_cfg = {},
    const ConsensusConfig& consensus_cfg = {});

}  // namespace llmap::mapper
