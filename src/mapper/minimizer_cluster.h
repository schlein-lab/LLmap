// LLmap — Minimizer-overlap CPU read clustering (no FAISS/GPU).
//
// The Stage-1 self-interference clustering needs FAISS-GPU, absent on the
// (no-GPU) compute path. This is the CPU fallback that detects "clearly
// overlapping reads" (the operator's merge precondition) so they can be
// assembled into a mini-contig and mapped once. Reads from the same locus share
// most CANONICAL minimizers (canonical ⇒ strand-agnostic: a read and its
// reverse-complement cluster together; the per-read orientation is recorded).
//
// Operator constraints honoured here:
//  * NO smoothing — this only GROUPS reads; the consensus/layout step (lossless)
//    is separate and decides representative-read vs multi-allele layout.
//  * min-overlap NOT hardcoded — exposed as config (min_shared + min_overlap_bp),
//    conservative default (under-merge rather than build chimeric clusters).
//  * strand-aware — per-read reverse flag stored, never lost.
//  * a frequency cap drops repetitive minimizers (paralog/repeat over-merge guard).
//
// The WaveCollapse-ambiguity second gate (don't merge when placement is blurry)
// lives in the caller; this module is the cheap candidate detector.

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace llmap::mapper {

struct MinimizerClusterConfig {
    std::uint32_t k{15};               // k-mer length
    std::uint32_t w{10};               // minimizer window (in k-mers)
    std::uint32_t min_shared{5};       // min shared minimizers to link (conservative)
    std::uint32_t min_overlap_bp{300}; // min span covered by shared minimizers
    double        max_freq_frac{0.30}; // drop minimizers in > this fraction of reads
};

// Per-read cluster assignment + orientation relative to the cluster.
struct ReadCluster {
    std::uint32_t cluster_id{0};
    bool          reverse{false};   // read matched the cluster in reverse-complement
};

// Cluster reads by shared canonical-minimizer overlap (union-find). Returns one
// ReadCluster per input read; singletons get their own id. Strand-aware, no
// FAISS/GPU, no smoothing (grouping only).
[[nodiscard]] std::vector<ReadCluster> ClusterByMinimizers(
    std::span<const std::string> reads, const MinimizerClusterConfig& cfg = {});

// Group read indices by cluster id (clusters[c] = read indices).
[[nodiscard]] std::vector<std::vector<std::size_t>> GroupClusters(
    const std::vector<ReadCluster>& assignment);

}  // namespace llmap::mapper
