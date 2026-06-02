// LLmap — Cluster-Anchor-Promotion: turn Stage-1 read clusters into
// computed anchors that join the AnchorStore alongside the DB-sourced
// entries (GENCODE, MANE, IMGT, …).
//
// Stage 1 of the LLmap pipeline runs Self-WaveCollapse over the input
// reads (Foundation-Model embeddings → FAISS k-NN → Leiden clustering
// → intra-cluster EM). The output is a list of (read_idx, cluster_id,
// confidence, anchor_read) tuples. This module takes that output PLUS
// the original read sequences and emits AnchorRecord values into the
// AnchorStore so the Multi-Signal Fusion engine can treat
// cluster-internal anchors as a first-class source.
//
// Two reasons this matters:
//
//   1. Dark transcriptome — reads that don't match any DB annotation
//      still form clusters when they share sequence. A cluster's
//      representative is a perfectly valid anchor for the OTHER reads
//      in that cluster, even if no GENCODE/IMGT/HPRC entry covers
//      what's being transcribed.
//
//   2. Paralog disambiguation under coverage signal — when two
//      sequence-identical paralog copies have asymmetric expression,
//      the cluster structure preserves the asymmetry. Promoting the
//      cluster representatives lets the Wave-Collapse kernel exploit
//      that asymmetry without us needing per-cluster DB anchors.
//
// This module does NOT re-implement Self-WaveCollapse — it takes the
// existing pipeline's output via a small ReadAssignment struct (mirror
// of what the self_interference module produces) and walks it. Keeps
// us decoupled from any churn in the Stage-1 internals.

#pragma once

#include "anchor/anchor_record.h"
#include "anchor/anchor_store.h"

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::anchor {

/// Cluster-promotion policy. Defaults match the values from
/// AnchorStore::ClusterPromotionPolicy but kept here as a copy so this
/// header doesn't need to drag in the full AnchorStore declaration for
/// callers that only want to compute the promotion.
struct ClusterPromotionConfig {
    std::uint32_t min_cluster_size{3};        ///< singletons not promoted
    float         min_cluster_confidence{0.7f};
    bool          subsample_large_clusters{true};
    std::uint32_t max_anchors_per_cluster{5};
    /// Optional namespace prefix for the synthesised anchor_id, default
    /// "CLUSTER". The final id is "CLUSTER:<cluster_id>:<rank>".
    std::string   anchor_id_namespace{"CLUSTER"};
};

/// Lightweight mirror of what self_interference::SelfWaveCollapseResult
/// produces per read. We declare it here so this module can be linked
/// from tests + downstream consumers without dragging in the FAISS/
/// embedding dep chain. The Stage-1 pipeline converts its native
/// struct → this one before calling PromoteClusterAnchors.
struct ReadClusterAssignment {
    std::uint32_t read_idx{0};        ///< index into the input read array
    std::uint32_t cluster_id{0};
    float         confidence{0.0f};   ///< Stage-1 EM convergence weight
    bool          collapsed{false};   ///< true if read converged
    std::uint32_t anchor_read{0};     ///< most-central read in cluster
                                       ///<   (i.e. the cluster representative)
};

/// Stats returned by PromoteClusterAnchors so the caller can log them
/// into the lossless summary.
struct PromotionStats {
    std::uint32_t clusters_considered{0};
    std::uint32_t clusters_promoted{0};
    std::uint32_t anchors_emitted{0};
    std::uint32_t clusters_skipped_too_small{0};
    std::uint32_t clusters_skipped_low_confidence{0};
};

/// Promote cluster representatives into the AnchorStore.
///
/// Walks `assignments`, groups reads by cluster_id, and for each cluster
/// that passes the policy thresholds emits up to max_anchors_per_cluster
/// AnchorRecord values into `store`. Anchor IDs follow the pattern
/// "<namespace>:<cluster_id>:<rank-from-most-central>".
///
/// Returns aggregated statistics.
PromotionStats PromoteClusterAnchors(
    const std::vector<ReadClusterAssignment>& assignments,
    std::span<const std::string> read_sequences,
    const ClusterPromotionConfig& cfg,
    AnchorStore& store);

}  // namespace llmap::anchor
