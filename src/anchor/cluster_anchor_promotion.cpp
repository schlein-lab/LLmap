// LLmap — Cluster-Anchor-Promotion implementation.
//
// Single linear pass through the assignments vector, bucketing by
// cluster_id. For each bucket we pick the central read (cluster
// representative reported by Stage-1 EM) plus up to N-1 additional
// near-central reads as ranked by EM confidence.
//
// The emitted AnchorRecord fields:
//   source            = AnchorSource::Computed_Cluster
//   kind              = TranscriptKind::Unknown
//                        (the classifier (Block 2.5) refines later;
//                         we don't pre-judge here)
//   anchor_id         = "<ns>:<cluster_id>:<rank>"
//   sequence          = read_sequences[ranked_read_idx]
//   ref_chrom/start/end = unset (cluster anchors are by definition not
//                          referable to a single linear assembly)
//   transcript_id     = "" (no DB transcript)
//   host_gene_id      = "" (unknown by design)
//   exon_boundaries   = {}
//   tags              = {"computed_cluster",
//                         "cluster_size:<N>",
//                         "cluster_conf:<rounded-conf>"}
//
// The "cluster_size:N" tag lets the Fusion engine read the cluster
// cardinality without re-walking the assignments vector.

#include "anchor/cluster_anchor_promotion.h"

#include <algorithm>
#include <cstdio>
#include <unordered_map>
#include <vector>

namespace llmap::anchor {

namespace {

struct ClusterMembers {
    std::uint32_t cluster_id{0};
    std::uint32_t representative_read_idx{0};  // most-central
    float total_confidence{0.0f};
    std::uint32_t collapsed_count{0};
    std::vector<std::pair<std::uint32_t, float>>
        member_reads;  // (read_idx, confidence) sorted desc by conf
};

}  // namespace

PromotionStats PromoteClusterAnchors(
    const std::vector<ReadClusterAssignment>& assignments,
    std::span<const std::string> read_sequences,
    const ClusterPromotionConfig& cfg,
    AnchorStore& store) {

    PromotionStats stats;

    // Bucket by cluster_id.
    std::unordered_map<std::uint32_t, ClusterMembers> buckets;
    for (const auto& a : assignments) {
        auto& b = buckets[a.cluster_id];
        b.cluster_id = a.cluster_id;
        b.total_confidence += a.confidence;
        if (a.collapsed) ++b.collapsed_count;
        b.member_reads.emplace_back(a.read_idx, a.confidence);
        // Track the most-central read as the cluster representative.
        // We trust Stage-1's anchor_read pointer when available; fall
        // back to highest-confidence if anchor_read is 0 (unset).
        if (a.anchor_read != 0) {
            b.representative_read_idx = a.anchor_read;
        }
    }

    for (auto& [cid, b] : buckets) {
        ++stats.clusters_considered;

        // ----- Policy gates ------------------------------------------------
        if (b.member_reads.size() < cfg.min_cluster_size) {
            ++stats.clusters_skipped_too_small;
            continue;
        }
        const float mean_conf =
            b.total_confidence
            / static_cast<float>(b.member_reads.size());
        if (mean_conf < cfg.min_cluster_confidence) {
            ++stats.clusters_skipped_low_confidence;
            continue;
        }

        // Fall-back: if the bucket has no representative read recorded,
        // use the highest-confidence member.
        if (b.representative_read_idx == 0) {
            auto best = std::max_element(
                b.member_reads.begin(), b.member_reads.end(),
                [](const auto& l, const auto& r) {
                    return l.second < r.second;
                });
            if (best != b.member_reads.end()) {
                b.representative_read_idx = best->first;
            }
        }

        // Sort members by confidence descending for rank assignment.
        std::sort(b.member_reads.begin(), b.member_reads.end(),
                  [](const auto& l, const auto& r) {
                      return l.second > r.second;
                  });

        // Cap to max_anchors_per_cluster (unless subsample disabled).
        std::uint32_t to_emit = static_cast<std::uint32_t>(
            b.member_reads.size());
        if (cfg.subsample_large_clusters
            && to_emit > cfg.max_anchors_per_cluster) {
            to_emit = cfg.max_anchors_per_cluster;
        }

        // ----- Emit anchors ----------------------------------------------
        const std::size_t cluster_size = b.member_reads.size();
        const float rounded_conf =
            static_cast<float>(static_cast<int>(mean_conf * 100.0f))
            / 100.0f;

        for (std::uint32_t rank = 0; rank < to_emit; ++rank) {
            const std::uint32_t rid = b.member_reads[rank].first;
            if (rid >= read_sequences.size()) continue;

            AnchorRecord rec;
            rec.anchor_id = cfg.anchor_id_namespace
                + ":" + std::to_string(cid)
                + ":" + std::to_string(rank);
            rec.source = AnchorSource::Computed_Cluster;
            rec.sequence = read_sequences[rid];
            rec.tags.push_back("computed_cluster");
            rec.tags.push_back("cluster_size:"
                + std::to_string(cluster_size));
            // Tag the confidence with a 2-decimal-place rounded value
            // so downstream filtering can do exact-string matching
            // without floating-point gymnastics.
            char conf_buf[32];
            std::snprintf(conf_buf, sizeof(conf_buf),
                          "cluster_conf:%.2f", rounded_conf);
            rec.tags.emplace_back(conf_buf);
            // Mark the rank-0 anchor as the representative so callers
            // can prefer it when multiple anchors of the same cluster
            // match a read.
            if (rank == 0) rec.tags.emplace_back("cluster_representative");

            store.AddAnchor(std::move(rec));
            ++stats.anchors_emitted;
        }

        ++stats.clusters_promoted;
    }

    // After bulk insertion, rebuild secondary indices for consistency.
    store.Reindex();
    return stats;
}

}  // namespace llmap::anchor
