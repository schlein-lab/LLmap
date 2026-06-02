// LLmap — AnchorStore: multi-source anchor catalog for Transcript-Mode.
//
// In-memory lookup index over AnchorRecord values from many sources:
//   GENCODE GFF3, MANE, RefSeq, IMGT/GENE-DB, IPD-IMGT/HLA (external-ref
//   only), FANTOM5-CAT, ChessDB, HPRC per-haplotype, BRANCH bubbles,
//   computed cluster anchors.
//
// Lookups exposed:
//   - by stable anchor_id
//   - by region (chrom + half-open [start,end))
//   - by tag ("IGH", "constant_region", "Sgamma4_switch_region", …)
//   - by transcript_id
//   - by source (filter by AnchorSource enum)
//
// Threadsafety: the store is read-only after the initial bulk Load
// calls have completed. Multiple reader threads may call any Lookup*
// method concurrently. Mutation methods (Load*, RegisterClusterAnchor)
// MUST run single-threaded; the caller is responsible for that
// ordering. We don't lock internally because the alignment hot loop
// would pay for synchronisation that isn't needed once the store is
// frozen.
//
// Memory cost: each anchor carries its sequence inline. For GENCODE v46
// (~ 250k transcripts × ~10 exons each × ~250 bp avg = ~600 MB) that's
// the dominant footprint. Acceptable on modern HPC hosts; if it ever
// blows budget we can switch to mmap-backed sequence storage.

#pragma once

#include "anchor/anchor_record.h"

#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::anchor {

/// Outcome of a single Load* call; aggregated by callers.
struct LoadStatus {
    bool ok{false};
    std::size_t records_loaded{0};
    std::size_t records_skipped{0};
    std::string error;
};

/// Policy for cluster-anchor promotion (Block 4).
/// Kept here so the AnchorStore type knows about it; the actual
/// promotion logic lives in src/anchor/cluster_anchor_promotion.{h,cpp}.
struct ClusterPromotionPolicy {
    std::uint32_t min_cluster_size{3};       ///< singletons not promoted
    float min_cluster_confidence{0.7f};
    bool subsample_large_clusters{true};
    std::uint32_t max_anchors_per_cluster{5};
};

class AnchorStore {
public:
    AnchorStore() = default;

    // ----- Load APIs (single-threaded; see threadsafety note above) -----

    /// Load a GENCODE GFF3 (.gz or .gff3). Filters to feature=exon and
    /// builds AnchorRecord(s) per (transcript_id, exon) tuple. The
    /// caller-supplied reference FASTA is required if `with_sequence`
    /// is true; otherwise sequences are left empty and filled lazily
    /// when the k-mer index needs them.
    LoadStatus LoadGencodeGff(const std::filesystem::path& gff,
                              const std::filesystem::path& ref_fa,
                              bool with_sequence = true);

    /// Load MANE Select summary TSV. Mostly cross-references GENCODE
    /// transcript IDs as canonical; we add tags={"mane_select"} to the
    /// corresponding GENCODE anchors rather than duplicating them.
    LoadStatus LoadMane(const std::filesystem::path& summary_tsv);

    /// Load IMGT/GENE-DB FASTA (IGH/IGK/IGL/TR germline segments).
    /// Each segment becomes an AnchorRecord with source=Imgt_GeneDb.
    LoadStatus LoadImgtGeneDb(const std::filesystem::path& fasta);

    /// Load HPRC R2 per-haplotype annotation GFFs from a directory tree.
    /// Walks subdirs matching the HPRC layout
    ///   <root>/<sample_id>/annotation/*.gff3.gz
    /// and pulls exon records into anchors tagged with the sample/hap id.
    LoadStatus LoadPangenomeAnnotations(
        const std::filesystem::path& hprc_root,
        const std::vector<std::string>& sample_ids);

    /// Import BRANCH bubble BED ([[branch_tooling_state]]) as anchors.
    /// Each bubble becomes one AnchorRecord with source=Branch_Bubble,
    /// kind=NovelUnclassified, tags=["branch_bubble", "vaf_<...>"].
    /// Dark-transcriptome bridge: BRANCH catches what no annotation does.
    LoadStatus ImportBranchBubbles(const std::filesystem::path& bed);

    /// Add an externally-built anchor (used by Stage-1 cluster promotion
    /// and by tests). Returns the assigned index into anchors_.
    std::uint32_t AddAnchor(AnchorRecord rec);

    // ----- Read APIs (threadsafe after Load* has returned) --------------

    /// All anchors in load order. Stable while no mutation pending.
    [[nodiscard]] const std::vector<AnchorRecord>& anchors() const noexcept {
        return anchors_;
    }
    [[nodiscard]] std::size_t size() const noexcept { return anchors_.size(); }

    /// Lookup by stable anchor_id. nullptr if not found.
    [[nodiscard]] const AnchorRecord* ById(std::string_view id) const;

    /// All anchors with `tag` present in their tags vector.
    /// Returns indices into anchors_ (cheap; caller dereferences).
    [[nodiscard]] std::vector<std::uint32_t>
    ByTag(std::string_view tag) const;

    /// All anchors whose half-open genomic interval intersects
    /// [start, end) on `chrom`. Linear over per-chrom sorted index;
    /// for genome-wide queries use ByTag instead.
    [[nodiscard]] std::vector<std::uint32_t>
    ByRegion(std::string_view chrom,
             std::int64_t start,
             std::int64_t end) const;

    /// All anchors belonging to `transcript_id` (cross-source — could
    /// be GENCODE + MANE + RefSeq for the same transcript).
    [[nodiscard]] std::vector<std::uint32_t>
    ByTranscriptId(std::string_view transcript_id) const;

    /// Iterate every anchor matching a predicate. Returns indices.
    void ForEach(const std::function<bool(const AnchorRecord&)>& pred,
                  const std::function<void(std::uint32_t,
                                           const AnchorRecord&)>& cb) const;

    /// Count by source — telemetry / sanity-check helper.
    [[nodiscard]] std::size_t CountBySource(AnchorSource s) const;

    // ----- Maintenance -------------------------------------------------

    /// Rebuild the secondary indices (by_id_, by_chrom_, by_tag_,
    /// by_transcript_) from scratch. Called automatically at the end
    /// of each successful Load* call.
    void Reindex();

    /// Clear everything; mainly for tests.
    void Clear();

private:
    std::vector<AnchorRecord> anchors_;

    // Secondary indices — all populated by Reindex().
    std::unordered_map<std::string, std::uint32_t>          by_id_;
    std::unordered_map<std::string,
                       std::vector<std::uint32_t>>          by_chrom_;
    std::unordered_map<std::string,
                       std::vector<std::uint32_t>>          by_tag_;
    std::unordered_map<std::string,
                       std::vector<std::uint32_t>>          by_transcript_;
};

}  // namespace llmap::anchor
