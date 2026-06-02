// LLmap — Transcript-Mode output schema: BAM tags + Parquet sidecars.
//
// The existing src/output/bam_writer.* + parquet_writer.* serialise
// AlignmentRecord values. This module adds Transcript-Mode-specific
// tag generators + sidecar row structs that the writers consume when
// `include_transcript_tags=true` is set on their config.
//
// BAM tags emitted (matching minimap2 / STAR conventions for ecosystem
// tool compatibility):
//
//   XS — A — transcript strand: '+' / '-' / '?'
//   jI — B,i — junction list: pairs of (donor_ref_pos, acceptor_ref_pos)
//   jM — B,f — junction confidence per junction (float, not binary)
//   TI — Z — best-matching transcript_id
//          (GENCODE:ENST… / MANE:… / IMGT:… / PANGEN:HG00329_hap1:… /
//           CLUSTER:42 / NOVEL)
//   XK — Z — category tag: MAPPED / MAPPED_STERILE / MAPPED_PREMRNA /
//          MAPPED_CIRCULAR / PARTIAL / CHIM_I / CHIM_X / CHIM_V /
//          DARK / TENT / UNMAPPED
//   XC — Z — cluster_id from Stage-1 (always set; "-" when no cluster)
//   XA — Z — anchor source list, comma-separated
//          ("GENCODE,MANE,IMGT,PANGEN_HG00329,CLUSTER")
//   XB — Z — BRANCH bubble IDs when integrated, else absent
//   XQ — Z — legacy MAPQ value (the lossless contract recalibrates
//          MAPQ; XQ preserves the upstream seed-extend value so
//          downstream pipelines can transparently migrate)
//   XM — Z — RNA modification calls (comma-sep: kind:pos:conf)
//   XF — Z — splicing-state name (Block 2.6 SplicingState enum string)
//
// Parquet sidecars (one per concept; each row independent):
//
//   <out>.align.parquet        — per-read main row (already exists in
//                                  parquet_writer; this commit just
//                                  adds new columns)
//   <out>.exon_coverage.parquet — per-(transcript, exon) row
//   <out>.cluster_summary.parquet — per-cluster row
//   <out>.chimeric.parquet      — per-chimeric-read row
//   <out>.rna_modifications.parquet — per-mod-call row
//   <out>.splicing_states.parquet   — per-read splicing-state row

#pragma once

#include "core/alignment_record.h"
#include "core/transcript_kind.h"

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace llmap::output::transcript_schema {

// ===========================================================================
// BAM tag builders — pure helpers; no SAM/BAM dep, just string concat.
// ===========================================================================

/// Build the XS tag value from a strand char.
[[nodiscard]] char XsTag(char strand) noexcept;

/// Encode the jI tag value as "<donor1>,<acc1>,<donor2>,<acc2>,…".
/// Empty when no junctions.
[[nodiscard]] std::string JiTag(
    std::span<const std::pair<std::uint64_t, std::uint64_t>> junctions);

/// Encode the jM tag as comma-separated floats.
[[nodiscard]] std::string JmTag(std::span<const float> per_junction_conf);

/// XK status string for a given AlignmentStatus (uses the existing
/// AlignmentStatusName helper).
[[nodiscard]] const char* XkTag(AlignmentStatus s) noexcept;

/// XC tag value (cluster_id). Returns "-" when cluster_id == 0.
[[nodiscard]] std::string XcTag(std::uint32_t cluster_id);

/// XA tag: comma-separated AnchorSource names.
[[nodiscard]] std::string XaTag(std::span<const std::string> sources);

/// XQ tag: legacy MAPQ value as decimal string ("-1" if unknown).
[[nodiscard]] std::string XqTag(std::int32_t legacy_mapq);

/// XM tag: rna-modification call list, comma-separated
/// "kind:pos:conf". Empty when no calls.
struct ModCallView {
    std::string_view kind;
    std::uint32_t pos;
    float confidence;
};
[[nodiscard]] std::string XmTag(std::span<const ModCallView> calls);

/// XF tag: splicing-state name (Block 2.6 SplicingStateName).
[[nodiscard]] std::string XfTag(std::string_view splicing_state_name);

// ===========================================================================
// Parquet sidecar row structs.
// ===========================================================================

struct ExonCoverageRow {
    std::string transcript_id;
    std::uint32_t exon_index{0};       // 1-based per GENCODE
    std::string  ref_name;
    std::uint64_t ref_start{0};
    std::uint64_t ref_end{0};
    std::uint64_t reads_covering{0};
    double mean_depth{0.0};
    double locus_median_depth{0.0};
    double frac_locus_median{0.0};
    bool is_uniformity_flag{false};    // true when frac < 0.10
};

struct ClusterSummaryRow {
    std::uint32_t cluster_id{0};
    std::uint32_t n_reads{0};
    std::string representative_seq_hash;
    std::string matched_anchor_sources;   // comma-sep
    std::string dominant_transcript_id;
    bool is_dark_cluster{false};
};

struct ChimericRow {
    std::string read_id;
    char chimeric_kind{'.'};               // 'I' / 'X' / 'V'
    std::string parts_anchor_ids;          // comma-sep
    std::string parts_ref_pos;             // comma-sep "chrom:pos"
    std::string parts_prob;                // comma-sep floats
    bool vdj_class_switch_detected{false};
    std::uint64_t genomic_distance_bp{0};
};

struct RnaModificationRow {
    std::string read_id;
    std::string chrom;
    std::uint64_t pos{0};
    std::string kind;
    float confidence{0.0f};
    char source_tool{'?'};
    std::string custom_label;              // empty unless NovelMod
};

struct SplicingStateRow {
    std::string read_id;
    std::string dominant_state;
    std::string additional_states;         // comma-sep
    std::uint32_t n_retained_introns{0};
    std::optional<std::uint32_t> recursive_splice_site_pos;
    std::string trans_donor_gene;
    std::string trans_acceptor_gene;
};

}  // namespace llmap::output::transcript_schema
