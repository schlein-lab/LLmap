// LLmap — AnchorRecord and supporting types.
//
// One AnchorRecord per (exonic region, transcript variant, source) tuple
// loaded into the AnchorStore. Anchors are the lookup units the
// Multi-Signal Fusion engine (Plan-Block 4.5) consults to compute
// L_sequence / L_modification / L_expression / L_pseudogene per read.
//
// An anchor is NOT a read — it's a reference sequence + metadata that
// reads are matched against. We deliberately mix DB-sourced anchors
// (GENCODE exons, IMGT germline segments, HPRC haplotype-specific
// annotations) and computed anchors (cluster representatives emerging
// from Stage-1 Self-WaveCollapse) under one type so downstream code
// doesn't care where each anchor came from — only `source` differs.
//
// Lifetime: AnchorRecord values live in AnchorStore::anchors_ for the
// duration of the run. Pointers handed out via Lookup APIs are stable
// until the next bulk-load call. Don't store anchor_id_idx values
// across reloads.

#pragma once

#include "core/transcript_kind.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace llmap::anchor {

// ===========================================================================
// AnchorSource — where this anchor came from. Stable wire-level integer
// codes so on-disk catalog files don't break when we add a source.
// ===========================================================================

enum class AnchorSource : std::uint16_t {
    Unknown = 0,
    Gencode,            ///< GENCODE v46+ GFF3 — primary annotation
    Mane,               ///< MANE Select (NCBI/EBI consensus)
    Refseq,             ///< NCBI RefSeq
    Imgt_GeneDb,        ///< IMGT/GENE-DB — IGH/IGK/IGL/TR germline
    Imgt_Hla,           ///< IPD-IMGT/HLA (CC-BY-ND: external reference only,
                        ///<                  never inlined into bundles)
    Fantom5_Cat,        ///< FANTOM5 CAT — lncRNA, eRNA
    ChessDb,            ///< CHESS — long-read merged transcripts
    Pangenome_PerHap,   ///< HPRC R2 per-haplotype annotation GFFs
    Branch_Bubble,      ///< BRANCH bubble bridge — dark-transcriptome
    Computed_Cluster,   ///< Stage-1 Self-WaveCollapse cluster representative
    Custom,             ///< User-supplied / project-specific
};

const char* AnchorSourceName(AnchorSource s) noexcept;
std::optional<AnchorSource> ParseAnchorSource(std::string_view s) noexcept;

// ===========================================================================
// ExonBoundary — a single exon-intron junction with splice-site evidence.
//
// Used both for (a) deciding where the k-mer-index may build
// junction-spanning k-mers without producing nonsense matches, and
// (b) feeding the Junction probability into the Wave-Collapse kernel.
// ===========================================================================

struct ExonBoundary {
    /// Position in transcript-coord (0-based) where the junction sits.
    /// pos_in_transcript ∈ [0, transcript_length).
    std::uint32_t pos_in_transcript{0};

    /// Genomic position of the donor base (last base of upstream exon).
    /// For back-splice (circRNA): donor > acceptor in this struct.
    std::uint64_t donor_genomic_pos{0};

    /// Genomic position of the acceptor base (first base of downstream
    /// exon). For canonical splice: acceptor > donor; for back-splice:
    /// acceptor < donor (defining property of a circle).
    std::uint64_t acceptor_genomic_pos{0};

    /// 2-bp motif at the donor (intron-side: usually "GT" for U2,
    /// "AT" for U12). Empty string ⇒ unknown / not yet scored.
    std::string donor_motif;

    /// 2-bp motif at the acceptor (intron-side: usually "AG" for U2,
    /// "AC" for U12).
    std::string acceptor_motif;

    /// Spliceosome class:
    ///   0 = U2 major (GT-AG / GC-AG)
    ///   1 = U12 minor (AT-AC)
    ///   2 = non-canonical (anything else)
    ///   3 = back-splice (back-spliced junction for circRNA)
    std::uint8_t spliceosome_class{0};

    /// Splice-site PWM scores ∈ [0,1]. 0 ⇒ not scored yet; populated by
    /// SpliceSiteDb::ScoreJunction in Block 2.5.
    float donor_score{0.0f};
    float acceptor_score{0.0f};

    /// Distance (bp) from acceptor to branch-point adenosine. -1 ⇒ not
    /// detected. Typical mammalian value is -20 to -50.
    std::int32_t branch_point_offset{-1};

    /// Length of the polypyrimidine tract before the acceptor AG.
    /// Useful as a side-feature when scoring questionable acceptor sites.
    std::uint8_t polypyrimidine_len{0};
};

// ===========================================================================
// AnchorRecord — the lookup unit.
// ===========================================================================

struct AnchorRecord {
    /// Stable identifier. Convention:
    ///   "<SRC>:<src-specific-id>[:<sub-id>]"
    ///   e.g. "GENCODE:ENST00000390557.4:exon1",
    ///        "IMGT:IGHG4*01:CH1",
    ///        "PANGEN:HG00329_hap1:G030290:CH1",
    ///        "CLUSTER:42"
    std::string anchor_id;

    AnchorSource source{AnchorSource::Unknown};

    /// Biological class. Filled by TranscriptKindClassifier
    /// (Block 2.5). For DB-sourced anchors usually known up front from
    /// the source's biotype; for cluster anchors inferred heuristically.
    core::TranscriptKind kind{core::TranscriptKind::Unknown};

    /// Free-form extension tag when kind == NovelUnclassified.
    std::optional<core::CustomKindTag> custom_kind;

    /// Raw nucleotide sequence the k-mer index hashes over.
    /// **Empty string allowed** — for IPD-IMGT/HLA where CC-BY-ND
    /// forbids inlining we keep only the metadata and the sequence
    /// is fetched on demand via the catalog adapter.
    std::string sequence;

    /// Genomic coordinates when the anchor is referable to a linear
    /// assembly. nullopt ⇒ anchor is annotation-/cluster-only.
    std::optional<std::string> ref_chrom;
    std::optional<std::int64_t> ref_start;   ///< 0-based inclusive
    std::optional<std::int64_t> ref_end;     ///< exclusive
    char strand{'.'};

    /// Owning transcript id when the anchor is an exon/UTR slice.
    std::string transcript_id;

    /// Host-gene id — set for snoRNA-in-intron, pseudogene-in-host,
    /// sterile germline (C-gene id), circRNA host etc.
    std::string host_gene_id;

    /// Exon-intron architecture (empty for non-spliced anchors).
    /// Ordered 5'→3' in transcript coordinates.
    std::vector<ExonBoundary> exon_boundaries;

    /// Free-form tags. Convention: snake_case, e.g.
    ///   "IGH", "constant_region", "CH1",
    ///   "Sgamma4_switch_region", "primary_transcript",
    ///   "polymorphic_pseudogene"
    std::vector<std::string> tags;

    /// License flag: 1 = inline-bundleable, 0 = external-reference only.
    /// Set to 0 for IPD-IMGT/HLA records per CC-BY-ND; bundle builders
    /// must respect this.
    std::int8_t license_inline_ok{1};

    /// Convenience predicates.
    [[nodiscard]] bool has_genomic_coords() const noexcept {
        return ref_chrom.has_value()
            && ref_start.has_value()
            && ref_end.has_value();
    }

    [[nodiscard]] bool is_spliced() const noexcept {
        return !exon_boundaries.empty();
    }

    [[nodiscard]] bool is_circular() const noexcept {
        for (const auto& b : exon_boundaries) {
            if (b.spliceosome_class == 3) return true;
        }
        return false;
    }
};

}  // namespace llmap::anchor
