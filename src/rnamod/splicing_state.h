// LLmap — SplicingState: orthogonal classification of how a read was
// (or wasn't) spliced.
//
// AlignmentStatus tells us the *kind* of match (Mapped vs Tentative
// vs Chimeric…), TranscriptKind tells us the *biological class*
// (mRNA vs lncRNA vs miRNA…). SplicingState is the third axis:
// "what did the spliceosome do here". Important because the same
// transcript can land in many states across reads of the same library
// (canonical vs intron-retained vs lariat vs recursive vs trans-splice
// vs back-splice), and we want to track all of them losslessly.

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace llmap::splicing {

enum class SplicingState : std::uint16_t {
    Unknown = 0,
    Canonical,              ///< all introns spliced (mature)
    Unspliced,              ///< no introns spliced (pre-mRNA / nascent)
    PartiallySpliced,       ///< some introns retained, others spliced
    IntronRetained,         ///< single specific intron retained (regulatory)
    Lariat,                 ///< lariat intermediate captured (rare)
    RecursiveSpliced,       ///< long intron broken by RS-exon (TITIN, DMD…)
    TransSpliced,           ///< two pre-mRNAs joined
    AlternativeCassetteIn,  ///< cassette exon included
    AlternativeCassetteOut, ///< cassette exon skipped
    Alt5ss,                 ///< alternative 5' splice site
    Alt3ss,                 ///< alternative 3' splice site
    MutuallyExclusive,      ///< A xor B (not both, not neither)
    BackSplicedCircular,    ///< circRNA generation
    HalfSplicedCotrans,     ///< co-transcriptional half-spliced
    NovelSplicingState,     ///< doesn't fit a known category
};

const char* SplicingStateName(SplicingState s) noexcept;

/// Per-read inference. A single read can manifest more than one state
/// (e.g. one retained intron + one recursive splice), so we record a
/// dominant plus a list of additional states.
struct SplicingStateInference {
    SplicingState dominant{SplicingState::Unknown};
    std::vector<SplicingState> additional;
    float confidence{0.0f};                              ///< [0,1]
    std::vector<std::uint32_t> retained_intron_indices;
    std::optional<std::uint32_t> recursive_splice_site_pos;
    std::optional<std::pair<std::string, std::string>>
        trans_splice_donor_acceptor_genes;
};

/// Minimal junction record consumed by the classifier. Same shape as the
/// AnchorRecord::exon_boundaries entry but free of the
/// llmap::anchor namespace dependency so this header can be linked from
/// the lossless aggregator without pulling in the anchor library.
struct ObservedJunction {
    std::uint64_t donor_genomic_pos{0};
    std::uint64_t acceptor_genomic_pos{0};
    std::uint8_t spliceosome_class{2};   ///< 0=U2 / 1=U12 / 2=non-canon / 3=back
    bool is_annotated{false};            ///< matched JunctionDb evidence
    bool is_retained{false};             ///< intron retained at this position
    std::uint32_t pos_in_read{0};        ///< 0-based read offset
};

class SplicingStateClassifier {
public:
    /// Classify a read from its observed junction list.
    [[nodiscard]] SplicingStateInference Classify(
        std::span<const ObservedJunction> junctions) const;
};

}  // namespace llmap::splicing
