// LLmap — SegDup catalog interface (stub for the IGH-mode classifier).
//
// The full SegDup catalog loader is being implemented in `src/catalog/` by a
// parallel agent. To unblock the IGH-mode upstream-anchor classifier in
// `src/classify/`, this header declares the *minimum surface* the classifier
// needs from the catalog: read-only lookup of an entry's mapping-strategy
// hints and its diagnostic promoter signature.
//
// When the real loader lands, it must satisfy this interface (either by
// implementing `ISegDupCatalog` directly or via a thin adapter). Nothing in
// `src/classify/` should depend on the concrete catalog implementation.
//
// Naming: `iface` (not `interface`) because MSVC reserves the latter.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace llmap::classify::igh {

/// @brief Mapping-strategy hints from `catalog/curated/<locus>.json`'s
///        `mapping_strategy.primary` block.
///
/// Only the fields consumed by the upstream-anchor classifier are exposed.
/// All fields are post-validation (i.e. the catalog loader has already
/// enforced the JSON-schema ranges).
struct PrimaryMappingHints {
    /// k-mer size for the primary chain (catalog: kmer_size).
    std::uint8_t kmer_size{25};

    /// Maximum k-mer mismatches tolerated (catalog: max_mismatch).
    std::uint8_t max_mismatch{2};

    /// How far up-/downstream of `include_flanking_anchor` to extend the
    /// classification window, in base pairs (catalog: include_flanking_bp).
    /// 0 disables flanking — the classifier then operates on the CH window
    /// implied by the catalog's anchor only.
    std::uint32_t include_flanking_bp{0};

    /// Free-form anchor name (catalog: include_flanking_anchor); e.g.
    /// "CH1_start". The IGH-mode classifier interprets this as a 5' edge to
    /// extend from; unknown anchors fall back to "treat read start as anchor"
    /// in the test path.
    std::string include_flanking_anchor;
};

/// @brief Catalog promoter-signature block.
///
/// Mirrors `diagnostic_features.promoter_signature` in the curated JSON. The
/// motif strings may contain a literal `...` (three dots) meaning
/// "ordered-subsequence within the window"; e.g.
/// `"CGGTTCTT...GTCTATCTGCGATGG"` requires the prefix then the suffix to
/// appear in 5'->3' order with no constraint on the spacer length (within the
/// `window_relative_to_anchor` span).
///
/// `canonical_motif` / `duplicate_motif` are `std::optional` so legacy
/// entries that only fill `expected_motif` can be detected and rejected by
/// the classifier (the binary canonical-vs-chimdup discriminant requires
/// *both* motifs to be set).
struct PromoterSignature {
    /// Window relative to the anchor, e.g. "-100..0" (5' UTR / promoter).
    std::string window_relative_to_anchor;

    /// Anchor name, e.g. "IGHG4_CH1_start".
    std::string anchor;

    /// Canonical-allele motif (catalog: canonical_motif). May contain `...`.
    std::optional<std::string> canonical_motif;

    /// Duplicate-allele motif (catalog: duplicate_motif). May contain `...`.
    std::optional<std::string> duplicate_motif;
};

/// @brief Read-only view of a catalog entry, restricted to the fields the
/// IGH-mode classifier reads.
///
/// Owning storage lives in the catalog implementation; the view is valid for
/// as long as the catalog object is alive. The view never aliases mutable
/// state — concurrent reads are safe.
struct CatalogEntryView {
    /// Stable locus id (catalog: locus_id), e.g. "IGHG4_chimdup_tandem".
    std::string locus_id;

    /// Haplotype-class tag (catalog: haplotype_class), e.g.
    /// "IGHG4_chimdup_homozygous". The classifier emits this verbatim as
    /// `haplotype_class_call` when the read matches.
    std::string haplotype_class;

    /// Primary mapping hints — controls flanking-window size + k-mer params.
    PrimaryMappingHints mapping_primary;

    /// Promoter signature — drives the binary canonical/chimdup motif test.
    /// Absent if the entry has no `diagnostic_features.promoter_signature`.
    std::optional<PromoterSignature> promoter_signature;
};

/// @brief Read-only interface to a SegDup catalog.
///
/// The IGH-mode classifier depends on this interface, *not* on the concrete
/// loader. The real loader (in `src/catalog/`) will implement this; tests in
/// `tests/classify/` use an in-memory mock.
class ISegDupCatalog {
public:
    virtual ~ISegDupCatalog() = default;

    /// @brief Look up a catalog entry by reference target + 0-based position.
    /// @return The entry view, or `std::nullopt` if no curated entry covers
    /// `(target_id, pos)`.
    [[nodiscard]] virtual std::optional<CatalogEntryView> LookupByPosition(
        std::string_view target_id, std::uint64_t pos) const = 0;

    /// @brief Look up a catalog entry by its haplotype-class tag.
    /// @return The entry view, or `std::nullopt` if the tag is unknown.
    ///
    /// Used by tests and by re-classification paths where the caller already
    /// knows the haplotype hypothesis and wants its diagnostic features.
    [[nodiscard]] virtual std::optional<CatalogEntryView>
    LookupByHaplotypeClass(std::string_view haplotype_class) const = 0;
};

}  // namespace llmap::classify::igh
