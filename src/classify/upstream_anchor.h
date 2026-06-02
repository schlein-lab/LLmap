// LLmap — IGH-mode upstream-anchor classifier.
//
// Canonical-vs-chimdup IGHG4 alleles are *not* discriminable by CH1 SNPs
// alone: 4 of the 5 known CH1 SNPs are polymorphic in the canonical
// background. The robust discriminant is the ~100 bp upstream of CH1, where
// canonical haps carry `CGGTTCTT...GTCTATCTGCGATGG` and chimdup haps carry
// `GGGTTCTT...AACTGTCCGCGAGG` (the chimdup motif is IGHG1-recombinant; see
// `[[ighg4_upstream_signatures]]`).
//
// This classifier:
//   1. Looks up the catalog entry for a mapping (by position or by
//      pre-supplied haplotype-class hypothesis).
//   2. Extends the classification window upstream by
//      `mapping_strategy.primary.include_flanking_bp` (default 200 bp).
//   3. Performs a k-mer mismatch scan over flanking + CH1 (returns int).
//   4. Tests the *exact* canonical / duplicate promoter motifs from
//      `diagnostic_features.promoter_signature` against the upstream window.
//
// All motif/anchor strings come from the catalog — *nothing* is hard-coded
// here.
//
// Concurrency: the classifier is stateless. `Classify` is reentrant as long
// as the supplied `ISegDupCatalog` is thread-safe for reads (the curated
// catalogs are immutable after load, so this is the expected case).

#pragma once

#include "classify/segdup_catalog_iface.h"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace llmap::classify::igh {

/// @brief Outcome of the promoter-motif test.
enum class PromoterMotifMatch : std::uint8_t {
    /// Neither canonical nor duplicate motif matched the upstream window.
    None = 0,
    /// Canonical motif matched, duplicate did not.
    Canonical = 1,
    /// Duplicate motif matched, canonical did not.
    Chimdup = 2,
    /// Both motifs matched — ambiguous read (e.g. chimeric spanning read).
    Mixed = 3,
};

/// @brief String form of `PromoterMotifMatch` (stable, lowercase tokens
/// suitable for serialization in BAM tags / Parquet columns).
[[nodiscard]] std::string_view ToString(PromoterMotifMatch m) noexcept;

/// @brief Per-read classification result.
///
/// Emitted as auxiliary fields on the read; downstream code may forward to
/// BAM tags or Parquet columns. All string fields are owning copies so the
/// result outlives the catalog view.
struct ClassificationResult {
    /// Matched `haplotype_class` from the catalog, or empty if no entry was
    /// located. The serialized representation uses "unknown" for empty.
    std::string haplotype_class_call;

    /// Outcome of the canonical/duplicate motif test.
    PromoterMotifMatch promoter_motif_match{PromoterMotifMatch::None};

    /// Total k-mer mismatches across flanking + CH1 (the scan window).
    /// `0` means an exact match against the catalog reference; higher values
    /// indicate divergence or sequencing error.
    std::uint32_t flanking_kmer_mismatches{0};

    /// True iff a promoter-signature block was available and used.
    /// Distinguishes "motif test ran and said None" from "no motifs to test".
    bool promoter_test_ran{false};

    /// True iff `include_flanking_bp > 0` and the read could be extended into
    /// the flanking window (i.e. enough upstream sequence was present).
    bool flank_window_used{false};

    /// Serialize `haplotype_class_call` with "unknown" fallback.
    [[nodiscard]] std::string CallOrUnknown() const {
        return haplotype_class_call.empty() ? "unknown"
                                            : haplotype_class_call;
    }
};

/// @brief Input read description.
///
/// `read_seq` is the forward-strand read sequence (uppercase, A/C/G/T/N).
/// `mapping_start` is the 0-based start of the *primary* mapping on the
/// reference; the classifier uses this only when no haplotype-class
/// hypothesis is supplied (the catalog is queried by position).
///
/// `upstream_offset_in_read` is the position *within* `read_seq` (0-based)
/// that corresponds to the catalog's `anchor` (typically CH1_start). The
/// caller is responsible for computing this from the alignment / CIGAR; the
/// classifier treats it as a contract. If unset, the read start is used and
/// `flank_window_used` is set to false.
struct ReadInput {
    std::string_view read_seq;
    std::string_view target_id;
    std::uint64_t mapping_start{0};
    std::optional<std::uint32_t> upstream_offset_in_read;
    bool is_reverse{false};
};

/// @brief Optional caller-supplied classification context.
///
/// When the upstream caller already knows the candidate haplotype class
/// (e.g. PSV pipeline pre-classification), supplying `hypothesis` skips the
/// position-based catalog lookup. Useful for the test path and for cases
/// where the read maps inside the locus but to a position that is *itself*
/// ambiguous between canonical/chimdup.
struct ClassifyOptions {
    std::optional<std::string> hypothesis;  // haplotype_class tag
};

/// @brief Classify a single read against an IGH catalog entry.
///
/// @param catalog Catalog interface; the classifier looks up the entry by
///                the hypothesis (if supplied) or by mapping position.
/// @param read    Read sequence + mapping metadata.
/// @param opts    Optional hypothesis override.
/// @return        Per-read classification result. If no catalog entry was
///                found, `haplotype_class_call` is empty and
///                `promoter_test_ran == false`.
[[nodiscard]] ClassificationResult Classify(const ISegDupCatalog& catalog,
                                            const ReadInput& read,
                                            const ClassifyOptions& opts = {});

// ---------------------------------------------------------------------------
// Lower-level utilities — exposed for unit testing.

/// @brief Parse a window spec of the form "-100..0" into `(start, end)` with
/// negative-relative-to-anchor convention. Returns `nullopt` on parse fail.
[[nodiscard]] std::optional<std::pair<int, int>> ParseRelativeWindow(
    std::string_view spec);

/// @brief Test whether a catalog motif matches inside a window.
///
/// The motif may contain a literal `...` (three ASCII dots) splitting it
/// into prefix and suffix. Match semantics:
///   - No `...`: simple substring search within `window`.
///   - With `...`: prefix and suffix must appear in 5'->3' order.
[[nodiscard]] bool MotifMatchesWindow(std::string_view motif,
                                      std::string_view window);

/// @brief Hamming-style k-mer mismatch count of `read_window` vs
/// `reference_window`, sliding by 1 bp. Returns the *minimum* per-k-mer
/// mismatch count summed over the aligned overlap. Length mismatch is
/// tolerated: the shorter side defines the overlap.
[[nodiscard]] std::uint32_t FlankingKmerMismatches(
    std::string_view read_window,
    std::string_view reference_window,
    std::uint8_t kmer_size,
    std::uint8_t max_mismatch_per_kmer);

}  // namespace llmap::classify::igh
