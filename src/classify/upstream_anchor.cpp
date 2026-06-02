// LLmap — IGH-mode upstream-anchor classifier (implementation).

#include "classify/upstream_anchor.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <string>
#include <string_view>

namespace llmap::classify::igh {

namespace {

constexpr std::string_view kEllipsis{"..."};

/// Uppercase a single byte without locale.
constexpr char Upper(char c) noexcept {
    return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
}

/// Case-insensitive substring search. Returns position or
/// `std::string_view::npos`.
[[nodiscard]] std::size_t FindCI(std::string_view hay,
                                 std::string_view needle,
                                 std::size_t from = 0) noexcept {
    if (needle.empty() || hay.size() < needle.size() + from) {
        return needle.empty() ? from : std::string_view::npos;
    }
    const std::size_t last = hay.size() - needle.size();
    for (std::size_t i = from; i <= last; ++i) {
        bool ok = true;
        for (std::size_t j = 0; j < needle.size(); ++j) {
            if (Upper(hay[i + j]) != Upper(needle[j])) {
                ok = false;
                break;
            }
        }
        if (ok) return i;
    }
    return std::string_view::npos;
}

/// Locate the upstream window inside the read.
///
/// Given the catalog window "-100..0" and an anchor offset within the read,
/// returns the [begin, end) slice covering bases relative to the anchor.
/// The window is *clipped* to read bounds; the caller learns whether
/// clipping happened via the optional `clipped` out-param.
[[nodiscard]] std::string_view SliceUpstream(
    std::string_view read,
    std::uint32_t anchor_offset_in_read,
    int rel_begin,
    int rel_end,
    bool* clipped) {
    long long lo = static_cast<long long>(anchor_offset_in_read) + rel_begin;
    long long hi = static_cast<long long>(anchor_offset_in_read) + rel_end;
    bool was_clipped = false;
    if (lo < 0) {
        lo = 0;
        was_clipped = true;
    }
    if (hi > static_cast<long long>(read.size())) {
        hi = static_cast<long long>(read.size());
        was_clipped = true;
    }
    if (hi < lo) hi = lo;
    if (clipped) *clipped = was_clipped;
    return read.substr(static_cast<std::size_t>(lo),
                       static_cast<std::size_t>(hi - lo));
}

/// Determine whether `s` looks like a signed integer ("-100", "0", "200").
[[nodiscard]] bool ParseSignedInt(std::string_view s, int& out) noexcept {
    if (s.empty()) return false;
    auto first = s.data();
    auto last = s.data() + s.size();
    int value{};
    auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc{} || ptr != last) return false;
    out = value;
    return true;
}

}  // namespace

std::string_view ToString(PromoterMotifMatch m) noexcept {
    switch (m) {
        case PromoterMotifMatch::None:      return "none";
        case PromoterMotifMatch::Canonical: return "canonical";
        case PromoterMotifMatch::Chimdup:   return "chimdup";
        case PromoterMotifMatch::Mixed:     return "mixed";
    }
    return "none";
}

std::optional<std::pair<int, int>> ParseRelativeWindow(std::string_view spec) {
    // Format: "<int>..<int>", e.g. "-100..0".
    auto sep = spec.find("..");
    if (sep == std::string_view::npos) return std::nullopt;
    int lo{};
    int hi{};
    if (!ParseSignedInt(spec.substr(0, sep), lo)) return std::nullopt;
    if (!ParseSignedInt(spec.substr(sep + 2), hi)) return std::nullopt;
    if (hi < lo) std::swap(lo, hi);
    return std::make_pair(lo, hi);
}

bool MotifMatchesWindow(std::string_view motif, std::string_view window) {
    if (motif.empty() || window.empty()) return false;
    auto sep = motif.find(kEllipsis);
    if (sep == std::string_view::npos) {
        return FindCI(window, motif) != std::string_view::npos;
    }
    std::string_view prefix = motif.substr(0, sep);
    std::string_view suffix = motif.substr(sep + kEllipsis.size());
    if (prefix.empty() && suffix.empty()) return false;
    std::size_t prefix_end = 0;
    if (!prefix.empty()) {
        auto p = FindCI(window, prefix);
        if (p == std::string_view::npos) return false;
        prefix_end = p + prefix.size();
    }
    if (suffix.empty()) return true;
    return FindCI(window, suffix, prefix_end) != std::string_view::npos;
}

std::uint32_t FlankingKmerMismatches(std::string_view read_window,
                                     std::string_view reference_window,
                                     std::uint8_t kmer_size,
                                     std::uint8_t max_mismatch_per_kmer) {
    if (kmer_size == 0) return 0;
    const std::size_t k = static_cast<std::size_t>(kmer_size);
    const std::size_t overlap = std::min(read_window.size(),
                                         reference_window.size());
    if (overlap < k) {
        // Too short to form any k-mer — fall back to plain Hamming on overlap.
        std::uint32_t mm = 0;
        for (std::size_t i = 0; i < overlap; ++i) {
            if (Upper(read_window[i]) != Upper(reference_window[i])) ++mm;
        }
        return mm;
    }
    std::uint32_t total = 0;
    for (std::size_t i = 0; i + k <= overlap; ++i) {
        std::uint32_t mm = 0;
        for (std::size_t j = 0; j < k; ++j) {
            if (Upper(read_window[i + j]) != Upper(reference_window[i + j])) {
                ++mm;
            }
        }
        // Clamp per-k-mer contribution by `max_mismatch_per_kmer`. Catalog
        // semantics: chains tolerate up to max_mismatch per k-mer; k-mers
        // exceeding that contribute the cap, not the raw count (avoids
        // hot-spot domination of the global score).
        if (max_mismatch_per_kmer > 0 && mm > max_mismatch_per_kmer) {
            mm = max_mismatch_per_kmer;
        }
        total += mm;
    }
    return total;
}

namespace {

/// Resolve the catalog entry to use for a read.
[[nodiscard]] std::optional<CatalogEntryView> ResolveEntry(
    const ISegDupCatalog& catalog,
    const ReadInput& read,
    const ClassifyOptions& opts) {
    if (opts.hypothesis.has_value() && !opts.hypothesis->empty()) {
        auto v = catalog.LookupByHaplotypeClass(*opts.hypothesis);
        if (v.has_value()) return v;
    }
    return catalog.LookupByPosition(read.target_id, read.mapping_start);
}

/// Per-motif test using the catalog's promoter signature.
[[nodiscard]] PromoterMotifMatch TestPromoter(
    const PromoterSignature& sig,
    std::string_view upstream_window) {
    const bool has_can = sig.canonical_motif.has_value()
                         && !sig.canonical_motif->empty();
    const bool has_dup = sig.duplicate_motif.has_value()
                         && !sig.duplicate_motif->empty();
    if (!has_can && !has_dup) return PromoterMotifMatch::None;
    const bool can_hit = has_can &&
                         MotifMatchesWindow(*sig.canonical_motif,
                                            upstream_window);
    const bool dup_hit = has_dup &&
                         MotifMatchesWindow(*sig.duplicate_motif,
                                            upstream_window);
    if (can_hit && dup_hit) return PromoterMotifMatch::Mixed;
    if (can_hit)            return PromoterMotifMatch::Canonical;
    if (dup_hit)            return PromoterMotifMatch::Chimdup;
    return PromoterMotifMatch::None;
}

}  // namespace

ClassificationResult Classify(const ISegDupCatalog& catalog,
                              const ReadInput& read,
                              const ClassifyOptions& opts) {
    ClassificationResult out;

    auto entry = ResolveEntry(catalog, read, opts);
    if (!entry.has_value()) {
        // No catalog coverage — leave defaults (everything blank/false).
        return out;
    }
    out.haplotype_class_call = entry->haplotype_class;

    // Locate the anchor inside the read. If the caller did not supply an
    // offset, the read start is treated as the anchor (used by the
    // CH1-only test case).
    const std::uint32_t anchor_offset =
        read.upstream_offset_in_read.value_or(0);

    // Build the upstream window. Catalog supplies the relative span; if it
    // is missing or unparsable, fall back to `[-include_flanking_bp, 0]`
    // so the flank window is still honoured.
    int rel_lo = -static_cast<int>(entry->mapping_primary.include_flanking_bp);
    int rel_hi = 0;
    if (entry->promoter_signature.has_value()) {
        auto parsed = ParseRelativeWindow(
            entry->promoter_signature->window_relative_to_anchor);
        if (parsed.has_value()) {
            rel_lo = parsed->first;
            rel_hi = parsed->second;
        }
    }
    bool clipped = false;
    std::string_view upstream =
        SliceUpstream(read.read_seq, anchor_offset, rel_lo, rel_hi, &clipped);

    // Flank window is "used" if (a) the catalog asked for flanking and
    // (b) we actually got >0 bp of upstream sequence to scan (i.e. the
    // anchor was not at the very read start).
    out.flank_window_used =
        entry->mapping_primary.include_flanking_bp > 0
        && !upstream.empty();

    // Promoter motif test runs whenever the catalog entry carries the
    // signature; it does *not* require flanking to be reserved (a read may
    // overlap the upstream window even at include_flanking_bp == 0, if the
    // alignment naturally extends into 5' UTR).
    if (entry->promoter_signature.has_value()) {
        out.promoter_test_ran = true;
        out.promoter_motif_match =
            TestPromoter(*entry->promoter_signature, upstream);
    }

    // K-mer mismatch scoring. For the binary discriminant test, the
    // reference window is the *canonical* motif (without the ellipsis).
    // We score the read window against the canonical motif's literal
    // content, so canonical reads return 0 and chimdup reads return a
    // positive count. This mirrors the alignment-free principle used
    // throughout LLmap (see `igh_anchor_catalog.h`).
    if (entry->promoter_signature.has_value()
        && entry->promoter_signature->canonical_motif.has_value()) {
        std::string ref_literal;
        const std::string& motif = *entry->promoter_signature->canonical_motif;
        for (char c : motif) {
            if (c == '.') continue;  // drop ellipsis bytes
            ref_literal.push_back(c);
        }
        out.flanking_kmer_mismatches = FlankingKmerMismatches(
            upstream, ref_literal,
            entry->mapping_primary.kmer_size,
            entry->mapping_primary.max_mismatch);
    }

    return out;
}

}  // namespace llmap::classify::igh
