// LLmap — Adapter implementation (catalog::SegDupCatalog → ISegDupCatalog).

#include "classify/segdup_catalog_adapter.h"

#include "catalog/segdup_catalog.h"

#include <utility>

namespace llmap::classify::igh {

namespace {

CatalogEntryView make_view(const catalog::SegDupCatalogEntry& e) {
    CatalogEntryView v;
    v.locus_id        = e.locus_id;
    v.haplotype_class = e.haplotype_class;

    v.mapping_primary.kmer_size =
        static_cast<std::uint8_t>(std::max(0, e.mapping_primary.kmer_size));
    v.mapping_primary.max_mismatch =
        static_cast<std::uint8_t>(std::max(0, e.mapping_primary.max_mismatch));
    v.mapping_primary.include_flanking_bp =
        static_cast<std::uint32_t>(std::max(0, e.mapping_primary.include_flanking_bp));
    v.mapping_primary.include_flanking_anchor =
        e.mapping_primary.include_flanking_anchor;

    if (e.promoter_signature.has_value()) {
        PromoterSignature p;
        p.window_relative_to_anchor =
            e.promoter_signature->window_relative_to_anchor;
        p.anchor = e.promoter_signature->anchor;
        p.canonical_motif = e.promoter_signature->canonical_motif;
        p.duplicate_motif = e.promoter_signature->duplicate_motif;
        v.promoter_signature = std::move(p);
    }
    return v;
}

}  // namespace

std::optional<CatalogEntryView>
SegDupCatalogAdapter::LookupByPosition(std::string_view target_id,
                                        std::uint64_t pos) const {
    // `lookup_by_coords` is the most-specific lookup (nested-aware).
    auto hit = cat_.lookup_by_coords(assembly_, target_id,
                                      static_cast<std::int64_t>(pos));
    if (!hit) return std::nullopt;
    return make_view(*hit);
}

std::optional<CatalogEntryView>
SegDupCatalogAdapter::LookupByHaplotypeClass(
    std::string_view haplotype_class) const {
    // The catalog has no haplotype-keyed index — we scan every curated
    // record linearly via for_each_curated(). The curated tier is expected
    // to stay well under a thousand entries (currently 24, target ~300 by
    // v2026.Q4), so a linear scan is cheap; this method is also never
    // called from the alignment hot loop.
    //
    // First match wins. If two curated entries happen to share the same
    // haplotype_class string the catalog is malformed — we surface the
    // first one as a best-effort hit rather than failing, but a future
    // catalog validator should reject duplicates upstream.
    std::optional<CatalogEntryView> found;
    cat_.for_each_curated(
        [&](const catalog::SegDupCatalogEntry& e) {
            if (found) return;                                  // first wins
            if (e.haplotype_class == haplotype_class) {
                found = make_view(e);
            }
        });
    return found;
}

}  // namespace llmap::classify::igh
