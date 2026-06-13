// LLmap — Parent / processed-pseudogene pair catalog.
//
// Feeds the `at_pseudogene_parent_locus` evidence of mapping_confusion: a read
// landing in a locus that belongs to a known parent-gene / processed-pseudogene
// pair, when it is INTRONLESS (no N ops) while the parent is intron-bearing, is
// pseudogene-derived — the classic confusion (GBAP1↔GBA1, PMS2↔PMS2CL,
// SMN1↔SMN2) that fakes intronless "expression" / mis-mapped variants.
//
// Two sources:
//   * a small built-in STARTER set (the canonical clinically-relevant pairs),
//     coordinates are approximate GRCh38 and flagged for verification;
//   * a loadable BED/TSV for a production-grade catalog (pseudogene.org /
//     GENCODE pseudogene biotype), same pattern as LLmap's other catalogs.
//
// Dependency-light (stdlib only); the integration layer queries Lookup() per
// read to set the evidence flag.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::provenance {

enum class PseudogeneRole : std::uint8_t { Parent, Pseudogene };

struct PseudogeneLocus {
    std::string   chrom;       // e.g. "chr1" (matches the reference naming)
    std::uint64_t start{0};    // 0-based half-open, GRCh38
    std::uint64_t end{0};
    std::string   parent;      // parent gene symbol (the intron-bearing copy)
    std::string   pseudogene;  // processed-pseudogene symbol
    PseudogeneRole role{PseudogeneRole::Parent};
    bool          verified{false};  // false ⇒ approximate starter coordinate
};

class PseudogeneCatalog {
public:
    // Built-in clinically-relevant starter pairs (approximate GRCh38; verified
    // flag false). For production, prefer LoadBed() with a curated source.
    void LoadBuiltinStarter();

    // Load from a BED-like TSV: chrom, start, end, parent, pseudogene, role
    // ("parent"|"pseudogene"). Returns false on parse failure.
    [[nodiscard]] bool LoadBed(const std::string& path);

    // The locus covering (chrom, pos), or nullopt. pos is 0-based.
    [[nodiscard]] const PseudogeneLocus* Lookup(std::string_view chrom,
                                                std::uint64_t pos) const;

    [[nodiscard]] std::size_t Size() const noexcept { return loci_.size(); }
    [[nodiscard]] const std::vector<PseudogeneLocus>& Entries() const noexcept {
        return loci_;
    }

private:
    std::vector<PseudogeneLocus> loci_;
};

}  // namespace llmap::provenance
