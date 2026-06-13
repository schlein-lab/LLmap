// LLmap — NUMT (nuclear-mitochondrial segment) catalog — Block 3.
//
// Feeds the `at_mt_homologous_locus` evidence of mapping_confusion: a read that
// lands at a nuclear locus homologous to mtDNA (a NUMT) is the classic source
// of fake low-frequency "heteroplasmy". The catalog marks those loci; the
// numt-vs-real-mt-heteroplasmy decision then uses identity_to_mt vs
// identity_to_nuclear_numt (computed in the pipeline) — Layer-1 `numt` artifact
// when nuclear wins, Layer-3 `bio:mthet` when the mt reference wins.
//
// Two sources (same pattern as pseudogene_catalog): a small built-in STARTER of
// well-known large NUMTs (approximate GRCh38, verified=false) + a loadable BED
// for a production catalog (UCSC NumtS track / Dayama et al. dinumt calls).
// Dependency-light (stdlib only); the integration layer queries Lookup() per read.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::provenance {

struct NumtLocus {
    std::string   chrom;           // nuclear chromosome, e.g. "chr1"
    std::uint64_t start{0};        // 0-based half-open, GRCh38
    std::uint64_t end{0};
    std::string   mt_region;       // homologous mt segment, e.g. "MT:1-1200" / "ND4"
    bool          verified{false}; // false ⇒ approximate starter coordinate
};

class NumtCatalog {
public:
    // Built-in well-known large NUMTs (approximate GRCh38; verified=false).
    // For production, prefer LoadBed() with the UCSC NumtS track.
    void LoadBuiltinStarter();

    // BED-like TSV: chrom, start, end, mt_region. Returns false on parse failure.
    [[nodiscard]] bool LoadBed(const std::string& path);

    // The NUMT locus covering (chrom, pos), or nullptr. pos is 0-based.
    [[nodiscard]] const NumtLocus* Lookup(std::string_view chrom,
                                          std::uint64_t pos) const;

    [[nodiscard]] std::size_t Size() const noexcept { return loci_.size(); }
    [[nodiscard]] const std::vector<NumtLocus>& Entries() const noexcept {
        return loci_;
    }

private:
    std::vector<NumtLocus> loci_;
};

}  // namespace llmap::provenance
