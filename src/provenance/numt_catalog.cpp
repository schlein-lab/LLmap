// LLmap — NUMT catalog implementation.

#include "provenance/numt_catalog.h"

#include <fstream>
#include <sstream>

namespace llmap::provenance {

namespace {

NumtLocus Make(std::string chrom, std::uint64_t s, std::uint64_t e,
               std::string mt_region) {
    NumtLocus l;
    l.chrom = std::move(chrom);
    l.start = s;
    l.end = e;
    l.mt_region = std::move(mt_region);
    l.verified = false;  // approximate starter coordinate
    return l;
}

}  // namespace

void NumtCatalog::LoadBuiltinStarter() {
    // A few well-documented large NUMTs (approximate GRCh38; verified=false).
    // NOT exhaustive — replace with LoadBed() of the UCSC NumtS track / dinumt
    // calls for production. They bound the locus for the at_mt_homologous flag.
    loci_ = {
        Make("chr1",  564'000,   570'000,   "MT:whole"),   // chr1 large NUMT cluster
        Make("chr5",  79'945'000, 79'951'000, "MT:ND2-COX1"),
        Make("chr17", 22'020'000, 22'027'000, "MT:whole"),  // chr17 large NUMT
        Make("chrX",  153'000'000, 153'003'000, "MT:CYTB"),
    };
}

bool NumtCatalog::LoadBed(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::vector<NumtLocus> loaded;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        NumtLocus l;
        if (!(ss >> l.chrom >> l.start >> l.end)) continue;  // chrom start end
        ss >> l.mt_region;     // optional 4th column
        l.verified = true;     // a curated catalog is authoritative
        loaded.push_back(std::move(l));
    }
    if (loaded.empty()) return false;
    loci_ = std::move(loaded);
    return true;
}

const NumtLocus* NumtCatalog::Lookup(std::string_view chrom,
                                     std::uint64_t pos) const {
    for (const auto& l : loci_) {
        if (l.chrom == chrom && pos >= l.start && pos < l.end) return &l;
    }
    return nullptr;
}

}  // namespace llmap::provenance
