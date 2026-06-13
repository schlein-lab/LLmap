// LLmap — Parent / processed-pseudogene pair catalog implementation.

#include "provenance/pseudogene_catalog.h"

#include <fstream>
#include <sstream>
#include <string>

namespace llmap::provenance {

namespace {

PseudogeneLocus Make(std::string chrom, std::uint64_t s, std::uint64_t e,
                     std::string parent, std::string pseudo, PseudogeneRole role) {
    PseudogeneLocus l;
    l.chrom = std::move(chrom);
    l.start = s;
    l.end = e;
    l.parent = std::move(parent);
    l.pseudogene = std::move(pseudo);
    l.role = role;
    l.verified = false;  // starter coordinates are approximate GRCh38
    return l;
}

}  // namespace

void PseudogeneCatalog::LoadBuiltinStarter() {
    // Approximate GRCh38 loci for the canonical, clinically-relevant pairs.
    // NOT precise gene models — they bound the locus for the per-read flag and
    // are flagged verified=false. Replace with LoadBed() for production.
    loci_ = {
        // GBA1 (parent, intron-bearing) ↔ GBAP1 (processed pseudogene), chr1q22.
        Make("chr1", 155'234'000, 155'245'000, "GBA1", "GBAP1",
             PseudogeneRole::Parent),
        Make("chr1", 155'213'000, 155'228'000, "GBA1", "GBAP1",
             PseudogeneRole::Pseudogene),
        // PMS2 (parent) ↔ PMS2CL (pseudogene), chr7p22.1.
        Make("chr7", 5'970'000, 5'991'000, "PMS2", "PMS2CL",
             PseudogeneRole::Parent),
        Make("chr7", 6'735'000, 6'748'000, "PMS2", "PMS2CL",
             PseudogeneRole::Pseudogene),
        // SMN1 ↔ SMN2 (near-identical paralogs / pseudogene-like), chr5q13.2.
        Make("chr5", 70'924'000, 70'953'000, "SMN1", "SMN2",
             PseudogeneRole::Parent),
        Make("chr5", 69'345'000, 69'374'000, "SMN1", "SMN2",
             PseudogeneRole::Pseudogene),
    };
}

bool PseudogeneCatalog::LoadBed(const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::vector<PseudogeneLocus> loaded;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        PseudogeneLocus l;
        std::string role;
        if (!(ss >> l.chrom >> l.start >> l.end >> l.parent >> l.pseudogene >>
              role)) {
            continue;  // skip malformed lines
        }
        l.role = (role == "pseudogene") ? PseudogeneRole::Pseudogene
                                        : PseudogeneRole::Parent;
        l.verified = true;  // a curated catalog is taken as authoritative
        loaded.push_back(std::move(l));
    }
    if (loaded.empty()) return false;
    loci_ = std::move(loaded);
    return true;
}

const PseudogeneLocus* PseudogeneCatalog::Lookup(std::string_view chrom,
                                                 std::uint64_t pos) const {
    for (const auto& l : loci_) {
        if (l.chrom == chrom && pos >= l.start && pos < l.end) return &l;
    }
    return nullptr;
}

}  // namespace llmap::provenance
