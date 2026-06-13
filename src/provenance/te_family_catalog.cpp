// LLmap — TE family catalog implementation.

#include "provenance/te_family_catalog.h"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace llmap::provenance {

const char* TeClassTag(TeClass c) noexcept {
    switch (c) {
        case TeClass::Alu:         return "alu";
        case TeClass::Line1:       return "l1";
        case TeClass::Sva:         return "sva";
        case TeClass::Herv:        return "herv";
        case TeClass::OtherRepeat: return "rep";
    }
    return "rep";
}

namespace {
TeClass ParseClass(std::string_view s) {
    if (s == "alu" || s == "Alu" || s == "SINE/Alu") return TeClass::Alu;
    if (s == "l1" || s == "L1" || s == "LINE/L1")    return TeClass::Line1;
    if (s == "sva" || s == "SVA")                    return TeClass::Sva;
    if (s == "herv" || s == "HERV" || s.rfind("LTR", 0) == 0) return TeClass::Herv;
    return TeClass::OtherRepeat;
}
void SortLoci(std::vector<TeLocus>& v) {
    std::sort(v.begin(), v.end(), [](const TeLocus& a, const TeLocus& b) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.start < b.start;
    });
}
}  // namespace

void TeFamilyCatalog::LoadBuiltinStarter() {
    // A handful of representative GRCh38 TE loci (approximate; for tests/smoke —
    // production uses LoadRepeatMasker/LoadBed). Young (low-divergence) elements
    // are the ambiguous ones (recent insertions, near-identical to consensus).
    loci_ = {
        {"chr1", 145000000, 145000300, TeClass::Alu,   "AluYa5", 0.5f, '+'},
        {"chr1", 145100000, 145106000, TeClass::Line1, "L1HS",   1.0f, '-'},
        {"chr1", 145200000, 145201800, TeClass::Sva,   "SVA_E",  2.0f, '+'},
        {"chr7", 142000000, 142000310, TeClass::Alu,   "AluY",   1.5f, '+'},
        {"chr19", 40000000, 40006000,  TeClass::Line1, "L1PA2",  6.0f, '+'},
    };
    SortLoci(loci_);
}

bool TeFamilyCatalog::LoadBed(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    loci_.clear();
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::vector<std::string> t;
        std::string c;
        while (std::getline(ss, c, '\t')) t.push_back(c);
        if (t.size() < 3) continue;
        TeLocus l;
        l.chrom = t[0];
        l.start = std::strtoull(t[1].c_str(), nullptr, 10);
        l.end   = std::strtoull(t[2].c_str(), nullptr, 10);
        if (t.size() > 3) l.cls = ParseClass(t[3]);
        if (t.size() > 4) l.subfamily = t[4];
        if (t.size() > 5) l.divergence = std::strtof(t[5].c_str(), nullptr);
        if (t.size() > 6 && !t[6].empty()) l.strand = t[6][0];
        loci_.push_back(std::move(l));
    }
    SortLoci(loci_);
    return true;
}

const TeLocus* TeFamilyCatalog::Lookup(std::string_view chrom,
                                       std::uint64_t pos) const {
    // loci_ sorted by (chrom, start). upper_bound → first locus AFTER (chrom,pos)
    // [i.e. chrom greater, or same chrom with start > pos]; the predecessor is the
    // last locus on this chrom with start <= pos. Non-overlapping starter/RM data
    // ⇒ that single predecessor is the only containment candidate.
    auto it = std::upper_bound(
        loci_.begin(), loci_.end(), pos,
        [chrom](std::uint64_t p, const TeLocus& l) {
            if (std::string_view(l.chrom) != chrom)
                return std::string_view(l.chrom) > chrom;
            return l.start > p;
        });
    if (it == loci_.begin()) return nullptr;
    --it;
    if (std::string_view(it->chrom) == chrom && it->start <= pos && pos < it->end)
        return &*it;
    return nullptr;
}

}  // namespace llmap::provenance
