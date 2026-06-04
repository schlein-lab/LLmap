// LLmap — junction_hunter: per-pair multi-k k-mer index implementation.

#include "junction_hunter/pair_kmer_index.h"

#include <cstdint>

namespace llmap::junction_hunter {

namespace {

inline char Comp(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

inline bool IsAcgt(char c) noexcept {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T'
        || c == 'a' || c == 'c' || c == 'g' || c == 't';
}

std::uint64_t HashKmer(std::string_view k) noexcept {
    std::uint64_t h_fwd = 1469598103934665603ULL;
    std::uint64_t h_rev = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    for (std::size_t i = 0; i < k.size(); ++i) {
        char cf = k[i];
        char cr = Comp(k[k.size() - 1 - i]);
        h_fwd = (h_fwd ^ static_cast<std::uint64_t>(cf)) * fnv_prime;
        h_rev = (h_rev ^ static_cast<std::uint64_t>(cr)) * fnv_prime;
    }
    return h_fwd < h_rev ? h_fwd : h_rev;
}

void IndexSegment(std::string_view seq,
                  std::uint8_t k,
                  LocusClass cls,
                  KmerClassMap& tbl,
                  std::size_t& unique_count,
                  std::size_t& amb_count) {
    if (seq.size() < k) return;
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        bool clean = true;
        for (std::size_t j = 0; j < k; ++j) {
            if (!IsAcgt(seq[i + j])) { clean = false; break; }
        }
        if (!clean) continue;
        std::uint64_t h = HashKmer(seq.substr(i, k));
        auto [it, inserted] = tbl.try_emplace(
            h, KmerLoc{cls, static_cast<std::uint32_t>(i)});
        if (inserted) {
            ++unique_count;
        } else if (it->second.cls != cls
                   && it->second.cls != LocusClass::Ambiguous) {
            it->second.cls = LocusClass::Ambiguous;
            ++amb_count;
        }
    }
}

}  // namespace

const char* LocusClassName(LocusClass c) noexcept {
    switch (c) {
        case LocusClass::None:      return "none";
        case LocusClass::LcrUp:     return "lcr_up";
        case LocusClass::LcrDown:   return "lcr_down";
        case LocusClass::Interior:  return "interior";
        case LocusClass::Ambiguous: return "ambiguous";
    }
    return "unknown";
}

PairIndexBuildStats BuildPairKmerIndex(const NahrPair& pair,
                                       std::string_view lcr_up_seq,
                                       std::string_view lcr_down_seq,
                                       std::string_view interior_seq,
                                       const MultiKConfig& cfg,
                                       PairKmerIndex& out) {
    PairIndexBuildStats stats;
    out.pair_id = pair.pair_id;
    out.k_values = cfg.k_values;
    for (auto& m : out.per_k_class) m.clear();

    for (std::size_t ki = 0; ki < cfg.k_values.size(); ++ki) {
        const std::uint8_t k = cfg.k_values[ki];
        std::size_t uniq = 0, amb = 0;
        IndexSegment(lcr_up_seq,   k, LocusClass::LcrUp,    out.per_k_class[ki], uniq, amb);
        IndexSegment(lcr_down_seq, k, LocusClass::LcrDown,  out.per_k_class[ki], uniq, amb);
        IndexSegment(interior_seq, k, LocusClass::Interior, out.per_k_class[ki], uniq, amb);
        stats.unique_kmers_per_k[ki] = uniq;
        stats.ambiguous_per_k[ki] = amb;
        stats.total_unique += uniq;
    }
    return stats;
}

}  // namespace llmap::junction_hunter
