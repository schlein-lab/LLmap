// LLmap — junction_hunter: variable-tier cascade pair-index impl.

#include "junction_hunter/cascade_pair_index.h"

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
    std::uint64_t hf = 1469598103934665603ULL;
    std::uint64_t hr = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    for (std::size_t i = 0; i < k.size(); ++i) {
        char cf = k[i];
        char cr = Comp(k[k.size() - 1 - i]);
        hf = (hf ^ static_cast<std::uint64_t>(cf)) * fnv_prime;
        hr = (hr ^ static_cast<std::uint64_t>(cr)) * fnv_prime;
    }
    return hf < hr ? hf : hr;
}

void IndexSegment(std::string_view seq,
                  std::uint8_t k,
                  LocusClass cls,
                  KmerClassMap& tbl,
                  std::size_t& unique_count) {
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
        }
    }
}

}  // namespace

std::size_t BuildCascadeTier1(const NahrPair& pair,
                              std::string_view lcr_up_seq,
                              std::string_view lcr_down_seq,
                              std::string_view interior_seq,
                              CascadePairIndex& out) {
    out.pair_meta = pair;
    if (out.tier_class.empty()) out.tier_class.resize(out.k_values.size());
    if (out.tier_built.empty()) out.tier_built.assign(out.k_values.size(), false);
    if (out.k_values.empty()) return 0;
    out.tier_class[0].clear();
    std::size_t uniq = 0;
    const std::uint8_t k = out.k_values[0];
    IndexSegment(lcr_up_seq,   k, LocusClass::LcrUp,    out.tier_class[0], uniq);
    IndexSegment(lcr_down_seq, k, LocusClass::LcrDown,  out.tier_class[0], uniq);
    IndexSegment(interior_seq, k, LocusClass::Interior, out.tier_class[0], uniq);
    out.tier_built[0] = true;
    return uniq;
}

std::size_t BuildCascadeTierN(std::size_t tier_idx,
                              std::string_view lcr_up_seq,
                              std::string_view lcr_down_seq,
                              std::string_view interior_seq,
                              CascadePairIndex& idx) {
    if (tier_idx >= idx.k_values.size()) return 0;
    if (idx.tier_built[tier_idx]) return 0;
    auto& tbl = idx.tier_class[tier_idx];
    tbl.clear();
    std::size_t uniq = 0;
    const std::uint8_t k = idx.k_values[tier_idx];
    IndexSegment(lcr_up_seq,   k, LocusClass::LcrUp,    tbl, uniq);
    IndexSegment(lcr_down_seq, k, LocusClass::LcrDown,  tbl, uniq);
    IndexSegment(interior_seq, k, LocusClass::Interior, tbl, uniq);
    idx.tier_built[tier_idx] = true;
    return uniq;
}

}  // namespace llmap::junction_hunter
