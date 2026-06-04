// LLmap — junction_hunter: variable-tier cascade caller impl.

#include "junction_hunter/cascade_caller.h"
#include "junction_hunter/pair_kmer_index.h"
#include "junction_hunter/read_tiler.h"
#include "junction_hunter/consensus_caller.h"

#include <algorithm>
#include <array>
#include <cstdint>

namespace llmap::junction_hunter {

CascadeConfig CascadeConfig::LongReadPreset() {
    CascadeConfig c;
    // Default cascade: small-k membership + PSV-grade discrimination.
    // Same range works for short reads (k_max=125 fits in 150 bp);
    // long reads may extend beyond 125 via the --k-values override.
    c.k_values = {11, 17, 25, 35, 51, 71, 101, 125};
    c.tier1_min_anchor_hits = 5;
    c.consensus_min = 3;
    c.monotonicity_min = 0.95f;
    return c;
}

CascadeConfig CascadeConfig::ShortReadPreset() {
    // Same as long-read preset — the cascade covers k=11..125 which is
    // the exact range a 150 bp short read needs when a NAHR junction
    // splits the read into a small fragment (e.g. 25 bp) from one
    // copy and a larger fragment (125 bp) from the other.
    return LongReadPreset();
}

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

struct Tier1HitCounts {
    std::uint32_t up{0};
    std::uint32_t dn{0};
    std::uint32_t inn{0};
    std::uint32_t amb{0};
    std::uint32_t signal() const { return up + dn + inn; }
};

Tier1HitCounts CountTier1Hits(std::string_view read_seq,
                              std::uint8_t k,
                              const KmerClassMap& tbl) {
    Tier1HitCounts c;
    if (read_seq.size() < k) return c;
    for (std::size_t i = 0; i + k <= read_seq.size(); ++i) {
        bool clean = true;
        for (std::size_t j = 0; j < k; ++j) {
            if (!IsAcgt(read_seq[i + j])) { clean = false; break; }
        }
        if (!clean) continue;
        auto h = HashKmer(read_seq.substr(i, k));
        auto it = tbl.find(h);
        if (it == tbl.end()) continue;
        switch (it->second.cls) {
            case LocusClass::LcrUp:     ++c.up;  break;
            case LocusClass::LcrDown:   ++c.dn;  break;
            case LocusClass::Interior:  ++c.inn; break;
            case LocusClass::Ambiguous: ++c.amb; break;
            default: break;
        }
    }
    return c;
}

/// Build a multi-k ReadTiling matching the cascade's k_values. Only
/// emits tilings at k for which the pair index has a built tier.
ReadTiling TileForCascade(std::string_view read_seq,
                          const CascadePairIndex& pair_index) {
    ReadTiling out;
    // ReadTiling fixed-size 5 array — we fit by taking the first 5
    // built tiers; the consensus caller votes over what we provide.
    // (Future refactor: ReadTiling itself to vector<vector<...>>.)
    std::array<std::uint8_t, 5> emit_k{};
    std::size_t n_emit = 0;
    for (std::size_t ti = 0;
         ti < pair_index.k_values.size() && n_emit < emit_k.size(); ++ti) {
        if (pair_index.tier_built[ti]) emit_k[n_emit++] = pair_index.k_values[ti];
    }
    while (n_emit < emit_k.size()) emit_k[n_emit++] = emit_k[0];

    out.k_values = emit_k;
    for (std::size_t ki = 0; ki < emit_k.size(); ++ki) {
        const std::uint8_t k = emit_k[ki];
        if (read_seq.size() < k) continue;
        for (std::size_t i = 0; i + k <= read_seq.size(); ++i) {
            bool clean = true;
            for (std::size_t j = 0; j < k; ++j) {
                if (!IsAcgt(read_seq[i + j])) { clean = false; break; }
            }
            if (!clean) continue;
            out.per_k_hashes[ki].push_back(
                {HashKmer(read_seq.substr(i, k)), static_cast<std::uint32_t>(i)});
        }
    }
    return out;
}

/// Adapt the variable-tier cascade pair-index to the flat fixed-5-tier
/// PairKmerIndex expected by the existing CallJunction. Copies up to
/// 5 built tiers in cascade-order; un-built tiers contribute empty maps.
PairKmerIndex AdaptToFlat(const CascadePairIndex& cpi,
                           const std::array<std::uint8_t, 5>& emit_k) {
    PairKmerIndex pki;
    pki.pair_id = cpi.pair_meta.pair_id;
    pki.k_values = emit_k;
    std::size_t emit_i = 0;
    for (std::size_t ti = 0;
         ti < cpi.k_values.size() && emit_i < pki.per_k_class.size(); ++ti) {
        if (cpi.tier_built[ti]) {
            pki.per_k_class[emit_i] = cpi.tier_class[ti];
            ++emit_i;
        }
    }
    return pki;
}

}  // namespace

JunctionRecord CascadeCall(std::string_view read_id,
                           std::string_view read_seq,
                           const CascadePairIndex& pair_index,
                           const CascadeConfig& cfg) {
    JunctionRecord rec;
    rec.read_id = read_id;
    rec.pair_id = pair_index.pair_meta.pair_id;

    if (cfg.k_values.empty() || pair_index.tier_class.empty()) {
        rec.call = JunctionCall::Unmapped;
        return rec;
    }

    // Tier 0 fast-path: count k-mers at smallest k against tier_class[0].
    auto t0 = CountTier1Hits(read_seq, cfg.k_values[0], pair_index.tier_class[0]);
    if (t0.signal() < cfg.tier1_min_anchor_hits) {
        rec.call = JunctionCall::Unmapped;
        rec.n_kmer_total = 0;
        return rec;
    }

    // Tier 1+ have been built externally before invocation. Tile + vote
    // via the existing consensus caller using whatever tiers exist.
    auto tiling = TileForCascade(read_seq, pair_index);
    auto flat = AdaptToFlat(pair_index, tiling.k_values);

    MultiKConfig mk;
    mk.k_values = tiling.k_values;
    mk.consensus_min = cfg.consensus_min;
    mk.monotonicity_min = cfg.monotonicity_min;
    mk.min_psv_switches = 3;

    return CallJunction(read_id, tiling, flat, pair_index.pair_meta, mk);
}

}  // namespace llmap::junction_hunter
