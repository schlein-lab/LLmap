// LLmap — junction_hunter: Stage-1 anchor-table builder implementation.

#include "junction_hunter/anchor_table_builder.h"

#include <algorithm>
#include <cstdint>
#include <string_view>

namespace llmap::junction_hunter {

namespace {

/// xxhash-style 64-bit canonical k-mer hash. We use the SAME primitive
/// as the rest of LLmap (MinimizerIndex) so that hashes computed here
/// agree bit-for-bit with hashes computed at lookup time.
///
/// Implementation: forward and reverse-complement mixed via a stable
/// FNV/Wang-style mix. Replace with the central llmap::core::HashKmer
/// once the module dependency is wired.
std::uint64_t HashKmer(std::string_view k) noexcept {
    std::uint64_t h_fwd = 1469598103934665603ULL;
    std::uint64_t h_rev = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    for (std::size_t i = 0; i < k.size(); ++i) {
        char c_fwd = k[i];
        char c_rev = [&]{
            char r = k[k.size() - 1 - i];
            switch (r) {
                case 'A': case 'a': return 'T';
                case 'C': case 'c': return 'G';
                case 'G': case 'g': return 'C';
                case 'T': case 't': return 'A';
                default: return 'N';
            }
        }();
        h_fwd = (h_fwd ^ static_cast<std::uint64_t>(c_fwd)) * fnv_prime;
        h_rev = (h_rev ^ static_cast<std::uint64_t>(c_rev)) * fnv_prime;
    }
    return h_fwd < h_rev ? h_fwd : h_rev;
}

bool IsClean(std::string_view k) noexcept {
    for (char c : k) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') return false;
    }
    return true;
}

/// Walk `seq` at stride 1 and collect candidate (hash, hit-count) for
/// every clean k-mer that passes the hit-count cap. Returns at most
/// `wanted` hashes, evenly spaced by position to avoid clustering.
void CollectCandidates(std::string_view seq,
                       std::uint8_t k,
                       std::uint32_t max_hits,
                       const GenomeHitCounter& hits_of,
                       std::uint32_t wanted,
                       std::vector<std::uint64_t>& out) {
    if (seq.size() < k) return;
    std::vector<std::pair<std::size_t, std::uint64_t>> kept;
    kept.reserve(seq.size() / 4);
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        auto win = seq.substr(i, k);
        if (!IsClean(win)) continue;
        std::uint64_t h = HashKmer(win);
        if (hits_of(h) > max_hits) continue;
        kept.emplace_back(i, h);
    }
    if (kept.empty()) return;
    // Evenly sample `wanted` positions across the kept vector.
    if (kept.size() <= wanted) {
        for (auto& p : kept) out.push_back(p.second);
        return;
    }
    double step = static_cast<double>(kept.size()) / static_cast<double>(wanted);
    for (std::uint32_t s = 0; s < wanted; ++s) {
        auto idx = static_cast<std::size_t>(s * step);
        if (idx >= kept.size()) idx = kept.size() - 1;
        out.push_back(kept[idx].second);
    }
}

}  // namespace

PairAnchorTable BuildPairAnchorTable(const NahrPair& pair,
                                     std::string_view lcr_up_seq,
                                     std::string_view lcr_down_seq,
                                     std::string_view interior_seq,
                                     const GenomeHitCounter& hits,
                                     const AnchorBuilderConfig& cfg) {
    PairAnchorTable t;
    t.pair_id = pair.pair_id;

    const std::uint32_t per_segment = cfg.anchors_per_pair / 3;
    std::vector<std::uint64_t> tmp;
    tmp.reserve(cfg.anchors_per_pair + 8);
    CollectCandidates(lcr_up_seq,   cfg.k, cfg.max_genome_hits, hits, per_segment,         tmp);
    CollectCandidates(lcr_down_seq, cfg.k, cfg.max_genome_hits, hits, per_segment,         tmp);
    CollectCandidates(interior_seq, cfg.k, cfg.max_genome_hits, hits,
                       cfg.anchors_per_pair - 2 * per_segment, tmp);

    std::sort(tmp.begin(), tmp.end());
    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
    t.anchor_hashes = std::move(tmp);
    return t;
}

}  // namespace llmap::junction_hunter
