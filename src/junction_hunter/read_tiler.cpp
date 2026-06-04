// LLmap — junction_hunter: read-tiler implementation.

#include "junction_hunter/read_tiler.h"

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

/// Canonical 64-bit hash (forward / reverse-complement min).
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

}  // namespace

ReadTiling TileRead(std::string_view read_seq, const MultiKConfig& cfg) {
    ReadTiling out;
    out.k_values = cfg.k_values;
    for (std::size_t ki = 0; ki < cfg.k_values.size(); ++ki) {
        const std::uint8_t k = cfg.k_values[ki];
        if (read_seq.size() < k) continue;
        // First-pass scan: skip windows containing non-ACGT bases.
        // For simplicity we re-validate each window — O(N·k); acceptable
        // because the multi-k pass is bounded by the longest k (101).
        const std::size_t n = read_seq.size();
        for (std::size_t i = 0; i + k <= n; ++i) {
            bool clean = true;
            for (std::size_t j = 0; j < k; ++j) {
                if (!IsAcgt(read_seq[i + j])) { clean = false; break; }
            }
            if (!clean) continue;
            out.per_k_hashes[ki].push_back({HashKmer(read_seq.substr(i, k)), static_cast<std::uint32_t>(i)});
        }
    }
    return out;
}

}  // namespace llmap::junction_hunter
