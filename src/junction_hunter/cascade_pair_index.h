// LLmap — junction_hunter: cascade pair k-mer index (variable-tier).
//
// Cascade architecture: small k first for cheap membership, climb to
// large k for PSV-grade paralog disambiguation. The number of tiers
// and the k value per tier are configurable rather than fixed at five.
//
//   Long-read preset (HiFi/ONT): k = {11, 17, 25, 35, 51, 71, 101, 125}
//   Short-read preset (Illumina): k = {11, 13, 15, 17, 19, 21, 25}
//
// Tier 0 (smallest k) is built for every pair up-front; it is the
// MEMBERSHIP filter and rarely disambiguates paralogs at LCR identity
// <99 %. Higher tiers are built lazily, per pair, the first time tier
// 0 nominates a read for promotion. Memory cost is dominated by tier 0
// (~30-80 MB for the full 3500-pair panel at k=11/17), with higher
// tiers added only for the small surviving set.

#pragma once

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/pair_kmer_index.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::junction_hunter {

struct CascadePairIndex {
    NahrPair pair_meta;
    /// k-mer class table per tier. tier_class[0] is always built;
    /// tier_class[i>=1] are built lazily on demand.
    std::vector<KmerClassMap> tier_class;
    std::vector<bool> tier_built;
    /// k values per tier — copied from the config so the index is
    /// self-describing.
    std::vector<std::uint8_t> k_values;
};

/// Build the bottom-tier (smallest k) table for a pair. Returns the
/// number of unique k-mers indexed at tier 0.
std::size_t BuildCascadeTier1(const NahrPair& pair,
                              std::string_view lcr_up_seq,
                              std::string_view lcr_down_seq,
                              std::string_view interior_seq,
                              CascadePairIndex& out);

/// Build a specific higher tier (1..N-1). Idempotent. Returns the
/// number of unique k-mers added.
std::size_t BuildCascadeTierN(std::size_t tier_idx,
                              std::string_view lcr_up_seq,
                              std::string_view lcr_down_seq,
                              std::string_view interior_seq,
                              CascadePairIndex& idx);

}  // namespace llmap::junction_hunter
