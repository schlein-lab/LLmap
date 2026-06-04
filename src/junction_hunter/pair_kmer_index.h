// LLmap — junction_hunter: per-pair multi-k k-mer index.
//
// For one NAHR-pair, holds segregated hash tables — one per k in
// {21, 31, 51, 71, 101} — mapping k-mer hashes to the genomic class
// {LcrUp, LcrDown, Interior} of the position where they occur. A k-mer
// that occurs in MULTIPLE classes (e.g. shared LCR_up / LCR_down k-mer,
// extremely common at LCR identity ≥98 %) is tagged as Ambiguous; the
// consensus stage will need a longer k at that position to disambiguate.

#pragma once

#include "junction_hunter/junction_hunter_types.h"

#include <array>
#include <cstdint>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::junction_hunter {

enum class LocusClass : std::uint8_t {
    None = 0,
    LcrUp,
    LcrDown,
    Interior,
    Ambiguous,   ///< k-mer found in more than one of LcrUp/LcrDown/Interior at this k
};

const char* LocusClassName(LocusClass c) noexcept;

/// Per-k lookup table. Values: class tag for that hash (most common
/// case is None when hash not present; full-bucket structure inflates
/// memory for marginal gain).
using KmerClassMap = std::unordered_map<std::uint64_t, LocusClass>;

/// All five per-k tables for one NAHR-pair.
struct PairKmerIndex {
    std::string pair_id;
    std::array<KmerClassMap, 5> per_k_class;
    std::array<std::uint8_t, 5> k_values{21, 31, 51, 71, 101};
};

struct PairIndexBuildStats {
    std::array<std::size_t, 5> unique_kmers_per_k{};
    std::array<std::size_t, 5> ambiguous_per_k{};
    std::size_t total_unique{0};
};

/// Build the multi-k index from the pair's LCR_up / LCR_down / interior
/// sequences. Each unique k-mer hash gets tagged with the FIRST class
/// it was seen in; if a later class re-hits it, the tag is upgraded
/// to Ambiguous.
PairIndexBuildStats BuildPairKmerIndex(const NahrPair& pair,
                                       std::string_view lcr_up_seq,
                                       std::string_view lcr_down_seq,
                                       std::string_view interior_seq,
                                       const MultiKConfig& cfg,
                                       PairKmerIndex& out);

}  // namespace llmap::junction_hunter
