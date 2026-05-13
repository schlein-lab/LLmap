#pragma once

// Internal implementation header for minimizer_index split files.
// Not for external use.

#include "classical/minimizer_index.h"

#include <unordered_map>
#include <vector>

namespace llmap::classical {

// Internal index storage
struct MinimizerIndex::Impl {
    // Hash table: minimizer hash -> list of (ref_id, position)
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> index;

    // Indexed sequences
    std::vector<IndexedSequence> sequences;

    // Statistics
    MinimizerStats stats;
};

// Builder internal storage
struct MinimizerIndex::Builder::Impl {
    MinimizerConfig config;
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};

namespace detail {

// Encode base to 2-bit
constexpr uint8_t EncodeBase(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;  // N and others treated as A
    }
}

// Check if base is valid (not N)
constexpr bool IsValidBase(char c) {
    return c == 'A' || c == 'a' ||
           c == 'C' || c == 'c' ||
           c == 'G' || c == 'g' ||
           c == 'T' || c == 't';
}

// Complement of 2-bit encoded base
constexpr uint8_t Complement(uint8_t base) {
    return 3 - base;
}

}  // namespace detail

}  // namespace llmap::classical
