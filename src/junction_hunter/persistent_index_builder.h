// LLmap — junction_hunter: persistent k-mer index builder.
//
// Build-once tool that reads the NAHR-pair panel + reference FASTA,
// extracts LCR_up / LCR_down / interior sequences for every pair,
// enumerates all k-mer hashes at each cascade tier, sorts globally,
// and writes one tier_k{N}.bin file per k value.
//
// Memory footprint during build is proportional to the number of
// entries at the largest tier; expect 30-80 GB for the full 3500-pair
// panel at k=21..125. Run on a bigmem node.

#pragma once

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/persistent_index.h"

#include <cstdint>
#include <string>
#include <vector>

namespace llmap::junction_hunter {

struct BuildOptions {
    std::string out_dir;                 ///< output directory; will hold tier_k{N}.bin
    std::vector<std::uint8_t> k_values;  ///< cascade tiers to build
    std::string panel_path;              ///< for sha256 stamp
    std::string reference_path;          ///< for sha256 stamp
    bool verbose{false};
};

struct BuildStats {
    std::vector<std::uint64_t> entries_per_tier;
    std::vector<std::uint64_t> bytes_per_tier;
    double seconds_total{0.0};
};

/// Build all tier files. Allocates one big std::vector<PersistentIndexEntry>
/// per tier (sequentially, not all at once) — sort in place, write to disk.
/// Pair sequences are streamed from `ref` once per tier to keep RSS bounded.
bool BuildPersistentIndex(const std::vector<NahrPair>& pairs,
                          const std::string& reference_path,
                          const BuildOptions& opts,
                          BuildStats& stats,
                          std::string& err);

}  // namespace llmap::junction_hunter
