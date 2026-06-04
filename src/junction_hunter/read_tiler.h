// LLmap — junction_hunter: read tiler.
//
// Tiles a read sequence at five k-lengths simultaneously, ANCHORED at
// the same read position. The multi-k design lets the consensus stage
// override a paralog-ambiguous shorter-k call with a longer-k unique
// hit (and vice versa, override a low-sensitivity longer-k drop with
// the shorter-k coverage).

#pragma once

#include "junction_hunter/junction_hunter_types.h"

#include <array>
#include <cstdint>
#include <string_view>
#include <vector>

namespace llmap::junction_hunter {

/// k-mer hash + the originating read-position. Position is in read
/// coordinates (0-based, 5'→3').
struct ReadKmerHit {
    std::uint64_t hash{0};
    std::uint32_t read_pos{0};
};

/// Per-read tiling output. Five parallel hash streams, one per k.
struct ReadTiling {
    std::array<std::vector<ReadKmerHit>, 5> per_k_hashes;
    std::array<std::uint8_t, 5> k_values{21, 31, 51, 71, 101};
};

/// Tile a read at the five configured k-lengths. For each k, emit one
/// (hash, position) per read position where the window is fully clean
/// (ACGT only). N-bases break the window — that window is skipped at
/// the affected k but the others may still emit, which is intentional:
/// a single sequencing N near the read end shouldn't kill all k.
ReadTiling TileRead(std::string_view read_seq,
                    const MultiKConfig& cfg = {});

}  // namespace llmap::junction_hunter
