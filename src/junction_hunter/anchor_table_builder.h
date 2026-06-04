// LLmap — junction_hunter: Stage-1 anchor-table builder.
//
// For each NahrPair, extract ~100 evenly-spaced k=51 k-mers from the
// pair's reference sequence (LCR_up ∪ LCR_down ∪ interior) and keep
// only those that occur at most a few times genome-wide. The resulting
// PairAnchorTable is tiny (~MB total across 3518 pairs) and lets Stage-1
// reject ~99.5 % of reads via cheap hash-membership tests.
//
// The "genome-wide hit count" check requires a reference k-mer index;
// to keep this module modular we accept a callable to test it instead
// of pulling in the classical minimizer index directly.

#pragma once

#include "junction_hunter/junction_hunter_types.h"

#include <cstdint>
#include <functional>
#include <string_view>
#include <vector>

namespace llmap::junction_hunter {

struct AnchorBuilderConfig {
    /// k for anchor k-mers. Must match Stage-1 lookup k.
    std::uint8_t k{51};
    /// Maximum genome-wide hits a k-mer may have to be eligible as
    /// anchor. 5 is a good default: tolerates segdup-internal near-
    /// duplicates within the pair itself, rejects general repeats.
    std::uint32_t max_genome_hits{5};
    /// Target number of anchors to keep per pair (evenly sampled across
    /// LCR_up ∪ LCR_down ∪ interior).
    std::uint32_t anchors_per_pair{100};
};

/// Callable signature: (kmer-hash) → genome-wide hit count.
/// Implemented by the classical minimizer index or any other module
/// that can answer this question.
using GenomeHitCounter = std::function<std::uint32_t(std::uint64_t)>;

/// Build an anchor table for one pair.
///
/// Required inputs:
///   pair              the NAHR-pair coordinates
///   lcr_up_seq        sequence covering [lcr_up_start, lcr_up_end)
///   lcr_down_seq      sequence covering [lcr_down_start, lcr_down_end)
///   interior_seq      sequence covering [interior_start, interior_end)
///                     — may be empty; the interior is often skipped
///                     for anchor sampling.
///   hits              genome-wide hit-count callable
///   cfg               builder configuration
///
/// Returns the populated PairAnchorTable; anchor_hashes is sorted.
PairAnchorTable BuildPairAnchorTable(const NahrPair& pair,
                                     std::string_view lcr_up_seq,
                                     std::string_view lcr_down_seq,
                                     std::string_view interior_seq,
                                     const GenomeHitCounter& hits,
                                     const AnchorBuilderConfig& cfg = {});

}  // namespace llmap::junction_hunter
