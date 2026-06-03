// LLmap — junction_hunter: shared types for Mode-5 NAHR-junction detection.
//
// Mode-5 detects NAHR breakpoints (read enters LCR_up, jumps to LCR_down
// without traversing the interior) using a k-mer-centric approach with
// multi-k consensus voting. The pipeline runs in two stages:
//
//   Stage-1 ("coarse-locate"): tile the read at k=51, count anchor hits
//   against each NAHR-pair's pre-extracted anchor table. Reads with ≥2
//   anchor hits for a pair are forwarded to Stage-2. ~99.5 % of reads
//   are dropped here (canonical genomic reads not in any LCR locus).
//
//   Stage-2 ("multi-k convergence"): at each read position, look up k=21,
//   31, 51, 71, 101 against the pair's LCR_up / LCR_down / interior
//   sub-indices. Majority-of-5 consensus per position; outside-in
//   monotonicity across LEFT (LCR_up) and RIGHT (LCR_down) halves;
//   breakpoint = the position where the consensus class flips.
//
// minimap2 is explicitly NOT used at any stage — see feedback memo
// `feedback_minimap2_nahr_forbidden`. The geometric pattern enforced
// here (monotonic outer→inner walk in each LCR copy + zero interior
// hits + multi-k agreement) cannot be produced by chimeric alignment
// artefacts, which is the failure mode minimap2 has at LCR_id <99 %.

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace llmap::junction_hunter {

/// NAHR-pair coordinates as loaded from the genome-wide pair TSV
/// (`genome_wide_nahr_2kb_pairs.tsv`). All coordinates are GRCh38
/// 0-based half-open intervals on the same chromosome.
struct NahrPair {
    std::string pair_id;        ///< Stable pair identifier (e.g. NAHR_02371_chr22_20M)
    std::string chrom;          ///< Chromosome (chr1..chr22, chrX, chrY)
    std::uint64_t lcr_up_start{0};
    std::uint64_t lcr_up_end{0};
    std::uint64_t lcr_down_start{0};
    std::uint64_t lcr_down_end{0};
    std::uint64_t interior_start{0};
    std::uint64_t interior_end{0};
    float lcr_identity{0.0f};   ///< Fraction identity between LCR_up and LCR_down
    std::uint32_t interior_kb{0};
};

/// Outcome of classifying a single read against a single NAHR-pair.
enum class JunctionCall : std::uint8_t {
    Unmapped = 0,        ///< No hits in pair regions
    CanonicalUp,         ///< Monotonic walk through LCR_up only
    CanonicalDown,       ///< Monotonic walk through LCR_down only
    CanonicalInterior,   ///< Hits in interior — locus is intact in this read
    JunctionReal,        ///< Outside-in convergence: real NAHR breakpoint
    ChimeraArtifact,     ///< Non-monotonic / mixed pattern — sequencer fluke
    ParalogAmbiguous,    ///< Multi-k consensus conflict — cannot decide
};

const char* JunctionCallName(JunctionCall c) noexcept;

/// Single record per (read, pair) test. Aggregable across cohorts.
struct JunctionRecord {
    std::string read_id;
    std::string pair_id;
    std::uint32_t n_kmer_total{0};        ///< Positions in the read evaluated
    std::uint32_t n_consensus_up{0};
    std::uint32_t n_consensus_dn{0};
    std::uint32_t n_consensus_in{0};
    std::uint32_t n_ambiguous{0};         ///< Multi-k conflict positions
    float up_monotonicity{0.0f};          ///< Spearman ρ; ≥0.95 = monotonic
    float dn_monotonicity{0.0f};
    std::uint32_t breakpoint_read_pos{0}; ///< Index where consensus flips
    float breakpoint_quality{0.0f};       ///< Multi-k agreement in ±50 bp window
    JunctionCall call{JunctionCall::Unmapped};
};

/// Multi-k configuration. The five k-lengths are anchored at the same
/// read position; consensus requires majority (≥3 of 5) agreement.
struct MultiKConfig {
    std::array<std::uint8_t, 5> k_values{21, 31, 51, 71, 101};
    /// Minimum number of k-values that must agree on the hit class at a
    /// given read position for it to count as a consensus position.
    std::uint8_t consensus_min{3};
    /// Monotonicity threshold (Spearman ρ on hit positions vs read
    /// positions); below this, the LCR-half is rejected as non-monotonic.
    float monotonicity_min{0.95f};
    /// Minimum PSV-switch count to call JunctionReal (per spec §1.4).
    std::uint32_t min_psv_switches{3};
};

/// Stage-1 anchor table — pre-extracted k=51 k-mers unique to a single
/// NAHR-pair (≤5 genome-wide hits). 100 anchors evenly sampled across
/// (LCR_up ∪ LCR_down ∪ interior); only their hashes are kept here for
/// fast counting in the coarse-locate phase.
struct PairAnchorTable {
    std::string pair_id;
    std::vector<std::uint64_t> anchor_hashes;
};

}  // namespace llmap::junction_hunter
