// LLmap — Spliced-chain joiner.
//
// Joins linear sub-chains across intron-sized gaps. The classical
// linear chainer (src/classical/) produces ordered chains of anchors
// with hard gap limits (max_gap_diff=500). This module operates ON
// the linear chainer's output: it walks pairs of adjacent sub-chains,
// checks whether the gap between them looks like an intron (size,
// canonical/annotated junction), and merges them into a SplicedChain.
//
// The classical chainer is NOT modified. DNA-mode pipelines that don't
// run this joiner produce identical output to before.
//
// Two-step structure:
//
//   1. JoinSplicedChains(input_chains, config) → SplicedChainResult
//        Takes a list of LinearSubChain values (one per linear-chainer
//        output chain), pairs adjacent ones whose gap geometry passes
//        the joiner criteria, and emits SplicedChain values with a
//        Junction record per merge boundary.
//
//   2. EmitSplicedCigar(spliced_chain) → CIGAR string with N ops for
//        confirmed junctions and D ops for ambiguous gaps. Wraps the
//        per-sub-chain CIGARs so downstream BAM writers don't need to
//        know about splicing — they get a complete CIGAR.
//
// All gap-acceptance logic is probabilistic, not binary: pass the
// JoinerConfig::min_junction_probability threshold to control how
// permissive the joiner is. Default 0.30 (Plan-Block 6).

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace llmap::mapping {

// ===========================================================================
// Input — one linear sub-chain.
//
// Mirror of the relevant fields from src/classical/chain.h::Chain
// (we don't include that header directly so this module builds
// without dragging the chainer's dep tree).
// ===========================================================================

struct LinearSubChain {
    /// Reference target (chrom or transcript-anchor) — opaque string.
    std::string ref_id;

    /// Genomic span: half-open [start, end).
    std::uint64_t ref_start{0};
    std::uint64_t ref_end{0};

    /// Read span: half-open [start, end).
    std::uint32_t query_start{0};
    std::uint32_t query_end{0};

    /// Chainer score (for tie-breaking).
    std::int32_t score{0};

    /// Strand: '+' or '-'.
    char strand{'+'};

    /// Per-sub-chain CIGAR (e.g. "100M2D50M"). Joiner concatenates
    /// these with intron N ops in EmitSplicedCigar.
    std::string cigar;
};

// ===========================================================================
// Junction — one merge boundary inside a SplicedChain.
// ===========================================================================

struct Junction {
    std::uint64_t donor_ref_pos{0};      // genomic 0-based
    std::uint64_t acceptor_ref_pos{0};
    std::uint32_t query_gap{0};          // gap in read coords (usually 0)

    /// Per-junction probability from KSpliced + JunctionProbability
    /// (Block 5). Range [0, 1].
    float probability{0.0f};

    /// True if the junction was confirmed by canonical splice motif
    /// or annotated junction-DB evidence. Drives CIGAR-N vs CIGAR-D.
    bool is_confirmed{false};
};

// ===========================================================================
// Output — one spliced chain.
// ===========================================================================

struct SplicedChain {
    std::vector<LinearSubChain> sub_chains;   // ordered along strand
    std::vector<Junction> junctions;          // size = sub_chains.size() - 1

    std::int32_t total_score{0};              // sum of sub-chain scores
    char strand{'+'};
    std::string ref_id;                        // shared by all sub_chains
};

struct SplicedChainResult {
    std::vector<SplicedChain> chains;
    /// Sub-chains that the joiner couldn't merge into anything else.
    /// They survive as singleton SplicedChains so the pipeline stays
    /// lossless (no sub-chain ever silently dropped).
    std::uint32_t n_singletons_kept{0};
};

// ===========================================================================
// Configuration.
// ===========================================================================

struct JoinerConfig {
    /// Minimum reference gap to consider as an intron (smaller gaps
    /// are deletions, not introns). Plan-Block 6 default 50 bp.
    std::uint64_t min_intron_bp{50};

    /// Maximum reference gap. Recursive splicing of TITIN-class
    /// genes (Sibley 2016) goes to ~363 kb; allow some headroom.
    std::uint64_t max_intron_bp{1'000'000};

    /// Maximum query gap to tolerate while merging (cDNA at a junction
    /// is contiguous; small gaps from chainer slop are ok).
    std::uint32_t max_query_gap_bp{30};

    /// Minimum joint junction probability (from KSpliced ×
    /// JunctionProbability) for a merge. Lower threshold ⇒ more
    /// aggressive joiner. Plan-Block 6 default 0.30.
    float min_junction_probability{0.30f};
};

// ===========================================================================
// Public API.
// ===========================================================================

/// Join a list of linear sub-chains into SplicedChains.
///
/// `junction_probs[i]` is the precomputed Block-5 junction probability
/// for the pair (input[i], input[i+1]). Must satisfy
/// junction_probs.size() == input.size() - 1, or 0 when input.size()==1.
/// Pass an empty vector to use the default 0.5 (no PWM / DB info).
[[nodiscard]] SplicedChainResult JoinSplicedChains(
    std::span<const LinearSubChain> input,
    std::span<const float> junction_probs,
    const JoinerConfig& cfg = {});

/// Emit a single CIGAR string for the whole spliced chain:
///   - per-sub-chain CIGAR ops concatenated
///   - between consecutive sub-chains: N<intron_len> if junction is
///     confirmed, D<intron_len> otherwise.
/// Always returns a well-formed CIGAR (no zero-length ops).
[[nodiscard]] std::string EmitSplicedCigar(const SplicedChain& sc);

}  // namespace llmap::mapping
