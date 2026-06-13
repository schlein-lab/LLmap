// LLmap — WaveCollapse merge-ambiguity guard (read-vs-read second merge gate).
//
// minimizer_cluster links two reads on shared-canonical-minimizer count + overlap
// span + a self-frequency cap. That is necessary but not sufficient: two reads
// from different copies of a repeat/paralog family can share many minimizers
// WITHOUT a real overlap — the shared anchors scatter across many offsets instead
// of agreeing on one. This guard is the WaveCollapse collapse-criterion applied
// read-vs-read: a TRUE overlap makes the anchor mass collapse onto ONE diagonal
// (a single consistent offset posB - posA); a spurious match leaves the mass
// spread (high offset entropy) → do NOT merge (keep the reads separate; the
// pipeline maps them individually — lossless, nothing dropped).
//
// This is a pure predicate over the per-anchor offset list, decoupled from the
// clustering itself so it can be tested in isolation and wired in opt-in (the
// default clustering behaviour is unchanged until a caller enables it).

#pragma once

#include <cstdint>
#include <span>

namespace llmap::mapper {

struct AmbiguityConfig {
    std::int64_t  offset_tolerance{20};        // anchors within ±tol ⇒ same diagonal
    double        min_dominant_fraction{0.60}; // collapse threshold (clear overlap)
    std::uint32_t min_anchors{3};              // fewer ⇒ cannot judge ⇒ not clear
};

struct DiagonalClarity {
    std::uint32_t shared{0};            // total shared anchors considered
    std::uint32_t on_diagonal{0};       // anchors on the dominant diagonal (±tol)
    double        dominant_fraction{0}; // on_diagonal / shared
    double        offset_entropy{0};    // normalised Shannon entropy of offsets [0,1]
    bool          clear{false};         // mass collapsed ⇒ safe to merge
};

// Assess whether the shared anchors collapse onto one diagonal. `offsets` is the
// list of per-anchor offsets (posB - posA) for the shared canonical minimizers of
// a read pair. `clear` ⇒ a real, unambiguous overlap (safe to merge).
[[nodiscard]] DiagonalClarity AssessDiagonal(std::span<const std::int64_t> offsets,
                                             const AmbiguityConfig& cfg = {});

}  // namespace llmap::mapper
