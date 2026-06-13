// LLmap — bridge: pangenome M(pos) prior → TE spread-mass placements.
//
// The Operator's directive (2026-06-13) made concrete: a TE-copy placement gains
// mass when it is the consistently-mappable one across the pangenome, and keeps an
// honest spread when blurriness is reproducible everywhere. te_spread_mass already
// consumes a per-placement `population_support`; this is the thin layer that fills
// it from the cross-sample mappability track M(pos) — the second free aggregate of
// the same HPRC pangenome run that produced the provenance baseline.
//
// Kept separate from te_spread_mass (pure stdlib resolver) because it depends on
// core::PangenomeMappability; provenance already links core, so this is the natural
// home for the bridge.

#pragma once

#include "core/pangenome_mappability.h"
#include "provenance/te_spread_mass.h"

#include <vector>

namespace llmap::provenance {

// Fill each candidate's population_support from the pangenome M(pos) track, then
// the caller runs ResolveTeSpreadMass with the prior active.
//
//   * window seen across samples → M ∈ [0,1] (reproducible mappability) → prior.
//   * window never seen → `neutral_for_unseen`. Pass a NEGATIVE value (default
//     -1) to leave such placements at the te_spread_mass "no baseline" sentinel,
//     so a locus absent from the pangenome falls back to reference-only behaviour
//     instead of being pushed toward a fabricated 0.5. Pass 0.5 to treat unseen as
//     neutral-but-present (a uniform 0.5 across all candidates is a no-op in the
//     softmax anyway; it only bites when SOME candidates are seen and others not).
void ApplyPopulationPrior(std::vector<TePlacement>& cands,
                          const core::PangenomeMappability& track,
                          float neutral_for_unseen = -1.0f);

// Convenience: fill the prior then resolve in one call.
[[nodiscard]] SpreadMass ResolveWithPopulationPrior(
    std::vector<TePlacement> cands,
    const core::PangenomeMappability& track,
    float temperature = 1.0f,
    float neutral_for_unseen = -1.0f);

// The anomaly signal the Operator asked for, kept distinct from the prior: a
// placement whose own alignment is weak exactly where the pangenome says the locus
// is reproducibly sharp (M high) is sample-specific — an artefact/contamination
// candidate, NOT a real ambiguity. Returns the indices of such anomalous
// placements (M(pos) ≥ m_sharp but this read failed to anchor there).
[[nodiscard]] std::vector<int> FlagBlurryOnlyHere(
    const std::vector<TePlacement>& cands,
    const core::PangenomeMappability& track,
    float m_sharp = 0.8f);

}  // namespace llmap::provenance
