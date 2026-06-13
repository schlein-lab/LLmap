// LLmap — pangenome M(pos) prior → TE spread-mass bridge implementation.

#include "provenance/te_population_prior.h"

#include <utility>

namespace llmap::provenance {

void ApplyPopulationPrior(std::vector<TePlacement>& cands,
                          const core::PangenomeMappability& track,
                          float neutral_for_unseen) {
    for (auto& c : cands) {
        const std::uint32_t seen = track.SampleCount(c.chrom, c.pos);
        if (seen > 0) {
            // Reproducible mappability across the pangenome at this copy.
            c.population_support = track.PopulationSupport(c.chrom, c.pos);
        } else {
            // Never observed in the pangenome run → caller's policy. A negative
            // value keeps the te_spread_mass "no baseline" sentinel (reference-only
            // fallback); a non-negative value treats unseen as that neutral mass.
            c.population_support = neutral_for_unseen;
        }
    }
}

SpreadMass ResolveWithPopulationPrior(std::vector<TePlacement> cands,
                                      const core::PangenomeMappability& track,
                                      float temperature,
                                      float neutral_for_unseen) {
    ApplyPopulationPrior(cands, track, neutral_for_unseen);
    return ResolveTeSpreadMass(cands, temperature);
}

std::vector<int> FlagBlurryOnlyHere(const std::vector<TePlacement>& cands,
                                    const core::PangenomeMappability& track,
                                    float m_sharp) {
    std::vector<int> anomalies;
    for (std::size_t i = 0; i < cands.size(); ++i) {
        // Only meaningful where the pangenome actually has an opinion.
        if (track.SampleCount(cands[i].chrom, cands[i].pos) == 0) continue;
        const float m = track.PopulationSupport(cands[i].chrom, cands[i].pos);
        // The pangenome says this locus is reproducibly sharp, yet THIS read did
        // not anchor here → sample-specific, not a real multi-copy ambiguity.
        if (m >= m_sharp && !cands[i].has_unique_anchor)
            anomalies.push_back(static_cast<int>(i));
    }
    return anomalies;
}

}  // namespace llmap::provenance
