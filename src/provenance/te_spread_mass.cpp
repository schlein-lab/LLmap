// LLmap — TE-scale WaveCollapse spread-mass implementation.

#include "provenance/te_spread_mass.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace llmap::provenance {

namespace {
std::uint8_t HonestMapq(float best_posterior) {
    // MAPQ = -10*log10(P(placement wrong)) = -10*log10(1 - best_posterior),
    // capped at 60. N near-equal copies ⇒ best≈1/N ⇒ small MAPQ (truthful), not a
    // faked 60. This is the whole point: honest ambiguity, not a fabricated hit.
    const float perr = std::max(1.0e-6f, 1.0f - best_posterior);
    const float q = -10.0f * std::log10(perr);
    return static_cast<std::uint8_t>(std::min(60.0f, std::max(0.0f, q)));
}
}  // namespace

SpreadMass ResolveTeSpreadMass(const std::vector<TePlacement>& cands, float temperature) {
    SpreadMass r;
    const std::size_t n = cands.size();
    if (n == 0) { r.is_ambiguous = true; return r; }
    r.posterior.assign(n, 0.0f);

    // 1. A unique anchor (a read flank reaching out of the TE) disambiguates →
    //    collapse to the best-scoring anchored placement. Confident, real locus.
    int best_anchor = -1;
    std::int32_t best_anchor_score = std::numeric_limits<std::int32_t>::min();
    for (std::size_t i = 0; i < n; ++i) {
        if (cands[i].has_unique_anchor && cands[i].align_score > best_anchor_score) {
            best_anchor = static_cast<int>(i);
            best_anchor_score = cands[i].align_score;
        }
    }
    if (best_anchor >= 0) {
        r.posterior[static_cast<std::size_t>(best_anchor)] = 1.0f;
        r.collapsed = best_anchor;
        r.confidence = 1.0f;
        r.mapq = 60;
        r.is_ambiguous = false;
        return r;
    }

    // 2. No unique anchor → keep mass spread over the copies (softmax of scores),
    //    weighted by the POPULATION PRIOR when a per-locus mappability baseline is
    //    available: a placement that is consistently the sensible one across the
    //    pangenome gains mass; one rarely right across the population loses it.
    //    Absent baseline (population_support < 0) → reference-only behaviour, so a
    //    locus that is blurry *everywhere* keeps an honest even spread (real
    //    ambiguity), while the caller separately flags blurry-only-here anomalies.
    const float T = std::max(1.0e-3f, temperature);
    std::int32_t max_s = std::numeric_limits<std::int32_t>::min();
    for (const auto& c : cands) max_s = std::max(max_s, c.align_score);
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double expo = (static_cast<double>(cands[i].align_score) -
                       static_cast<double>(max_s)) / T;
        if (cands[i].population_support >= 0.0f)
            expo += std::log(std::max(1.0e-3f, cands[i].population_support));
        const double e = std::exp(expo);
        r.posterior[i] = static_cast<float>(e);
        sum += e;
    }
    float best_post = 0.0f;
    for (std::size_t i = 0; i < n; ++i) {
        r.posterior[i] = static_cast<float>(r.posterior[i] / sum);
        best_post = std::max(best_post, r.posterior[i]);
    }

    // confidence = 1 - normalised entropy (1 = one copy dominates; 0 = N equal).
    double entropy = 0.0;
    for (const float p : r.posterior) if (p > 0.0f) entropy -= p * std::log(p);
    const double norm = (n > 1) ? std::log(static_cast<double>(n)) : 1.0;
    r.confidence = (n > 1) ? static_cast<float>(1.0 - entropy / norm) : 1.0f;
    r.mapq = HonestMapq(best_post);
    r.is_ambiguous = true;
    return r;
}

}  // namespace llmap::provenance
