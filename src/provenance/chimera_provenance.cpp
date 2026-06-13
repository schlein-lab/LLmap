// LLmap — Chimera provenance detector implementation.

#include "provenance/chimera_provenance.h"

#include <cmath>
#include <vector>

namespace llmap::provenance {

namespace {

// Map the chimera detector's log-prior to a [0,1] confidence (logistic). The
// scale is heuristic; the decision (a hypothesis exists at all) is what the
// detector already thresholded — confidence just orders the calls.
float Confidence(float log_prior) {
    return 1.0f / (1.0f + std::exp(-log_prior));
}

}  // namespace

ChimeraCall ClassifyChimera(std::span<const chimera::AlignedPart> parts,
                            const chimera::VdjLocusMask& vdj_mask,
                            const chimera::ChimeraConfig& cfg) {
    const std::vector<chimera::ChimericHypothesis> hyps =
        chimera::Analyze(parts, vdj_mask, cfg);
    if (hyps.empty()) return ChimeraCall{};  // host / partial — no chimera

    // Best (highest log-prior) hypothesis.
    const chimera::ChimericHypothesis* best = &hyps.front();
    for (const auto& h : hyps) {
        if (h.log_prior > best->log_prior) best = &h;
    }

    ChimeraCall call;
    call.kind = best->kind;
    call.confidence = Confidence(best->log_prior);

    if (best->kind == 'V' || best->vdj_class_switch_detected) {
        // Real immune biology — Layer-3 bioconfounder `bio:vdj`, NOT an artifact.
        call.is_vdj_recombination = true;
    } else if (best->kind == 'I') {
        // Intra-chromosomal split — a long read crossing a real structural
        // breakpoint (deletion/inversion/dup). Host biology → Layer-3 `bio:sv`,
        // not a ligation artifact (matched HG002: HiFi splits >> Illumina).
        call.is_sv_spanning = true;
    } else {
        // Cross-chromosomal split ('X') — ligation/PCR chimera artifact
        // (a true translocation is the rarer case). → Layer-1 `chim`.
        call.is_chimera = true;
    }
    return call;
}

}  // namespace llmap::provenance
