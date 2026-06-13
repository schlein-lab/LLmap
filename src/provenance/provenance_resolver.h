// LLmap — Provenance resolver: detector evidence → one ReadProvenance.
//
// The integration spine of the provenance mode. Each detector (mapping_confusion
// = paralog/numt/pseudo/rdna, the contaminant-panel pass = exogenous, the
// chimera detector, the duplicate flagger, the reference-artefact check) reports
// its evidence into `ReadEvidence`; the resolver collapses it into a single
// `ReadProvenance` — exactly one Layer-1 origin class (the Σ-partition) plus the
// Layer-3 bioconfounder overlay carried through unchanged.
//
// Host-conservative decision (Agent 1's calibration caveat): a confounder class
// claims the read ONLY if its posterior beats the host posterior by a margin
// (optionally weighted by per-class prior prevalence). Otherwise the read stays
// `Host`. This is the flag-but-don't-over-flag balance — a real low-VAF somatic
// must not be absorbed into an artefact bucket, and a host read must not be
// stripped just because a contaminant bucket exists.

#pragma once

#include "provenance/mapping_confusion.h"
#include "provenance/provenance_class.h"

#include <cstdint>
#include <string>

namespace llmap::provenance {

// Per-read evidence gathered from the detector passes. Absent detectors leave
// their fields at the inert defaults (posterior 0 / false) → no false claim.
struct ReadEvidence {
    // --- Layer-1 origin candidates ---
    MappingConfusionCall mapping{};            // paralog/numt/pseudo/rdna + post.
    float        exogenous_posterior{0.0f};    // contaminant-panel hit (opt-in)
    std::string  exogenous_detail;             // "exo:mouse", "exo:phix"
    float        cross_individual_posterior{0.0f};  // 2nd mammalian genome
    std::string  cross_individual_detail;      // "xindiv:mouse","xindiv:swap"
    float        ref_artefact_posterior{0.0f}; // reference-side artefact
    std::string  ref_artefact_detail;          // "ref:ancestry","ref:assembly"
    bool         is_chimera{false};            // Block-7 chimera detector
    bool         is_duplicate{false};          // PCR/optical duplicate

    // --- read context ---
    float         host_posterior{1.0f};        // how well the read fits host
    std::uint32_t aligned_bases{0};

    // --- Layer-3 overlay (carried through unchanged; never partitions) ---
    std::uint16_t bioconfounder{0};            // OR of BioConfounder flags
};

struct ResolveConfig {
    // A confounder claims the read only if posterior > host_posterior + margin.
    float min_margin_over_host{0.15f};
    // Chimera / duplicate are structural/technical and claim regardless of the
    // host margin (they are not a "better host fit" question).
    bool  chimera_partition{true};
    bool  duplicate_partition{true};
};

// Resolve one read's evidence into its provenance (Layer-1 origin + Layer-3).
[[nodiscard]] ReadProvenance ResolveProvenance(const ReadEvidence& ev,
                                               const ResolveConfig& cfg = {});

}  // namespace llmap::provenance
