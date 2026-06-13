// LLmap — Provenance resolver implementation.

#include "provenance/provenance_resolver.h"

namespace llmap::provenance {

namespace {

ProvenanceClass MappingToProvenance(MappingConfusion m) noexcept {
    switch (m) {
        case MappingConfusion::Paralog:    return ProvenanceClass::Paralog;
        case MappingConfusion::Numt:       return ProvenanceClass::Numt;
        case MappingConfusion::Pseudogene: return ProvenanceClass::Pseudogene;
        case MappingConfusion::Rdna:       return ProvenanceClass::Rdna;
        case MappingConfusion::None:       return ProvenanceClass::Host;
    }
    return ProvenanceClass::Host;
}

}  // namespace

ReadProvenance ResolveProvenance(const ReadEvidence& ev, const ResolveConfig& cfg) {
    ReadProvenance out;
    out.aligned_bases = ev.aligned_bases;
    out.bioconfounder = ev.bioconfounder;   // Layer-3 carried through unchanged

    // Gather the best competing-origin confounder (the only ones that ask the
    // "better host fit than host?" question). Each candidate = (class, posterior,
    // detail). mapping_confusion already arbitrated paralog/numt/pseudo/rdna.
    ProvenanceClass best_cls = ProvenanceClass::Host;
    float           best_post = 0.0f;
    std::string     best_detail;

    auto consider = [&](ProvenanceClass c, float post, const std::string& detail) {
        if (post > best_post) { best_post = post; best_cls = c; best_detail = detail; }
    };

    if (ev.mapping.kind != MappingConfusion::None) {
        consider(MappingToProvenance(ev.mapping.kind), ev.mapping.confidence,
                 MappingConfusionTag(ev.mapping.kind));
    }
    consider(ProvenanceClass::Exogenous, ev.exogenous_posterior, ev.exogenous_detail);
    consider(ProvenanceClass::CrossIndividual, ev.cross_individual_posterior,
             ev.cross_individual_detail);
    consider(ProvenanceClass::RefArtefact, ev.ref_artefact_posterior,
             ev.ref_artefact_detail);

    // Host-conservative: an origin confounder claims the read only if it beats
    // the host fit by the margin. Otherwise the read stays Host (a real low-VAF
    // somatic must not be absorbed; a host read must not be over-stripped).
    if (best_cls != ProvenanceClass::Host &&
        best_post > ev.host_posterior + cfg.min_margin_over_host) {
        out.origin = best_cls;
        out.posterior = best_post;
        out.detail = best_detail;
        return out;
    }

    // Origin is host-like → apply the structural / technical partition flags.
    if (ev.is_chimera && cfg.chimera_partition) {
        out.origin = ProvenanceClass::Chimera;
        out.posterior = 1.0f;
        out.detail = "chim";
        return out;
    }
    if (ev.is_duplicate && cfg.duplicate_partition) {
        out.origin = ProvenanceClass::Multiplicity;
        out.posterior = 1.0f;
        out.detail = "dup";
        return out;
    }

    out.origin = ProvenanceClass::Host;
    out.posterior = ev.host_posterior;
    return out;
}

}  // namespace llmap::provenance
