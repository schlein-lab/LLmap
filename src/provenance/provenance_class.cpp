// LLmap — Read provenance class implementation.

#include "provenance/provenance_class.h"

#include <array>
#include <string_view>

namespace llmap::provenance {

namespace {
// Index by ProvenanceClass value; keep in lock-step with the enum order.
constexpr std::array<std::string_view,
                     static_cast<std::size_t>(ProvenanceClass::Count)> kTags = {
    "host", "exo", "xindiv", "para", "numt", "pseudo", "rdna", "mei",
    "refartefact", "chim", "xsample", "dup",
};
}  // namespace

const char* ProvenanceClassTag(ProvenanceClass c) noexcept {
    const auto i = static_cast<std::size_t>(c);
    return i < kTags.size() ? kTags[i].data() : "host";
}

std::optional<ProvenanceClass> ParseProvenanceClass(std::string_view tag) noexcept {
    for (std::size_t i = 0; i < kTags.size(); ++i) {
        if (tag == kTags[i]) return static_cast<ProvenanceClass>(i);
    }
    return std::nullopt;
}

bool IsOriginBucket(ProvenanceClass c) noexcept {
    switch (c) {
        // A: competing-reference origin buckets (a *different sequence*).
        case ProvenanceClass::Exogenous:
        case ProvenanceClass::CrossIndividual:
        case ProvenanceClass::Paralog:
        case ProvenanceClass::Numt:
        case ProvenanceClass::Pseudogene:
        case ProvenanceClass::Rdna:
        case ProvenanceClass::MobileElement:
            return true;
        default:
            // Host, RefArtefact, Chimera, CrossSample, Multiplicity: partition
            // members but not "another sequence" buckets.
            return false;
    }
}

const char* BioConfounderTag(BioConfounder f) noexcept {
    switch (f) {
        case BioConfounder::Vdj:            return "bio:vdj";
        case BioConfounder::Shm:            return "bio:shm";
        case BioConfounder::ClassSwitch:    return "bio:classswitch";
        case BioConfounder::GeneConversion: return "bio:geneconv";
        case BioConfounder::MtHeteroplasmy: return "bio:mthet";
        case BioConfounder::Mosaicism:      return "bio:mosaic";
        case BioConfounder::Imprinting:     return "bio:imprint";
        case BioConfounder::SvSpanning:     return "bio:sv";
        case BioConfounder::None:           return "bio:none";
    }
    return "bio:none";
}

}  // namespace llmap::provenance
