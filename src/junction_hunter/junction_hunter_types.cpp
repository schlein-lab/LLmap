// LLmap — junction_hunter: type-helper implementations.

#include "junction_hunter/junction_hunter_types.h"

namespace llmap::junction_hunter {

const char* JunctionCallName(JunctionCall c) noexcept {
    switch (c) {
        case JunctionCall::Unmapped:           return "unmapped";
        case JunctionCall::CanonicalUp:        return "canonical_up";
        case JunctionCall::CanonicalDown:      return "canonical_down";
        case JunctionCall::CanonicalInterior:  return "canonical_interior";
        case JunctionCall::JunctionReal:       return "junction_real";
        case JunctionCall::ChimeraArtifact:    return "chimera_artifact";
        case JunctionCall::ParalogAmbiguous:   return "paralog_ambiguous";
    }
    return "unknown";
}

}  // namespace llmap::junction_hunter
