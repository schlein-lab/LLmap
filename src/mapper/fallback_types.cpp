// LLmap — SegDup mapper fallback hierarchy: shared-type definitions.

#include "mapper/fallback_types.h"

namespace llmap::mapper::fallback {

const char* StageIdName(StageId id) noexcept {
    switch (id) {
        case StageId::relaxed_mismatch:     return "relaxed_mismatch";
        case StageId::chain_only:           return "chain_only";
        case StageId::multi_position:       return "multi_position";
        case StageId::llm_checkpoint:       return "llm_checkpoint";
        case StageId::novel_haplotype_flag: return "novel_haplotype_flag";
    }
    return "unknown";
}

}  // namespace llmap::mapper::fallback
