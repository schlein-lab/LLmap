// LLmap — Transcript-Mode dispatch implementation.

#include "core/transcript_mode.h"

#include <cstdint>

namespace llmap::core {

const char* TranscriptModeName(TranscriptMode m) noexcept {
    switch (m) {
        case TranscriptMode::Auto:            return "auto";
        case TranscriptMode::Transcript:      return "transcript";
        case TranscriptMode::GenomeReads:     return "reads";
        case TranscriptMode::Assembly:        return "assembly";
        case TranscriptMode::ReadsVsAssembly: return "reads_vs_assembly";
    }
    return "auto";
}

std::optional<TranscriptMode> ParseTranscriptMode(std::string_view s) noexcept {
    if (s == "auto")              return TranscriptMode::Auto;
    if (s == "transcript")        return TranscriptMode::Transcript;
    if (s == "reads" || s == "genome_reads")
                                  return TranscriptMode::GenomeReads;
    if (s == "assembly")          return TranscriptMode::Assembly;
    if (s == "reads_vs_assembly") return TranscriptMode::ReadsVsAssembly;
    return std::nullopt;
}

}  // namespace llmap::core
