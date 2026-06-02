// LLmap — Transcript-Mode dispatch enum.
//
// CLI flag --mode {auto,transcript,reads,assembly,reads_vs_assembly}
// resolves to one of these values. The Stage-2 pipeline + the
// Multi-Signal Fusion engine + the spliced-chain joiner all branch
// on this.
//
// `Auto` triggers the input sniffer (Block 9 CLI surface): a quick
// pass over the first ~10 k FASTQ records picks Transcript if filename
// matches FLNC/isoseq/cdna OR header tokens / read-length distribution
// / polyA signal suggest mRNA input. Otherwise resolves to GenomeReads
// for backward compatibility with existing DNA mode.

#pragma once

#include <cstdint>
#include <optional>
#include <string_view>

namespace llmap::core {

enum class TranscriptMode : std::uint8_t {
    Auto = 0,
    Transcript,
    GenomeReads,
    Assembly,
    ReadsVsAssembly,
};

[[nodiscard]] const char* TranscriptModeName(TranscriptMode m) noexcept;
[[nodiscard]] std::optional<TranscriptMode> ParseTranscriptMode(
    std::string_view s) noexcept;

}  // namespace llmap::core
