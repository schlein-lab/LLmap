// LLmap — cmd_align reporting module.
// Extracted from cmd_align.cpp to keep each module < 400 LOC.

#include "cli/cmd_align_internal.h"

#include <cstdio>

namespace llmap::cli::align_internal {

void PrintAlignmentSummary(
    const AlignArgs& args,
    const classical::ClassicalPipelineStats& stats,
    std::size_t total_reads,
    std::size_t n_mapped,
    std::size_t n_unmapped,
    float align_time_ms,
    float total_time_ms) {

    float mapping_rate = total_reads > 0
        ? static_cast<float>(n_mapped) / static_cast<float>(total_reads)
        : 0.0f;

    std::printf("Alignment complete:\n");
    std::printf("  Input reads:    %zu\n", total_reads);
    std::printf("  Mapped:         %zu (%.1f%%)\n",
                n_mapped, 100.0f * mapping_rate);
    std::printf("  Unmapped:       %zu\n", n_unmapped);
    std::printf("  Pipeline drop breakdown:\n");
    std::printf("    minimizer hits found: %zu (avg %.1f/read)\n",
                stats.total_hits,
                static_cast<float>(stats.total_hits) /
                    std::max<size_t>(1, total_reads));
    std::printf("    chains formed:        %zu (avg %.2f/read)\n",
                stats.total_chains,
                static_cast<float>(stats.total_chains) /
                    std::max<size_t>(1, total_reads));
    std::printf("    extensions accepted:  %zu\n", stats.total_extensions);
    std::printf("    rejected: identity:   %zu\n",
                stats.alignments_filtered_by_identity);
    std::printf("    rejected: length:     %zu\n",
                stats.alignments_filtered_by_length);
    std::printf("    avg identity:         %.3f\n", stats.avg_identity);
    std::printf("  Phase breakdown (sum across threads):\n");
    std::printf("    seeding:     %.2f s\n", stats.seeding_time_ms / 1000.0f);
    std::printf("    chaining:    %.2f s\n", stats.chaining_time_ms / 1000.0f);
    std::printf("    extension:   %.2f s\n", stats.extension_time_ms / 1000.0f);
    std::printf("  Align time:     %.2f s\n", align_time_ms / 1000.0f);
    std::printf("  Total time:     %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Throughput:     %.1f reads/s\n",
                total_reads / (align_time_ms / 1000.0f));
    std::printf("  Output:         %s\n", args.output.c_str());
    if (!args.parquet_output.empty()) {
        std::printf("  Parquet:        %s\n", args.parquet_output.c_str());
    }
}

bool ShouldRunLlmDiagnostics(const AlignArgs& args, float mapping_rate) {
    bool llm_enabled = args.enable_llm && !args.classical_only;
    return llm_enabled && mapping_rate < args.llm_threshold;
}

bool IsLlmEnabledButSkipped(const AlignArgs& args, float mapping_rate) {
    bool llm_enabled = args.enable_llm && !args.classical_only;
    return llm_enabled && mapping_rate >= args.llm_threshold;
}

}  // namespace llmap::cli::align_internal
