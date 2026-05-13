// LLmap — cmd_align LLM diagnostics.

#include "cli/cmd_align_internal.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>

namespace llmap::cli::align_internal {

std::string GetApiKey(const AlignArgs& args) {
    if (!args.llm_api_key.empty()) {
        return args.llm_api_key;
    }
    const char* env_key = std::getenv("ANTHROPIC_API_KEY");
    return env_key ? std::string(env_key) : "";
}

std::optional<claude_agent::PipelineAgent> CreatePipelineAgent(
    const AlignArgs& args, bool verbose) {

    std::string api_key = GetApiKey(args);
    if (api_key.empty()) {
        if (verbose) {
            std::fprintf(stderr, "Warning: --llm enabled but no API key provided. "
                         "Set ANTHROPIC_API_KEY or use --llm-api-key\n");
        }
        return std::nullopt;
    }

    claude_agent::PipelineAgentConfig config;
    config.agent.api_key = api_key;
    config.agent.model = "claude-sonnet-4-20250514";
    config.agent.max_tokens = 4096;
    config.enable_diagnostics = true;
    config.enable_reporter = true;
    config.enable_cuda_codegen = false;

    if (!args.llm_work_dir.empty()) {
        config.work_dir = args.llm_work_dir;
    } else {
        config.work_dir = std::filesystem::temp_directory_path() / "llmap_llm";
    }

    std::filesystem::create_directories(config.work_dir);

    return claude_agent::PipelineAgent(std::move(config));
}

void RunLlmDiagnostics(
    claude_agent::PipelineAgent& agent,
    const AlignArgs& args,
    const std::vector<classical::ReadAlignmentResult>& results,
    std::size_t n_mapped, std::size_t n_unmapped,
    float avg_identity) {

    float mapping_rate = static_cast<float>(n_mapped) /
        static_cast<float>(n_mapped + n_unmapped);

    std::printf("\n--- LLM Diagnostics ---\n");
    std::printf("Session ID: %s\n", agent.SessionId().c_str());

    claude_agent::StallMetrics stall;
    stall.type = claude_agent::StallType::NoProgress;
    stall.reads_affected = n_unmapped;

    std::vector<std::pair<std::string, double>> bucket_probs;
    for (std::size_t i = 0; i < std::min(results.size(), std::size_t{100}); ++i) {
        const auto& r = results[i];
        double prob = r.HasAlignment() ? 1.0 : 0.0;
        bucket_probs.emplace_back(r.query_name, prob);
    }

    auto work_dir = agent.GetConfig().work_dir;
    auto wave_state = claude_agent::WriteWaveStateJson(
        work_dir, 1, 1.0 - mapping_rate, bucket_probs);

    if (wave_state.empty()) {
        std::fprintf(stderr, "Warning: failed to write wave state for diagnostics\n");
        return;
    }

    claude_agent::DiagnosticContext ctx;
    ctx.stall = stall;
    ctx.wave_state_path = wave_state;
    ctx.output_dir = work_dir / "diagnostics";

    std::printf("Analyzing alignment issues...\n");
    auto resolution = agent.DiagnoseAndResolve(ctx);

    std::printf("\nDiagnostic Report:\n");
    std::printf("  Stall Pattern:    %s\n", resolution.report.stall_pattern.c_str());
    std::printf("  Root Cause:       %s\n", resolution.report.root_cause.c_str());
    std::printf("  Resolution:       %s\n", resolution.report.resolution.c_str());
    std::printf("  Latency:          %ld ms\n", resolution.latency.count());

    if (resolution.report.custom_kernel_path) {
        std::printf("  Custom Kernel:    %s\n", resolution.report.custom_kernel_path->c_str());
    }
    if (resolution.report.kernel_hot_loaded) {
        std::printf("  Kernel Loaded:    yes\n");
    }

    std::printf("-----------------------\n");
}

}  // namespace llmap::cli::align_internal
