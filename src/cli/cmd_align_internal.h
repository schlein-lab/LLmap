// LLmap — cmd_align internal header.
// Shared types and declarations for split cmd_align modules.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "classical/classical_pipeline.h"
#include "claude_agent/pipeline_agent.h"

namespace llmap::cli::align_internal {

struct AlignArgs {
    std::string reads;
    std::string reference;
    std::string output;
    std::string parquet_output;
    bool use_bam = false;
    bool use_sam = true;
    int kmer_size = 15;
    int window_size = 10;
    int min_chain = 30;
    float min_identity = 0.70f;
    int threads = 1;
    int max_chains = 10;
    bool verbose = false;
    bool help = false;

    bool enable_llm = false;
    std::string llm_api_key;
    float llm_threshold = 0.50f;
    std::string llm_work_dir;
};

void PrintAlignUsage();
bool ParseAlignArgs(int argc, char** argv, AlignArgs& args);

std::string GetApiKey(const AlignArgs& args);

std::optional<claude_agent::PipelineAgent> CreatePipelineAgent(
    const AlignArgs& args, bool verbose);

void RunLlmDiagnostics(
    claude_agent::PipelineAgent& agent,
    const AlignArgs& args,
    const std::vector<classical::ReadAlignmentResult>& results,
    std::size_t n_mapped, std::size_t n_unmapped,
    float avg_identity);

}  // namespace llmap::cli::align_internal
