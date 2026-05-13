// LLmap — `llmap align` CLI command.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "classical/classical_pipeline.h"
#include "claude_agent/pipeline_agent.h"
#include "core/alignment_record.h"
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"

namespace llmap::cli {

namespace {

void print_align_usage() {
    std::puts(
        "Usage: llmap align [options]\n"
        "\n"
        "Align reads to a reference genome.\n"
        "\n"
        "Required:\n"
        "  -r, --reads FILE        Input reads (FASTQ)\n"
        "  -x, --reference FILE    Reference genome (FASTA)\n"
        "  -o, --output FILE       Output alignment file (SAM/BAM)\n"
        "\n"
        "Output format:\n"
        "  --bam                   Output BAM format (requires htslib)\n"
        "  --sam                   Output SAM format (default)\n"
        "  --parquet FILE          Also output probabilistic Parquet\n"
        "\n"
        "Alignment parameters:\n"
        "  -k, --kmer INT          Minimizer k-mer size [15]\n"
        "  -w, --window INT        Minimizer window [10]\n"
        "  --min-chain INT         Minimum chain score [30]\n"
        "  --min-identity FLOAT    Minimum alignment identity [0.70]\n"
        "\n"
        "Performance:\n"
        "  -t, --threads INT       Number of threads [1]\n"
        "  --max-chains INT        Max chains to extend per read [10]\n"
        "\n"
        "LLM-assisted diagnostics:\n"
        "  --llm                   Enable Claude LLM for alignment diagnostics\n"
        "  --llm-api-key KEY       Anthropic API key (or set ANTHROPIC_API_KEY)\n"
        "  --llm-threshold FLOAT   Mapping rate threshold to trigger diagnostics [0.50]\n"
        "  --llm-work-dir DIR      Working directory for LLM artifacts\n"
        "\n"
        "Other:\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Example:\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.bam --bam --parquet out.parquet\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam --llm\n"
    );
}

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

bool parse_align_args(int argc, char** argv, AlignArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if ((arg == "-r" || arg == "--reads") && i + 1 < argc) {
            args.reads = argv[++i];
        } else if ((arg == "-x" || arg == "--reference") && i + 1 < argc) {
            args.reference = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            args.output = argv[++i];
        } else if (arg == "--parquet" && i + 1 < argc) {
            args.parquet_output = argv[++i];
        } else if (arg == "--bam") {
            args.use_bam = true;
            args.use_sam = false;
        } else if (arg == "--sam") {
            args.use_sam = true;
            args.use_bam = false;
        } else if ((arg == "-k" || arg == "--kmer") && i + 1 < argc) {
            args.kmer_size = std::stoi(argv[++i]);
        } else if ((arg == "-w" || arg == "--window") && i + 1 < argc) {
            args.window_size = std::stoi(argv[++i]);
        } else if (arg == "--min-chain" && i + 1 < argc) {
            args.min_chain = std::stoi(argv[++i]);
        } else if (arg == "--min-identity" && i + 1 < argc) {
            args.min_identity = std::stof(argv[++i]);
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            args.threads = std::stoi(argv[++i]);
        } else if (arg == "--max-chains" && i + 1 < argc) {
            args.max_chains = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            args.verbose = true;
        } else if (arg == "--llm") {
            args.enable_llm = true;
        } else if (arg == "--llm-api-key" && i + 1 < argc) {
            args.llm_api_key = argv[++i];
        } else if (arg == "--llm-threshold" && i + 1 < argc) {
            args.llm_threshold = std::stof(argv[++i]);
        } else if (arg == "--llm-work-dir" && i + 1 < argc) {
            args.llm_work_dir = argv[++i];
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

std::string CigarToString(const std::vector<classical::CigarElement>& cigar) {
    std::ostringstream oss;
    for (const auto& elem : cigar) {
        oss << elem.ToString();
    }
    return oss.str();
}

AlignmentRecord ConvertAlignment(const classical::ClassicalAlignment& aln,
                                  std::uint32_t read_len) {
    AlignmentHit hit;
    hit.target_id = aln.ref_name;
    hit.start = static_cast<std::uint64_t>(aln.ref_start);
    hit.end = static_cast<std::uint64_t>(aln.ref_end);
    hit.cigar = CigarString{CigarToString(aln.cigar)};
    hit.score = aln.score;
    hit.nm = static_cast<std::uint32_t>((1.0f - aln.identity) * aln.AlignedBases());

    return make_mapped(aln.query_name, read_len, std::move(hit));
}

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

}  // namespace

int run_align(int argc, char** argv) {
    AlignArgs args;
    if (!parse_align_args(argc, argv, args)) {
        print_align_usage();
        return 1;
    }

    if (args.help) {
        print_align_usage();
        return 0;
    }

    if (args.reads.empty()) {
        std::fprintf(stderr, "Error: --reads is required\n");
        print_align_usage();
        return 1;
    }

    if (args.reference.empty()) {
        std::fprintf(stderr, "Error: --reference is required\n");
        print_align_usage();
        return 1;
    }

    if (args.output.empty()) {
        std::fprintf(stderr, "Error: --output is required\n");
        print_align_usage();
        return 1;
    }

    if (!std::filesystem::exists(args.reads)) {
        std::fprintf(stderr, "Error: reads file not found: %s\n", args.reads.c_str());
        return 1;
    }

    if (!std::filesystem::exists(args.reference)) {
        std::fprintf(stderr, "Error: reference file not found: %s\n", args.reference.c_str());
        return 1;
    }

    auto start_time = std::chrono::steady_clock::now();

    // Load reference sequences
    if (args.verbose) {
        std::fprintf(stderr, "Loading reference: %s\n", args.reference.c_str());
    }

    io::FastaReader ref_reader(args.reference);
    std::vector<std::string> ref_names;
    std::vector<std::string> ref_seqs;

    while (ref_reader.HasMore()) {
        auto record = ref_reader.Next();
        if (record.IsValid()) {
            ref_names.push_back(record.name);
            ref_seqs.push_back(record.sequence);
        }
    }

    if (ref_names.empty()) {
        std::fprintf(stderr, "Error: no sequences in reference\n");
        return 1;
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loaded %zu reference sequences\n", ref_names.size());
    }

    // Build minimizer index
    if (args.verbose) {
        std::fprintf(stderr, "Building minimizer index (k=%d, w=%d)\n",
                     args.kmer_size, args.window_size);
    }

    classical::MinimizerConfig mini_cfg;
    mini_cfg.k = args.kmer_size;
    mini_cfg.w = args.window_size;

    classical::MinimizerIndex::Builder idx_builder(mini_cfg);
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        idx_builder.AddSequence(ref_names[i], ref_seqs[i]);
    }
    auto index = idx_builder.Build();

    if (args.verbose) {
        std::fprintf(stderr, "Index built: %zu minimizers\n", index->Size());
    }

    // Configure pipeline
    classical::ClassicalPipelineConfig pipe_cfg;
    pipe_cfg.minimizer_config = mini_cfg;
    pipe_cfg.chain_config.min_chain_score = args.min_chain;
    pipe_cfg.extension_config.match_score = 2;
    pipe_cfg.extension_config.mismatch_penalty = 4;
    pipe_cfg.extension_config.gap_open = 4;
    pipe_cfg.extension_config.gap_extend = 2;
    pipe_cfg.min_identity = args.min_identity;
    pipe_cfg.max_chains_to_extend = args.max_chains;

    classical::ClassicalPipeline pipeline(pipe_cfg);
    pipeline.SetIndex(std::move(index));

    // Load reads
    if (args.verbose) {
        std::fprintf(stderr, "Loading reads: %s\n", args.reads.c_str());
    }

    auto read_reader = io::FastqReader::Open(args.reads);
    if (!read_reader) {
        std::fprintf(stderr, "Error: failed to open reads file: %s\n", args.reads.c_str());
        return 1;
    }

    std::vector<std::string> read_names;
    std::vector<std::string> read_seqs;
    std::vector<std::uint32_t> read_lens;

    while (read_reader->HasMore()) {
        auto record = read_reader->Next();
        if (record && record->IsValid()) {
            read_names.push_back(record->id);
            read_seqs.push_back(record->sequence);
            read_lens.push_back(static_cast<std::uint32_t>(record->sequence.size()));
        }
    }

    if (read_names.empty()) {
        std::fprintf(stderr, "Error: no reads in input\n");
        return 1;
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loaded %zu reads\n", read_names.size());
    }

    // Align reads
    if (args.verbose) {
        std::fprintf(stderr, "Aligning...\n");
    }

    auto align_start = std::chrono::steady_clock::now();
    auto results = pipeline.AlignReads(read_names, read_seqs);
    auto align_end = std::chrono::steady_clock::now();

    float align_time_ms = std::chrono::duration<float, std::milli>(
        align_end - align_start).count();

    // Convert to AlignmentRecords
    std::vector<AlignmentRecord> records;
    records.reserve(results.size());

    std::size_t n_mapped = 0;
    std::size_t n_unmapped = 0;

    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];

        if (result.HasAlignment()) {
            const auto* primary = result.PrimaryAlignment();
            if (primary) {
                records.push_back(ConvertAlignment(*primary, read_lens[i]));
                n_mapped++;
            }
        } else {
            records.push_back(make_unmapped(
                result.query_name, read_lens[i], RejectionReason::NoSeeds));
            n_unmapped++;
        }
    }

    // Write SAM/BAM output
    if (args.verbose) {
        std::fprintf(stderr, "Writing output: %s\n", args.output.c_str());
    }

    std::vector<output::ReferenceSequence> ref_info;
    ref_info.reserve(ref_names.size());
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        ref_info.push_back({ref_names[i], ref_seqs[i].size()});
    }

    output::BamWriterConfig bam_cfg;
    bam_cfg.format = args.use_bam ? output::BamOutputFormat::BAM
                                   : output::BamOutputFormat::SAM;
    bam_cfg.include_wavecollapse_tags = false;
    bam_cfg.include_paralog_tags = false;

    auto bam_writer = output::BamWriter::Create(args.output, ref_info, bam_cfg);
    if (!bam_writer) {
        std::fprintf(stderr, "Error: failed to create output file: %s\n",
                     args.output.c_str());
        return 1;
    }

    if (!bam_writer->WriteBatch(records)) {
        std::fprintf(stderr, "Error: failed to write alignments: %s\n",
                     bam_writer->LastError().c_str());
        return 1;
    }

    bam_writer->Close();

    // Optionally write Parquet output
    if (!args.parquet_output.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Writing Parquet: %s\n", args.parquet_output.c_str());
        }

        output::ParquetWriterConfig parquet_cfg;
        parquet_cfg.min_probability = 0.0f;

        auto parquet_writer = output::ParquetWriter::Create(
            args.parquet_output, parquet_cfg);
        if (!parquet_writer) {
            std::fprintf(stderr, "Warning: failed to create Parquet output\n");
        } else {
            parquet_writer->WriteBatch(records);
            parquet_writer->Close();
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    // Summary
    float mapping_rate = static_cast<float>(n_mapped) /
        static_cast<float>(read_names.size());

    std::printf("Alignment complete:\n");
    std::printf("  Input reads:    %zu\n", read_names.size());
    std::printf("  Mapped:         %zu (%.1f%%)\n",
                n_mapped, 100.0f * mapping_rate);
    std::printf("  Unmapped:       %zu\n", n_unmapped);
    std::printf("  Align time:     %.2f s\n", align_time_ms / 1000.0f);
    std::printf("  Total time:     %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Throughput:     %.1f reads/s\n",
                read_names.size() / (align_time_ms / 1000.0f));
    std::printf("  Output:         %s\n", args.output.c_str());
    if (!args.parquet_output.empty()) {
        std::printf("  Parquet:        %s\n", args.parquet_output.c_str());
    }

    if (args.enable_llm && mapping_rate < args.llm_threshold) {
        if (args.verbose) {
            std::fprintf(stderr, "\nMapping rate %.1f%% below threshold %.1f%%, "
                         "running LLM diagnostics...\n",
                         100.0f * mapping_rate, 100.0f * args.llm_threshold);
        }

        auto agent = CreatePipelineAgent(args, args.verbose);
        if (agent) {
            float avg_identity = pipeline.Stats().avg_identity;
            RunLlmDiagnostics(*agent, args, results, n_mapped, n_unmapped, avg_identity);
        }
    } else if (args.enable_llm && args.verbose) {
        std::fprintf(stderr, "Mapping rate %.1f%% meets threshold, "
                     "skipping LLM diagnostics\n", 100.0f * mapping_rate);
    }

    return 0;
}

}  // namespace llmap::cli
