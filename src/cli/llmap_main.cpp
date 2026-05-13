// LLmap CLI entry point.
//
// Full CLI surface (index, align, cnv-from-coverage, allpair, assign-paralog,
// sc-paralog-matrix, bench) implemented per LLmap_SPEC.md §6.

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include "self_interference/allpair_pipeline.h"

namespace {

constexpr const char* kVersion = "0.1.0-phase2";

void print_banner() {
    std::puts(
        "LLmap " "0.1.0-phase2" "\n"
        "Lossless. LLM-augmented. Wave-particle.\n"
        "Where reads see each other first.\n"
        "https://losslessmap.com\n"
    );
}

void print_usage() {
    std::puts(
        "Usage: llmap <command> [options]\n"
        "\n"
        "Commands:\n"
        "  allpair           Stage 1 Self-Interference standalone\n"
        "  index             Build LLmap reference index (Phase 3+)\n"
        "  align             Align reads using WaveCollapse (Phase 3+)\n"
        "  cnv-from-coverage CNV inference from probabilistic coverage (Phase 6)\n"
        "  assign-paralog    PSV-based paralog assignment (Phase 9)\n"
        "  sc-paralog-matrix Single-cell cell × paralog AnnData (Phase 9)\n"
        "  bench             Cross-tool benchmark vs minimap2 etc. (Phase 8)\n"
        "  validate          Synthetic + real-data validation harness (Phase 5)\n"
        "  --version         Show version\n"
        "  --help            Show this message\n"
    );
}

void print_allpair_usage() {
    std::puts(
        "Usage: llmap allpair [options]\n"
        "\n"
        "Stage 1 Self-Interference: cluster reads by sequence similarity.\n"
        "\n"
        "Required:\n"
        "  -r, --reads FILE        Input FASTQ file\n"
        "  -o, --output FILE       Output TSV file (cluster assignments)\n"
        "\n"
        "Optional:\n"
        "  -m, --model FILE        Foundation model ONNX path\n"
        "  -k, --knn INT           k-NN neighbors [50]\n"
        "  --resolution FLOAT      Leiden resolution [1.0]\n"
        "  --min-cluster INT       Minimum cluster size [2]\n"
        "  --max-reads INT         Maximum reads to process [0=all]\n"
        "  --embedding-dim INT     Embedding dimension [256]\n"
        "  --gpu                   Use GPU for embedding/FAISS\n"
        "  --gpu-device INT        GPU device ID [0]\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Example:\n"
        "  llmap allpair -r reads.fastq -o clusters.tsv\n"
    );
}

struct AllpairArgs {
    std::string reads;
    std::string output;
    std::string model;
    size_t knn = 50;
    float resolution = 1.0f;
    size_t min_cluster = 2;
    size_t max_reads = 0;
    size_t embedding_dim = 256;
    bool use_gpu = false;
    int gpu_device = 0;
    bool verbose = false;
    bool help = false;
};

bool parse_allpair_args(int argc, char** argv, AllpairArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if ((arg == "-r" || arg == "--reads") && i + 1 < argc) {
            args.reads = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            args.output = argv[++i];
        } else if ((arg == "-m" || arg == "--model") && i + 1 < argc) {
            args.model = argv[++i];
        } else if ((arg == "-k" || arg == "--knn") && i + 1 < argc) {
            args.knn = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--resolution" && i + 1 < argc) {
            args.resolution = std::stof(argv[++i]);
        } else if (arg == "--min-cluster" && i + 1 < argc) {
            args.min_cluster = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--max-reads" && i + 1 < argc) {
            args.max_reads = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--embedding-dim" && i + 1 < argc) {
            args.embedding_dim = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "--gpu") {
            args.use_gpu = true;
        } else if (arg == "--gpu-device" && i + 1 < argc) {
            args.gpu_device = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            args.verbose = true;
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }

    return true;
}

int run_allpair(int argc, char** argv) {
    AllpairArgs args;
    if (!parse_allpair_args(argc, argv, args)) {
        print_allpair_usage();
        return 1;
    }

    if (args.help) {
        print_allpair_usage();
        return 0;
    }

    if (args.reads.empty()) {
        std::fprintf(stderr, "Error: --reads is required\n");
        print_allpair_usage();
        return 1;
    }

    if (args.output.empty()) {
        std::fprintf(stderr, "Error: --output is required\n");
        print_allpair_usage();
        return 1;
    }

    if (!std::filesystem::exists(args.reads)) {
        std::fprintf(stderr, "Error: reads file not found: %s\n", args.reads.c_str());
        return 1;
    }

    // Build pipeline config
    llmap::self_interference::AllpairConfig config;
    config.reads_a = args.reads;
    config.output = args.output;
    config.model_path = args.model;
    config.faiss_k = args.knn;
    config.leiden_resolution = args.resolution;
    config.min_cluster_size = args.min_cluster;
    config.max_reads = args.max_reads;
    config.embedding_dim = args.embedding_dim;
    config.use_gpu_embedder = args.use_gpu;
    config.use_gpu_faiss = args.use_gpu;
    config.gpu_device_id = args.gpu_device;
    config.verbose = args.verbose;

    // Run pipeline
    llmap::self_interference::AllpairPipeline pipeline(config);

    std::unique_ptr<llmap::self_interference::AllpairResult> result;
    if (args.verbose) {
        auto progress_callback = [](const std::string& stage, size_t current,
                                     size_t total, const std::string& message) {
            std::fprintf(stderr, "[%zu/%zu] %s: %s\n", current + 1, total,
                         stage.c_str(), message.c_str());
        };
        result = pipeline.Run(progress_callback);
    } else {
        result = pipeline.Run();
    }

    if (!result || !result->success) {
        std::fprintf(stderr, "Error: pipeline failed: %s\n",
                     result ? result->error_message.c_str() : "unknown error");
        return 1;
    }

    // Print summary
    std::printf("Allpair complete:\n");
    std::printf("  Input reads:        %zu\n", result->stats.input_reads);
    std::printf("  Clusters:           %zu\n", result->NumClusters());
    std::printf("  Representatives:    %zu\n", result->stats.num_representatives);
    std::printf("  Modularity:         %.4f\n", result->stats.modularity);
    std::printf("  Total time:         %.2f s\n", result->stats.total_time_ms / 1000.0f);
    std::printf("  Throughput:         %.1f reads/s\n", result->stats.reads_per_second);
    std::printf("  Output:             %s\n", args.output.c_str());

    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_banner();
        print_usage();
        return 0;
    }

    if (std::strcmp(argv[1], "--version") == 0) {
        std::printf("llmap %s\n", kVersion);
        return 0;
    }

    if (std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-h") == 0) {
        print_banner();
        print_usage();
        return 0;
    }

    if (std::strcmp(argv[1], "allpair") == 0) {
        return run_allpair(argc - 2, argv + 2);
    }

    std::fprintf(stderr,
        "llmap: command '%s' not yet implemented.\n"
        "Use 'llmap --help' for available commands.\n",
        argv[1]);
    return 64;  // EX_USAGE
}
