// LLmap — `llmap index` CLI command.
//
// Builds a minimizer index from a reference FASTA file. The index can be
// loaded by `llmap align` for faster alignment initialization.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "classical/minimizer_index.h"
#include "io/fasta_reader.h"

namespace llmap::cli {

namespace {

struct IndexArgs {
    std::string reference;
    std::string output;
    int kmer_size = 19;
    int window_size = 19;
    std::size_t max_occ = 500;
    bool verbose = false;
    bool help = false;
};

void PrintIndexUsage() {
    std::puts(
        "Usage: llmap index [options]\n"
        "\n"
        "Build a minimizer index from a reference FASTA file.\n"
        "\n"
        "Required:\n"
        "  -r, --reference FILE    Reference genome (FASTA)\n"
        "  -o, --output FILE       Output index file (.llmi)\n"
        "\n"
        "Index parameters:\n"
        "  -k, --kmer INT          Minimizer k-mer size [19]\n"
        "  -w, --window INT        Minimizer window [19]\n"
        "  --max-occ INT           Skip minimizers with count > N [500]\n"
        "\n"
        "Other:\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Example:\n"
        "  llmap index -r hg38.fa -o hg38.llmi\n"
        "  llmap index -r ref.fasta -o ref.llmi -k 15 -w 10\n"
    );
}

bool ParseIndexArgs(int argc, char** argv, IndexArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if ((arg == "-r" || arg == "--reference") && i + 1 < argc) {
            args.reference = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            args.output = argv[++i];
        } else if ((arg == "-k" || arg == "--kmer") && i + 1 < argc) {
            args.kmer_size = std::stoi(argv[++i]);
        } else if ((arg == "-w" || arg == "--window") && i + 1 < argc) {
            args.window_size = std::stoi(argv[++i]);
        } else if (arg == "--max-occ" && i + 1 < argc) {
            args.max_occ = static_cast<std::size_t>(std::stoi(argv[++i]));
        } else if (arg == "-v" || arg == "--verbose") {
            args.verbose = true;
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_index(int argc, char** argv) {
    IndexArgs args;
    if (!ParseIndexArgs(argc, argv, args)) {
        PrintIndexUsage();
        return 1;
    }

    if (args.help) {
        PrintIndexUsage();
        return 0;
    }

    if (args.reference.empty()) {
        std::fprintf(stderr, "Error: --reference is required\n");
        PrintIndexUsage();
        return 1;
    }

    if (args.output.empty()) {
        std::fprintf(stderr, "Error: --output is required\n");
        PrintIndexUsage();
        return 1;
    }

    if (!std::filesystem::exists(args.reference)) {
        std::fprintf(stderr, "Error: reference file not found: %s\n",
                     args.reference.c_str());
        return 1;
    }

    if (args.kmer_size < 5 || args.kmer_size > 31) {
        std::fprintf(stderr, "Error: k-mer size must be between 5 and 31\n");
        return 1;
    }

    if (args.window_size < 1 || args.window_size > 255) {
        std::fprintf(stderr, "Error: window size must be between 1 and 255\n");
        return 1;
    }

    auto start_time = std::chrono::steady_clock::now();

    if (args.verbose) {
        std::fprintf(stderr, "Loading reference: %s\n", args.reference.c_str());
    }

    io::FastaReader reader(args.reference);
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    std::uint64_t total_length = 0;

    while (reader.HasMore()) {
        auto record = reader.Next();
        if (record.IsValid()) {
            total_length += record.sequence.size();
            names.push_back(std::move(record.name));
            sequences.push_back(std::move(record.sequence));
        }
    }

    if (names.empty()) {
        std::fprintf(stderr, "Error: no sequences in reference\n");
        return 1;
    }

    auto load_time = std::chrono::steady_clock::now();
    float load_time_ms = std::chrono::duration<float, std::milli>(
        load_time - start_time).count();

    if (args.verbose) {
        std::fprintf(stderr, "Loaded %zu sequences (%.2f MB) in %.2f s\n",
                     names.size(),
                     static_cast<float>(total_length) / (1024.0f * 1024.0f),
                     load_time_ms / 1000.0f);
        std::fprintf(stderr, "Building index (k=%d, w=%d, max_occ=%zu)...\n",
                     args.kmer_size, args.window_size, args.max_occ);
    }

    classical::MinimizerConfig config;
    config.k = static_cast<std::uint8_t>(args.kmer_size);
    config.w = static_cast<std::uint8_t>(args.window_size);
    config.max_occ = args.max_occ;

    classical::MinimizerIndex::Builder builder(config);
    for (std::size_t i = 0; i < names.size(); ++i) {
        builder.AddSequence(names[i], sequences[i]);
    }

    auto index = builder.Build();

    auto build_time = std::chrono::steady_clock::now();
    float build_time_ms = std::chrono::duration<float, std::milli>(
        build_time - load_time).count();

    if (args.verbose) {
        std::fprintf(stderr, "Index built in %.2f s\n", build_time_ms / 1000.0f);
        std::fprintf(stderr, "  Minimizers: %zu\n", index->Size());
        const auto& stats = index->GetStats();
        std::fprintf(stderr, "  Total k-mers: %zu\n", stats.total_kmers);
        std::fprintf(stderr, "  Unique minimizers: %zu\n", stats.unique_minimizers);
        std::fprintf(stderr, "  Suppressed (high-occ): %zu\n", stats.suppressed_high_occ);
        std::fprintf(stderr, "  Avg spacing: %.2f\n", stats.avg_minimizer_spacing);
    }

    if (args.verbose) {
        std::fprintf(stderr, "Saving index: %s\n", args.output.c_str());
    }

    if (!index->Save(args.output)) {
        std::fprintf(stderr, "Error: failed to save index to %s\n",
                     args.output.c_str());
        return 1;
    }

    auto save_time = std::chrono::steady_clock::now();
    float save_time_ms = std::chrono::duration<float, std::milli>(
        save_time - build_time).count();

    auto file_size = std::filesystem::file_size(args.output);
    float total_time_ms = std::chrono::duration<float, std::milli>(
        save_time - start_time).count();

    std::printf("Index built successfully:\n");
    std::printf("  Reference:     %s\n", args.reference.c_str());
    std::printf("  Sequences:     %zu\n", names.size());
    std::printf("  Total length:  %.2f MB\n",
                static_cast<float>(total_length) / (1024.0f * 1024.0f));
    std::printf("  k-mer size:    %d\n", args.kmer_size);
    std::printf("  Window size:   %d\n", args.window_size);
    std::printf("  Minimizers:    %zu\n", index->Size());
    std::printf("  Index size:    %.2f MB\n",
                static_cast<float>(file_size) / (1024.0f * 1024.0f));
    std::printf("  Build time:    %.2f s\n", build_time_ms / 1000.0f);
    std::printf("  Save time:     %.2f s\n", save_time_ms / 1000.0f);
    std::printf("  Total time:    %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Output:        %s\n", args.output.c_str());

    return 0;
}

}  // namespace llmap::cli
