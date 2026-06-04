// LLmap — `llmap junction-index build` CLI command.
//
// Build the persistent NAHR-pair k-mer cache once, then re-use it
// across many junction-hunt invocations via `--index-dir`.

#include "cli/commands.h"

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/nahr_pair_loader.h"
#include "junction_hunter/persistent_index_builder.h"
#include "junction_hunter/cascade_caller.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

struct Args {
    std::string pairs_tsv;
    std::string reference;
    std::string out_dir;
    std::string preset = "long";  // long | short
    std::vector<std::uint8_t> k_values;  // override preset if non-empty
    int max_pairs{0};
    bool verbose{false};
    bool help{false};
};

void Usage() {
    std::puts(
        "Usage: llmap junction-index build [options]\n"
        "\n"
        "Build the persistent multi-tier k-mer cache for the NAHR-pair panel.\n"
        "Memory-intensive: run on a bigmem node (~256-500 GB RAM).\n"
        "\n"
        "Required:\n"
        "  --pairs FILE          NAHR-pair TSV (linter output)\n"
        "  --reference FILE      Reference FASTA (GRCh38)\n"
        "  --out-dir DIR         Output directory for tier_k*.bin files\n"
        "\n"
        "Optional:\n"
        "  --preset NAME         long|short (default: long; HiFi/ONT k=11..125)\n"
        "  --k-values 11,17,...  Override preset, comma-separated k list\n"
        "  --max-pairs N         Limit panel for testing\n"
        "  -v, --verbose         Progress to stderr\n"
        "  -h, --help            Show this help\n"
    );
}

bool ParseArgs(int argc, char** argv, Args& a) {
    for (int i = 0; i < argc; ++i) {
        std::string x = argv[i];
        auto take = [&](const char* n) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "error: %s requires value\n", n); return nullptr; }
            return argv[++i];
        };
        if (x == "-h" || x == "--help") { a.help = true; return true; }
        else if (x == "-v" || x == "--verbose") a.verbose = true;
        else if (x == "--pairs")     { auto v = take("--pairs");     if (!v) return false; a.pairs_tsv = v; }
        else if (x == "--reference"|| x == "-r") { auto v = take("--reference"); if (!v) return false; a.reference = v; }
        else if (x == "--out-dir")   { auto v = take("--out-dir");   if (!v) return false; a.out_dir = v; }
        else if (x == "--preset")    { auto v = take("--preset");    if (!v) return false; a.preset = v; }
        else if (x == "--max-pairs") { auto v = take("--max-pairs"); if (!v) return false; a.max_pairs = std::stoi(v); }
        else if (x == "--k-values")  {
            auto v = take("--k-values"); if (!v) return false;
            std::string s = v;
            std::string cur;
            for (char c : s) {
                if (c == ',') {
                    if (!cur.empty()) { a.k_values.push_back(static_cast<std::uint8_t>(std::stoi(cur))); cur.clear(); }
                } else cur.push_back(c);
            }
            if (!cur.empty()) a.k_values.push_back(static_cast<std::uint8_t>(std::stoi(cur)));
        }
        else { std::fprintf(stderr, "error: unknown arg %s\n", x.c_str()); return false; }
    }
    if (a.pairs_tsv.empty() || a.reference.empty() || a.out_dir.empty()) {
        std::fprintf(stderr, "error: --pairs, --reference, --out-dir required\n");
        return false;
    }
    return true;
}

}  // namespace

int run_junction_index_build(int argc, char** argv) {
    Args args;
    if (!ParseArgs(argc, argv, args)) { Usage(); return 2; }
    if (args.help) { Usage(); return 0; }

    std::vector<junction_hunter::NahrPair> pairs;
    auto st = junction_hunter::LoadNahrPairsTsv(args.pairs_tsv, pairs);
    if (!st.ok) {
        std::fprintf(stderr, "error: pair-tsv load failed: %s\n", st.error.c_str());
        return 1;
    }
    if (args.max_pairs > 0 && static_cast<int>(pairs.size()) > args.max_pairs)
        pairs.resize(args.max_pairs);

    junction_hunter::BuildOptions opts;
    opts.out_dir        = args.out_dir;
    opts.panel_path     = args.pairs_tsv;
    opts.reference_path = args.reference;
    opts.verbose        = args.verbose;
    if (!args.k_values.empty()) {
        opts.k_values = args.k_values;
    } else if (args.preset == "short") {
        opts.k_values = junction_hunter::CascadeConfig::ShortReadPreset().k_values;
    } else {
        opts.k_values = junction_hunter::CascadeConfig::LongReadPreset().k_values;
    }

    junction_hunter::BuildStats stats;
    std::string err;
    if (!junction_hunter::BuildPersistentIndex(
            pairs, args.reference, opts, stats, err)) {
        std::fprintf(stderr, "error: build failed: %s\n", err.c_str());
        return 1;
    }

    std::fprintf(stderr, "[build] done in %.1fs — %zu pairs, %zu tiers\n",
                  stats.seconds_total, pairs.size(), opts.k_values.size());
    for (std::size_t t = 0; t < opts.k_values.size(); ++t) {
        std::fprintf(stderr, "  k=%-3u  %llu entries  %.2f GB\n",
                      static_cast<unsigned>(opts.k_values[t]),
                      static_cast<unsigned long long>(stats.entries_per_tier[t]),
                      stats.bytes_per_tier[t] / 1e9);
    }
    return 0;
}

int run_junction_index(int argc, char** argv) {
    if (argc < 1) {
        std::fprintf(stderr, "Usage: llmap junction-index <build> [options]\n");
        return 2;
    }
    std::string sub = argv[0];
    if (sub == "build") return run_junction_index_build(argc - 1, argv + 1);
    if (sub == "-h" || sub == "--help") {
        std::puts("Usage: llmap junction-index <subcommand> [options]\n"
                  "Subcommands:\n"
                  "  build    Build the persistent NAHR-pair k-mer cache\n");
        return 0;
    }
    std::fprintf(stderr, "error: unknown subcommand '%s'\n", sub.c_str());
    return 2;
}

}  // namespace llmap::cli
