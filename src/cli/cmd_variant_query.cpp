// LLmap -- `llmap variant-query` CLI command (Layer 4 stub).
//
// Spot-check tool for variant prior tables. Dumps prior records for a
// region to stdout, used to verify ingest correctness and to inspect what
// the runtime will see at a given locus.
//
// Usage:
//   llmap variant-query --priors dbsnp.priors
//                       --chr chr14
//                       --pos 105580000
//                       [--window 5000]
//
// Multiple --priors paths can be passed (comma-separated) to inspect the
// aggregated view exactly as `llmap align` will see it.
//
// Not wired into llmap_main.cpp yet — this is a scaffold for the upcoming
// variant-priors feature. Run-functions print "not yet implemented".
//
// TODO(layer4):
//   - mmap-based reader for llmap_priors_v1 format
//   - per-chromosome offset sidecar (.idx) reader
//   - multi-file aggregation per schema.json:aggregation_rules
//   - TSV / JSON output modes

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

struct VariantQueryArgs {
    std::vector<std::string> priors_paths;  // comma-split list
    std::string chr;
    uint64_t pos = 0;
    uint64_t window = 0;
    std::string output_mode = "tsv";        // tsv | json
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap variant-query --priors PRIORS[,PRIORS...] \\\n"
        "                           --chr CHR --pos POS [--window WIN]\n"
        "                           [--output tsv|json]\n"
        "\n"
        "Dumps prior records for a region.  WIN is half-window in bp around\n"
        "POS (default 0 = single bucket).  Multiple --priors files are\n"
        "aggregated according to knowledge/variants/schema.json.\n");
}

void SplitCsv(const std::string& s, std::vector<std::string>& out) {
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == ',') {
            if (i > start) out.emplace_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
}

bool ParseArgs(int argc, char** argv, VariantQueryArgs& a) {
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--priors" && i + 1 < argc) SplitCsv(argv[++i], a.priors_paths);
        else if (arg == "--chr"    && i + 1 < argc) a.chr = argv[++i];
        else if (arg == "--pos"    && i + 1 < argc) a.pos = std::stoull(argv[++i]);
        else if (arg == "--window" && i + 1 < argc) a.window = std::stoull(argv[++i]);
        else if (arg == "--output" && i + 1 < argc) a.output_mode = argv[++i];
        else if (arg == "--help" || arg == "-h")    { a.help = true; return true; }
        else {
            std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_variant_query(int argc, char** argv) {
    VariantQueryArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 1; }
    if (args.help) { PrintUsage(); return 0; }

    if (args.priors_paths.empty() || args.chr.empty() || args.pos == 0) {
        std::fprintf(stderr, "variant-query: --priors, --chr, --pos are required.\n");
        PrintUsage();
        return 1;
    }

    std::fprintf(stderr,
        "[variant-query] not yet implemented.\n"
        "  priors files  = %zu\n",
        args.priors_paths.size());
    for (const auto& p : args.priors_paths) {
        std::fprintf(stderr, "    - %s\n", p.c_str());
    }
    std::fprintf(stderr,
        "  chr           = %s\n"
        "  pos           = %llu\n"
        "  window        = %llu\n"
        "  output_mode   = %s\n"
        "\n"
        "Scaffold only; see TODO(layer4) markers in src/cli/cmd_variant_query.cpp\n",
        args.chr.c_str(),
        static_cast<unsigned long long>(args.pos),
        static_cast<unsigned long long>(args.window),
        args.output_mode.c_str());
    return 2;  // ENOSYS-style: command recognised, not yet implemented.
}

}  // namespace llmap::cli
