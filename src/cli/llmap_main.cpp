// LLmap CLI entry point.
//
// Full CLI surface (index, align, cnv-from-coverage, allpair, assign-paralog,
// sc-paralog-matrix, bench) implemented per LLmap_SPEC.md §6.
// Each command is in its own cmd_*.cpp file.

#include <cstdio>
#include <cstring>

#include "cli/commands.h"

namespace {

constexpr const char* kVersion = "0.1.0-phase6";

void print_banner() {
    std::puts(
        "LLmap " "0.1.0-phase6" "\n"
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
        "  generate-synth    Generate synthetic IGH locus data\n"
        "  validate-real     Real reference validation (hg38 IGH locus)\n"
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
        return llmap::cli::run_allpair(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "generate-synth") == 0) {
        return llmap::cli::run_generate_synth(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "validate-real") == 0) {
        return llmap::cli::run_validate_real(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "align") == 0) {
        return llmap::cli::run_align(argc - 2, argv + 2);
    }

    std::fprintf(stderr,
        "llmap: command '%s' not yet implemented.\n"
        "Use 'llmap --help' for available commands.\n",
        argv[1]);
    return 64;  // EX_USAGE
}
