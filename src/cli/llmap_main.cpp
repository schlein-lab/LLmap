// LLmap CLI entry point.
//
// Full CLI surface (index, align, cnv-from-coverage, allpair, assign-paralog,
// sc-paralog-matrix, bench) implemented per LLmap_SPEC.md §6.
// Each command is in its own cmd_*.cpp file.

#include <cstdio>
#include <cstring>

#include "cli/commands.h"
#include "core/version.h"

namespace {

void print_banner() {
    std::printf(
        "LLmap %.*s\n"
        "Lossless. LLM-augmented. Wave-particle.\n"
        "Where reads see each other first.\n"
        "%.*s\n\n",
        static_cast<int>(llmap::kVersion.size()), llmap::kVersion.data(),
        static_cast<int>(llmap::kHomepageUrl.size()), llmap::kHomepageUrl.data()
    );
}

void print_version_full() {
    std::printf("llmap %.*s\n",
        static_cast<int>(llmap::kVersion.size()), llmap::kVersion.data());
    std::printf("  commit:   %.*s\n",
        static_cast<int>(llmap::kGitCommit.size()), llmap::kGitCommit.data());
    std::printf("  built:    %.*s\n",
        static_cast<int>(llmap::kBuildDate.size()), llmap::kBuildDate.data());
    std::printf("  type:     %.*s\n",
        static_cast<int>(llmap::kBuildType.size()), llmap::kBuildType.data());
    std::printf("  compiler: %.*s %.*s\n",
        static_cast<int>(llmap::kCompilerId.size()), llmap::kCompilerId.data(),
        static_cast<int>(llmap::kCompilerVersion.size()), llmap::kCompilerVersion.data());
    std::printf("  features: %.*s\n",
        static_cast<int>(llmap::kFeatures.size()), llmap::kFeatures.data());
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
        "  sc-qc-report      Single-cell QC report with filtering (Phase 9)\n"
        "  bench             Cross-tool benchmark vs minimap2 etc. (Phase 8)\n"
        "  validate          Synthetic + real-data validation harness (Phase 5)\n"
        "  check             V1.0 readiness check (Phase 10)\n"
        "  annotate-ref      Compute region annotations for a reference (Phase 12)\n"
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

    if (std::strcmp(argv[1], "--version") == 0 || std::strcmp(argv[1], "-V") == 0) {
        print_version_full();
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

    if (std::strcmp(argv[1], "index") == 0) {
        return llmap::cli::run_index(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "align") == 0) {
        return llmap::cli::run_align(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "sc-paralog-matrix") == 0) {
        return llmap::cli::run_sc_paralog_matrix(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "sc-qc-report") == 0) {
        return llmap::cli::run_sc_qc_report(argc - 2, argv + 2);
    }

    if (std::strcmp(argv[1], "check") == 0) {
        return llmap::cli::run_check(argc, argv);
    }

    if (std::strcmp(argv[1], "annotate-ref") == 0) {
        return llmap::cli::run_annotate_ref(argc - 2, argv + 2);
    }

    std::fprintf(stderr,
        "llmap: command '%s' not yet implemented.\n"
        "Use 'llmap --help' for available commands.\n",
        argv[1]);
    return 64;  // EX_USAGE
}
