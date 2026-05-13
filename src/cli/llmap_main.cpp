// LLmap CLI entry point — Phase 0 stub.
// Full CLI surface (index, align, cnv-from-coverage, allpair, assign-paralog,
// sc-paralog-matrix, bench) implemented in later phases per LLmap_SPEC.md §6.

#include <cstdio>
#include <cstring>
#include <string>

namespace {
constexpr const char* kVersion = "0.0.1-bootstrap";

void print_banner() {
    std::puts(
        "LLmap " "0.0.1-bootstrap" "\n"
        "Lossless. LLM-augmented. Wave-particle.\n"
        "Where reads see each other first.\n"
        "https://losslessmap.com\n"
    );
}

void print_usage() {
    std::puts(
        "Usage: llmap <command> [options]\n"
        "\n"
        "Commands (V1.0 roadmap):\n"
        "  index             Build LLmap reference index (Phase 3+)\n"
        "  align             Align reads using WaveCollapse (Phase 3+)\n"
        "  allpair           Stage 1 Self-Interference standalone (Phase 2)\n"
        "  cnv-from-coverage CNV inference from probabilistic coverage (Phase 6)\n"
        "  assign-paralog    PSV-based paralog assignment (Phase 9)\n"
        "  sc-paralog-matrix Single-cell cell × paralog AnnData (Phase 9)\n"
        "  bench             Cross-tool benchmark vs minimap2 etc. (Phase 8)\n"
        "  validate          Synthetic + real-data validation harness (Phase 5)\n"
        "  --version         Show version\n"
        "  --help            Show this message\n"
        "\n"
        "Status: V1.0 in active autonomous build. See STATE.md.\n"
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
    std::fprintf(stderr,
        "llmap: command '%s' not yet implemented (Phase 0 bootstrap).\n"
        "See STATE.md for build progress.\n",
        argv[1]);
    return 64;  // EX_USAGE
}
