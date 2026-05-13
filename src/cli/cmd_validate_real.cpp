// LLmap — `llmap validate-real` CLI command.

#include "cli/commands.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "validation/real_reference.h"

namespace llmap::cli {

namespace {

void print_validate_real_usage() {
    std::puts(
        "Usage: llmap validate-real [options]\n"
        "\n"
        "Validate LLmap against real reference data (hg38 IGH locus).\n"
        "\n"
        "Required:\n"
        "  --reference FILE    Reference FASTA (hg38 chr14 IGH locus)\n"
        "  --reads FILE        Input reads FASTQ\n"
        "\n"
        "Optional:\n"
        "  --baseline FILE     Minimap2 BAM for baseline comparison\n"
        "  --truth FILE        Ground truth BED file\n"
        "  --output-dir DIR    Output directory for results\n"
        "  --k INT             Minimizer k-mer size [15]\n"
        "  --w INT             Minimizer window [10]\n"
        "  --min-identity FLOAT Minimum alignment identity [0.70]\n"
        "  --gpu               Enable GPU processing\n"
        "  --gpu-device INT    GPU device ID [0]\n"
        "  -h, --help          Show this help\n"
        "\n"
        "Example:\n"
        "  llmap validate-real --reference chr14_igh.fa --reads reads.fq\n"
    );
}

struct ValidateRealArgs {
    std::string reference;
    std::string reads;
    std::string baseline;
    std::string truth;
    std::string output_dir;
    int k = 15;
    int w = 10;
    float min_identity = 0.70f;
    bool use_gpu = false;
    int gpu_device = 0;
    bool help = false;
};

bool parse_validate_real_args(int argc, char** argv, ValidateRealArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if (arg == "--reference" && i + 1 < argc) {
            args.reference = argv[++i];
        } else if (arg == "--reads" && i + 1 < argc) {
            args.reads = argv[++i];
        } else if (arg == "--baseline" && i + 1 < argc) {
            args.baseline = argv[++i];
        } else if (arg == "--truth" && i + 1 < argc) {
            args.truth = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            args.output_dir = argv[++i];
        } else if (arg == "--k" && i + 1 < argc) {
            args.k = std::stoi(argv[++i]);
        } else if (arg == "--w" && i + 1 < argc) {
            args.w = std::stoi(argv[++i]);
        } else if (arg == "--min-identity" && i + 1 < argc) {
            args.min_identity = std::stof(argv[++i]);
        } else if (arg == "--gpu") {
            args.use_gpu = true;
        } else if (arg == "--gpu-device" && i + 1 < argc) {
            args.gpu_device = std::stoi(argv[++i]);
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_validate_real(int argc, char** argv) {
    ValidateRealArgs args;
    if (!parse_validate_real_args(argc, argv, args)) {
        print_validate_real_usage();
        return 1;
    }

    if (args.help) {
        print_validate_real_usage();
        return 0;
    }

    if (args.reference.empty()) {
        std::fprintf(stderr, "Error: --reference is required\n");
        print_validate_real_usage();
        return 1;
    }

    if (args.reads.empty()) {
        std::fprintf(stderr, "Error: --reads is required\n");
        print_validate_real_usage();
        return 1;
    }

    validation::RealReferenceConfig config;
    config.reference_fasta = args.reference;
    config.reads_fastq = args.reads;
    config.minimap2_bam = args.baseline;
    config.ground_truth_bed = args.truth;
    config.output_dir = args.output_dir;
    config.minimizer_k = static_cast<uint8_t>(args.k);
    config.minimizer_w = static_cast<uint8_t>(args.w);
    config.min_identity = args.min_identity;
    config.use_gpu = args.use_gpu;
    config.gpu_device = args.gpu_device;

    if (!config.Validate()) {
        std::fprintf(stderr, "Configuration errors:\n%s",
                     config.ValidationErrors().c_str());
        return 1;
    }

    std::fprintf(stderr, "Running real reference validation...\n");
    std::fprintf(stderr, "  Reference: %s\n", args.reference.c_str());
    std::fprintf(stderr, "  Reads:     %s\n", args.reads.c_str());
    if (!args.baseline.empty()) {
        std::fprintf(stderr, "  Baseline:  %s\n", args.baseline.c_str());
    }
    std::fprintf(stderr, "\n");

    auto result = validation::RunRealReferenceValidation(config);

    std::puts(result.Summary().c_str());

    return result.overall_pass ? 0 : 1;
}

}  // namespace llmap::cli
