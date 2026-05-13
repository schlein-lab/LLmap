// LLmap — `llmap generate-synth` CLI command.

#include "cli/commands.h"

#include <cstdio>
#include <cstring>
#include <string>

#include "synthetic/igh_locus_generator.h"

namespace llmap::cli {

namespace {

void print_generate_synth_usage() {
    std::puts(
        "Usage: llmap generate-synth [options]\n"
        "\n"
        "Generate synthetic IGH locus reference and reads with ground truth.\n"
        "\n"
        "Required:\n"
        "  --out-prefix PREFIX   Output prefix for files\n"
        "\n"
        "Optional:\n"
        "  --preset NAME         Preset: tiny, canonical, balanced, dup5, dup10,\n"
        "                        dup30, dup50, dup100, seq-identical [canonical]\n"
        "  --seed INT            Random seed [42]\n"
        "  --locus-length INT    Locus length in bp [50000]\n"
        "  --read-length INT     Mean read length [10000]\n"
        "  --n-psvs INT          Number of PSVs [20]\n"
        "  --dup-fraction FLOAT  Duplication fraction [0.0]\n"
        "  --depth INT           Coverage depth [30]\n"
        "  -h, --help            Show this help\n"
        "\n"
        "Output files:\n"
        "  PREFIX_reference.fa   Reference FASTA\n"
        "  PREFIX_reads.fq       Reads FASTQ\n"
        "  PREFIX_truth.tsv      Ground truth TSV\n"
        "\n"
        "Example:\n"
        "  llmap generate-synth --out-prefix test --preset balanced\n"
    );
}

struct GenerateSynthArgs {
    std::string out_prefix;
    std::string preset{"canonical"};
    uint64_t seed{42};
    uint64_t locus_length{50000};
    uint32_t read_length{10000};
    uint32_t n_psvs{20};
    float dup_fraction{0.0f};
    uint32_t depth{30};
    bool help{false};
};

bool parse_generate_synth_args(int argc, char** argv, GenerateSynthArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if (arg == "--out-prefix" && i + 1 < argc) {
            args.out_prefix = argv[++i];
        } else if (arg == "--preset" && i + 1 < argc) {
            args.preset = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            args.seed = std::stoull(argv[++i]);
        } else if (arg == "--locus-length" && i + 1 < argc) {
            args.locus_length = std::stoull(argv[++i]);
        } else if (arg == "--read-length" && i + 1 < argc) {
            args.read_length = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--n-psvs" && i + 1 < argc) {
            args.n_psvs = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--dup-fraction" && i + 1 < argc) {
            args.dup_fraction = std::stof(argv[++i]);
        } else if (arg == "--depth" && i + 1 < argc) {
            args.depth = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_generate_synth(int argc, char** argv) {
    GenerateSynthArgs args;
    if (!parse_generate_synth_args(argc, argv, args)) {
        print_generate_synth_usage();
        return 1;
    }

    if (args.help) {
        print_generate_synth_usage();
        return 0;
    }

    if (args.out_prefix.empty()) {
        std::fprintf(stderr, "Error: --out-prefix is required\n");
        print_generate_synth_usage();
        return 1;
    }

    synthetic::IGHLocusConfig config;
    config.seed = args.seed;

    if (args.preset == "tiny") {
        config = synthetic::presets::tiny_test(args.seed);
    } else if (args.preset == "canonical") {
        config = synthetic::presets::canonical_only(args.seed);
    } else if (args.preset == "balanced") {
        config = synthetic::presets::balanced_mosaic(args.seed);
    } else if (args.preset == "dup5") {
        config = synthetic::presets::dup_fraction_5(args.seed);
    } else if (args.preset == "dup10") {
        config = synthetic::presets::dup_fraction_10(args.seed);
    } else if (args.preset == "dup30") {
        config = synthetic::presets::dup_fraction_30(args.seed);
    } else if (args.preset == "dup50") {
        config = synthetic::presets::dup_fraction_50(args.seed);
    } else if (args.preset == "dup100") {
        config = synthetic::presets::dup_fraction_100(args.seed);
    } else if (args.preset == "seq-identical") {
        config = synthetic::presets::seq_identical_stress(args.seed);
    } else {
        config.locus_length = args.locus_length;
        config.read_length = args.read_length;
        config.n_psvs = args.n_psvs;
        config.mosaic.dup_fraction = args.dup_fraction;
        config.mosaic.canonical_depth = args.depth;
    }

    std::fprintf(stderr, "Generating synthetic IGH locus data...\n");
    std::fprintf(stderr, "  Preset:       %s\n", args.preset.c_str());
    std::fprintf(stderr, "  Seed:         %zu\n", static_cast<size_t>(config.seed));
    std::fprintf(stderr, "  Locus length: %zu bp\n",
                 static_cast<size_t>(config.locus_length));
    std::fprintf(stderr, "  Read length:  %u bp\n", config.read_length);
    std::fprintf(stderr, "  PSVs:         %u\n", config.n_psvs);
    std::fprintf(stderr, "  Dup fraction: %.2f\n", config.mosaic.dup_fraction);
    std::fprintf(stderr, "\n");

    synthetic::IGHLocusGenerator generator(config);
    auto dataset = generator.generate();

    std::string ref_path = args.out_prefix + "_reference.fa";
    std::string reads_path = args.out_prefix + "_reads.fq";
    std::string truth_path = args.out_prefix + "_truth.tsv";

    {
        std::FILE* fref = std::fopen(ref_path.c_str(), "w");
        if (!fref) {
            std::fprintf(stderr, "Error: cannot open %s for writing\n",
                        ref_path.c_str());
            return 1;
        }
        std::fprintf(fref, ">synthetic_igh_locus\n%s\n",
                    dataset.locus.canonical_sequence.c_str());
        std::fclose(fref);
    }

    synthetic::IGHLocusGenerator::write_fastq(dataset.reads, reads_path);
    synthetic::IGHLocusGenerator::write_ground_truth(dataset, truth_path);

    std::fprintf(stderr, "Generated:\n");
    std::fprintf(stderr, "  Reference:    %s (%zu bp)\n", ref_path.c_str(),
                 dataset.locus.canonical_sequence.size());
    std::fprintf(stderr, "  Reads:        %s (%zu reads)\n", reads_path.c_str(),
                 dataset.reads.size());
    std::fprintf(stderr, "  Ground truth: %s\n", truth_path.c_str());
    std::fprintf(stderr, "  Canonical:    %zu reads\n",
                 static_cast<size_t>(dataset.n_canonical_reads));
    std::fprintf(stderr, "  Duplicate:    %zu reads\n",
                 static_cast<size_t>(dataset.n_dup_reads));
    std::fprintf(stderr, "  PSV obs:      %zu\n",
                 static_cast<size_t>(dataset.total_psv_observations));

    return 0;
}

}  // namespace llmap::cli
