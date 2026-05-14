// LLmap — `llmap generate-synth` CLI command.

#include "cli/commands.h"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

#include "synthetic/benchmark_truth.h"
#include "synthetic/igh_locus_generator.h"

namespace llmap::cli {

namespace {

void print_generate_synth_usage() {
    std::puts(
        "Usage: llmap generate-synth [options]\n"
        "\n"
        "Generate synthetic data with ground truth for validation.\n"
        "\n"
        "Modes:\n"
        "  --task t1             Benchmark T1: synthetic WGS with positional truth\n"
        "  --task t2             Benchmark T2: paralog stress test with assignment truth\n"
        "  (no --task)           Legacy IGH locus generation (default)\n"
        "\n"
        "Required:\n"
        "  --output DIR          Output directory (created if needed)\n"
        "  --out-prefix PREFIX   Output prefix (legacy mode only)\n"
        "\n"
        "T1/T2 options:\n"
        "  --reads N             Number of reads [1000000 for T1, 500000 for T2]\n"
        "  --coverage N          Target coverage [30]\n"
        "  --seed N              Random seed [42]\n"
        "  --read-length N       Mean read length [10000]\n"
        "  --error-rate F        Error rate [0.001]\n"
        "\n"
        "T2 paralog options:\n"
        "  --paralog FAMILIES    Comma-separated: igh,mhc,nphp1 [all]\n"
        "\n"
        "Legacy options:\n"
        "  --preset NAME         Preset: tiny, canonical, balanced, dup5, dup10,\n"
        "                        dup30, dup50, dup100, seq-identical [canonical]\n"
        "  --locus-length INT    Locus length in bp [50000]\n"
        "  --n-psvs INT          Number of PSVs [20]\n"
        "  --dup-fraction FLOAT  Duplication fraction [0.0]\n"
        "  --depth INT           Coverage depth [30]\n"
        "\n"
        "  -h, --help            Show this help\n"
        "\n"
        "Output files (T1/T2 mode):\n"
        "  DIR/reads.fastq       Reads FASTQ\n"
        "  DIR/reference.fa      Reference FASTA\n"
        "  DIR/truth.tsv         Ground truth (read_id<TAB>chrom<TAB>pos)\n"
        "  DIR/truth_paralog.tsv Paralog truth (T2 only)\n"
        "\n"
        "Examples:\n"
        "  llmap generate-synth --task t1 --output /tmp/synth_t1 --reads 10000\n"
        "  llmap generate-synth --task t2 --output /tmp/synth_t2 --paralog igh,nphp1\n"
        "  llmap generate-synth --out-prefix test --preset balanced\n"
    );
}

struct GenerateSynthArgs {
    std::string task;                   // "t1", "t2", or empty for legacy
    std::string output_dir;             // --output for t1/t2
    std::string out_prefix;             // --out-prefix for legacy
    std::string preset{"canonical"};
    std::string paralog_families;       // comma-separated: igh,mhc,nphp1
    std::uint64_t seed{42};
    std::uint64_t n_reads{0};           // 0 = use default
    std::uint64_t locus_length{50000};
    std::uint32_t read_length{10000};
    std::uint32_t n_psvs{20};
    std::uint32_t coverage{30};
    float dup_fraction{0.0f};
    float error_rate{0.001f};
    std::uint32_t depth{30};
    bool help{false};
};

bool parse_generate_synth_args(int argc, char** argv, GenerateSynthArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if (arg == "--task" && i + 1 < argc) {
            args.task = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            args.output_dir = argv[++i];
        } else if (arg == "--out-prefix" && i + 1 < argc) {
            args.out_prefix = argv[++i];
        } else if (arg == "--preset" && i + 1 < argc) {
            args.preset = argv[++i];
        } else if (arg == "--paralog" && i + 1 < argc) {
            args.paralog_families = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            args.seed = std::stoull(argv[++i]);
        } else if (arg == "--reads" && i + 1 < argc) {
            args.n_reads = std::stoull(argv[++i]);
        } else if (arg == "--locus-length" && i + 1 < argc) {
            args.locus_length = std::stoull(argv[++i]);
        } else if ((arg == "--read-length" || arg == "--readlen") && i + 1 < argc) {
            args.read_length = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--n-psvs" && i + 1 < argc) {
            args.n_psvs = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--dup-fraction" && i + 1 < argc) {
            args.dup_fraction = std::stof(argv[++i]);
        } else if (arg == "--depth" && i + 1 < argc) {
            args.depth = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--coverage" && i + 1 < argc) {
            args.coverage = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        } else if (arg == "--error-rate" && i + 1 < argc) {
            args.error_rate = std::stof(argv[++i]);
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

int run_benchmark_t1(const GenerateSynthArgs& args) {
    synthetic::T1Config config;
    config.seed = args.seed;
    config.n_reads = args.n_reads > 0 ? args.n_reads : 1000000;
    config.coverage = args.coverage;
    config.read_length = args.read_length;
    config.error_rate = args.error_rate;

    std::fprintf(stderr, "Generating T1 (WGS) benchmark dataset...\n");
    std::fprintf(stderr, "  Seed:        %zu\n", static_cast<size_t>(config.seed));
    std::fprintf(stderr, "  Reads:       %zu\n", static_cast<size_t>(config.n_reads));
    std::fprintf(stderr, "  Coverage:    %u×\n", config.coverage);
    std::fprintf(stderr, "  Read length: %u bp\n", config.read_length);
    std::fprintf(stderr, "  Error rate:  %.4f\n", config.error_rate);
    std::fprintf(stderr, "\n");

    auto dataset = synthetic::BenchmarkGenerator::generate_t1(config);

    std::filesystem::create_directories(args.output_dir);
    std::string reads_path = args.output_dir + "/reads.fastq";
    std::string ref_path = args.output_dir + "/reference.fa";
    std::string truth_path = args.output_dir + "/truth.tsv";

    synthetic::BenchmarkGenerator::write_fastq(dataset, reads_path);
    synthetic::BenchmarkGenerator::write_reference(dataset, ref_path);
    synthetic::BenchmarkGenerator::write_truth_tsv(dataset, truth_path);

    std::fprintf(stderr, "Generated:\n");
    std::fprintf(stderr, "  Reads:     %s (%zu reads, %zu bp)\n",
                 reads_path.c_str(), dataset.total_reads, dataset.total_bases);
    std::fprintf(stderr, "  Reference: %s (%zu sequences)\n",
                 ref_path.c_str(), dataset.references.size());
    std::fprintf(stderr, "  Truth:     %s\n", truth_path.c_str());

    return 0;
}

int run_benchmark_t2(const GenerateSynthArgs& args) {
    synthetic::T2Config config;
    config.seed = args.seed;
    config.n_reads = args.n_reads > 0 ? args.n_reads : 500000;
    config.coverage = args.coverage;
    config.read_length = args.read_length;
    config.error_rate = args.error_rate;

    // Parse paralog families.
    if (args.paralog_families.empty() || args.paralog_families == "all") {
        config.families.push_back(synthetic::paralog_presets::igh_constant());
        config.families.push_back(synthetic::paralog_presets::nphp1());
        config.families.push_back(synthetic::paralog_presets::mhc_class1());
    } else {
        std::string families = args.paralog_families;
        std::size_t pos = 0;
        while ((pos = families.find(',')) != std::string::npos) {
            std::string fam = families.substr(0, pos);
            families = families.substr(pos + 1);
            if (fam == "igh") {
                config.families.push_back(synthetic::paralog_presets::igh_constant());
            } else if (fam == "nphp1") {
                config.families.push_back(synthetic::paralog_presets::nphp1());
            } else if (fam == "mhc") {
                config.families.push_back(synthetic::paralog_presets::mhc_class1());
            }
        }
        if (!families.empty()) {
            if (families == "igh") {
                config.families.push_back(synthetic::paralog_presets::igh_constant());
            } else if (families == "nphp1") {
                config.families.push_back(synthetic::paralog_presets::nphp1());
            } else if (families == "mhc") {
                config.families.push_back(synthetic::paralog_presets::mhc_class1());
            }
        }
    }

    if (config.families.empty()) {
        std::fprintf(stderr, "Error: no valid paralog families specified\n");
        return 1;
    }

    std::fprintf(stderr, "Generating T2 (paralog stress) benchmark dataset...\n");
    std::fprintf(stderr, "  Seed:        %zu\n", static_cast<size_t>(config.seed));
    std::fprintf(stderr, "  Reads:       %zu\n", static_cast<size_t>(config.n_reads));
    std::fprintf(stderr, "  Coverage:    %u×\n", config.coverage);
    std::fprintf(stderr, "  Read length: %u bp\n", config.read_length);
    std::fprintf(stderr, "  Families:    %zu\n", config.families.size());
    for (const auto& fam : config.families) {
        std::fprintf(stderr, "    - %s (%zu members)\n",
                     fam.name.c_str(), fam.members.size());
    }
    std::fprintf(stderr, "\n");

    auto dataset = synthetic::BenchmarkGenerator::generate_t2(config);

    std::filesystem::create_directories(args.output_dir);
    std::string reads_path = args.output_dir + "/reads.fastq";
    std::string ref_path = args.output_dir + "/reference.fa";
    std::string truth_path = args.output_dir + "/truth.tsv";
    std::string paralog_truth_path = args.output_dir + "/truth_paralog.tsv";

    synthetic::BenchmarkGenerator::write_fastq(dataset, reads_path);
    synthetic::BenchmarkGenerator::write_reference(dataset, ref_path);
    synthetic::BenchmarkGenerator::write_truth_tsv(dataset, truth_path);
    synthetic::BenchmarkGenerator::write_paralog_truth_tsv(dataset, paralog_truth_path);

    std::fprintf(stderr, "Generated:\n");
    std::fprintf(stderr, "  Reads:         %s (%zu reads, %zu bp)\n",
                 reads_path.c_str(), dataset.total_reads, dataset.total_bases);
    std::fprintf(stderr, "  Reference:     %s (%zu sequences)\n",
                 ref_path.c_str(), dataset.references.size());
    std::fprintf(stderr, "  Truth:         %s\n", truth_path.c_str());
    std::fprintf(stderr, "  Paralog truth: %s\n", paralog_truth_path.c_str());

    return 0;
}

int run_legacy_synth(const GenerateSynthArgs& args) {
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

    // Dispatch based on --task.
    if (args.task == "t1" || args.task == "T1") {
        if (args.output_dir.empty()) {
            std::fprintf(stderr, "Error: --output is required for --task t1\n");
            print_generate_synth_usage();
            return 1;
        }
        return run_benchmark_t1(args);
    } else if (args.task == "t2" || args.task == "T2") {
        if (args.output_dir.empty()) {
            std::fprintf(stderr, "Error: --output is required for --task t2\n");
            print_generate_synth_usage();
            return 1;
        }
        return run_benchmark_t2(args);
    } else if (!args.task.empty()) {
        std::fprintf(stderr, "Error: unknown task '%s' (use t1 or t2)\n",
                     args.task.c_str());
        print_generate_synth_usage();
        return 1;
    }

    // Legacy mode.
    if (args.out_prefix.empty()) {
        std::fprintf(stderr, "Error: --out-prefix is required (or use --task t1/t2)\n");
        print_generate_synth_usage();
        return 1;
    }

    return run_legacy_synth(args);
}

}  // namespace llmap::cli
