// LLmap — cmd_align argument parsing.

#include "cli/cmd_align_internal.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>

namespace llmap::cli::align_internal {

Preset ParsePreset(std::string_view name) {
    if (name == "map-hifi" || name == "hifi") return Preset::MapHifi;
    if (name == "map-ont" || name == "ont") return Preset::MapOnt;
    if (name == "map-pb" || name == "pb") return Preset::MapPb;
    if (name == "sr") return Preset::Sr;
    return Preset::None;
}

const char* PresetName(Preset preset) {
    switch (preset) {
        case Preset::MapHifi: return "map-hifi";
        case Preset::MapOnt:  return "map-ont";
        case Preset::MapPb:   return "map-pb";
        case Preset::Sr:      return "sr";
        case Preset::None:    return "none";
    }
    return "none";
}

void ApplyPreset(Preset preset, AlignArgs& args) {
    switch (preset) {
        case Preset::MapHifi:
            // PacBio HiFi. minimap2 -x map-hifi does not hard-filter by identity;
            // the empirical identity distribution on real HiFi vs an imperfect
            // reference (e.g., paralog-rich IGH) sits at mean ~0.72 because the
            // reference itself is collapsed. We follow the same policy: keep a
            // low floor and let MAPQ separate good from bad downstream.
            args.kmer_size = 19;
            args.window_size = 19;
            args.min_identity = 0.50f;
            args.min_chain = 40;
            break;

        case Preset::MapOnt:
        case Preset::MapPb:
            // ONT/legacy PacBio: higher error rate (~5-15%); same policy.
            args.kmer_size = 15;
            args.window_size = 10;
            args.min_identity = 0.40f;
            args.min_chain = 20;
            break;

        case Preset::Sr:
            // Short reads (Illumina): very high accuracy, short length.
            args.kmer_size = 21;
            args.window_size = 11;
            args.min_identity = 0.80f;   // SR alignments below this are noise.
            args.min_chain = 50;
            break;

        case Preset::None:
            // Keep defaults
            break;
    }
}

void PrintAlignUsage() {
    std::puts(
        "Usage: llmap align [options]\n"
        "\n"
        "Align reads to a reference genome.\n"
        "\n"
        "Required:\n"
        "  -r, --reads FILE        Input reads (FASTQ)\n"
        "  --reference FILE        Reference genome (FASTA)\n"
        "  -o, --output FILE       Output alignment file (SAM/BAM)\n"
        "\n"
        "Presets (recommended for best results):\n"
        "  -x PRESET               Read type preset:\n"
        "                            map-hifi: PacBio HiFi (k=19, w=19, identity=0.90)\n"
        "                            map-ont:  Oxford Nanopore (k=15, w=10, identity=0.70)\n"
        "                            map-pb:   Legacy PacBio CLR (same as map-ont)\n"
        "                            sr:       Short reads/Illumina (k=21, w=11, identity=0.95)\n"
        "\n"
        "Index caching:\n"
        "  -i, --index FILE        Use pre-built .llmi index (from `llmap index`)\n"
        "                          When provided, skips index building for faster startup\n"
        "\n"
        "Output format:\n"
        "  --bam                   Output BAM format (requires htslib)\n"
        "  --sam                   Output SAM format (default)\n"
        "  --parquet FILE          Also output probabilistic Parquet\n"
        "\n"
        "Alignment parameters (override preset values):\n"
        "  -k, --kmer INT          Minimizer k-mer size [15] (ignored if --index)\n"
        "  -w, --window INT        Minimizer window [10] (ignored if --index)\n"
        "  --min-chain INT         Minimum chain score [30]\n"
        "  --min-identity FLOAT    Minimum alignment identity [0.80]\n"
        "\n"
        "Performance:\n"
        "  -t, --threads INT       Number of threads [1]\n"
        "  --max-chains INT        Max chains to extend per read [10]\n"
        "\n"
        "PSV-based paralog disambiguation:\n"
        "  --psv-catalog FILE      PSV catalog (BED/VCF) for paralog assignment\n"
        "  --psv-weight FLOAT      Weight for PSV vs probabilistic assignment [0.50]\n"
        "  --psv-min-posterior F   Min posterior for confident paralog call [0.90]\n"
        "  --psv-only              Use only PSV-based assignment (skip probabilistic)\n"
        "\n"
        "LLM-assisted diagnostics:\n"
        "  --llm                   Enable Claude LLM for alignment diagnostics\n"
        "  --llm-api-key KEY       Anthropic API key (or set ANTHROPIC_API_KEY)\n"
        "  --llm-threshold FLOAT   Mapping rate threshold to trigger diagnostics [0.50]\n"
        "  --llm-work-dir DIR      Working directory for LLM artifacts\n"
        "\n"
        "Pipeline mode:\n"
        "  --classical-only        Pure seed-chain-extend mode (no probabilistic framework)\n"
        "                          Reduces memory footprint; disables --llm if set\n"
        "\n"
        "Other:\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Examples:\n"
        "  llmap align -x map-hifi -r reads.fastq --reference ref.fasta -o out.sam\n"
        "  llmap align -x map-ont -r reads.fastq --reference ref.fasta -o out.bam --bam\n"
        "  llmap align -r reads.fastq --reference ref.fasta -o out.sam  # default params\n"
        "  llmap align -r reads.fastq --reference ref.fasta -o out.sam -i ref.llmi\n"
        "  llmap align -r reads.fastq --reference ref.fasta -o out.sam --psv-catalog psv.bed\n"
    );
}

bool ParseAlignArgs(int argc, char** argv, AlignArgs& args) {
    // Track which args were explicitly set (for preset override detection)
    bool explicit_kmer = false, explicit_window = false;
    bool explicit_identity = false, explicit_chain = false;

    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if ((arg == "-r" || arg == "--reads") && i + 1 < argc) {
            args.reads = argv[++i];
        } else if (arg == "--reference" && i + 1 < argc) {
            args.reference = argv[++i];
        } else if (arg == "-x" && i + 1 < argc) {
            args.preset = ParsePreset(argv[++i]);
            if (args.preset == Preset::None) {
                std::fprintf(stderr, "Unknown preset: %s\n", argv[i]);
                std::fprintf(stderr, "Valid presets: map-hifi, map-ont, map-pb, sr\n");
                return false;
            }
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            args.output = argv[++i];
        } else if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            args.index = argv[++i];
        } else if (arg == "--parquet" && i + 1 < argc) {
            args.parquet_output = argv[++i];
        } else if (arg == "--bam") {
            args.use_bam = true;
            args.use_sam = false;
        } else if (arg == "--sam") {
            args.use_sam = true;
            args.use_bam = false;
        } else if ((arg == "-k" || arg == "--kmer") && i + 1 < argc) {
            args.kmer_size = std::stoi(argv[++i]);
            explicit_kmer = true;
        } else if ((arg == "-w" || arg == "--window") && i + 1 < argc) {
            args.window_size = std::stoi(argv[++i]);
            explicit_window = true;
        } else if (arg == "--min-chain" && i + 1 < argc) {
            args.min_chain = std::stoi(argv[++i]);
            explicit_chain = true;
        } else if (arg == "--min-identity" && i + 1 < argc) {
            args.min_identity = std::stof(argv[++i]);
            explicit_identity = true;
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            args.threads = std::stoi(argv[++i]);
        } else if (arg == "--max-chains" && i + 1 < argc) {
            args.max_chains = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            args.verbose = true;
        } else if (arg == "--llm") {
            args.enable_llm = true;
        } else if (arg == "--llm-api-key" && i + 1 < argc) {
            args.llm_api_key = argv[++i];
        } else if (arg == "--llm-threshold" && i + 1 < argc) {
            args.llm_threshold = std::stof(argv[++i]);
        } else if (arg == "--llm-work-dir" && i + 1 < argc) {
            args.llm_work_dir = argv[++i];
        } else if (arg == "--psv-catalog" && i + 1 < argc) {
            args.psv_catalog = argv[++i];
        } else if (arg == "--psv-weight" && i + 1 < argc) {
            args.psv_weight = std::stof(argv[++i]);
        } else if (arg == "--psv-min-posterior" && i + 1 < argc) {
            args.psv_min_posterior = std::stof(argv[++i]);
        } else if (arg == "--psv-only") {
            args.psv_only = true;
        } else if (arg == "--classical-only") {
            args.classical_only = true;
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }

    // Apply preset defaults, then restore explicit overrides
    if (args.preset != Preset::None) {
        int saved_kmer = args.kmer_size;
        int saved_window = args.window_size;
        int saved_chain = args.min_chain;
        float saved_identity = args.min_identity;

        ApplyPreset(args.preset, args);

        // Explicit CLI args override preset values
        if (explicit_kmer)     args.kmer_size = saved_kmer;
        if (explicit_window)   args.window_size = saved_window;
        if (explicit_chain)    args.min_chain = saved_chain;
        if (explicit_identity) args.min_identity = saved_identity;
    }

    return true;
}

}  // namespace llmap::cli::align_internal
