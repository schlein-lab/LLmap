// LLmap — cmd_align argument parsing.

#include "cli/cmd_align_internal.h"

#include <cstdio>
#include <string>

namespace llmap::cli::align_internal {

void PrintAlignUsage() {
    std::puts(
        "Usage: llmap align [options]\n"
        "\n"
        "Align reads to a reference genome.\n"
        "\n"
        "Required:\n"
        "  -r, --reads FILE        Input reads (FASTQ)\n"
        "  -x, --reference FILE    Reference genome (FASTA)\n"
        "  -o, --output FILE       Output alignment file (SAM/BAM)\n"
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
        "Alignment parameters:\n"
        "  -k, --kmer INT          Minimizer k-mer size [15] (ignored if --index)\n"
        "  -w, --window INT        Minimizer window [10] (ignored if --index)\n"
        "  --min-chain INT         Minimum chain score [30]\n"
        "  --min-identity FLOAT    Minimum alignment identity [0.70]\n"
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
        "Other:\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Example:\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.bam --bam --parquet out.parquet\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam -i ref.llmi\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam --psv-catalog psv.bed\n"
        "  llmap align -r reads.fastq -x ref.fasta -o out.sam --llm\n"
    );
}

bool ParseAlignArgs(int argc, char** argv, AlignArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            args.help = true;
            return true;
        } else if ((arg == "-r" || arg == "--reads") && i + 1 < argc) {
            args.reads = argv[++i];
        } else if ((arg == "-x" || arg == "--reference") && i + 1 < argc) {
            args.reference = argv[++i];
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
        } else if ((arg == "-w" || arg == "--window") && i + 1 < argc) {
            args.window_size = std::stoi(argv[++i]);
        } else if (arg == "--min-chain" && i + 1 < argc) {
            args.min_chain = std::stoi(argv[++i]);
        } else if (arg == "--min-identity" && i + 1 < argc) {
            args.min_identity = std::stof(argv[++i]);
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
        } else if (arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace llmap::cli::align_internal
