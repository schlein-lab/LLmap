// LLmap -- `llmap variant-ingest` CLI command (Layer 4 stub).
//
// Ingests a population variant catalog (VCF/BCF/BED/TSV) and writes a
// per-position prior table in the llmap_priors_v1 format documented at
// knowledge/variants/PRIOR_FORMAT.md.
//
// Usage:
//   llmap variant-ingest --source dbsnp
//                        --vcf input.vcf.gz
//                        --output dbsnp.priors
//                        [--bucket-bp 100]
//                        [--min-af 0.0]
//                        [--kind snv|sv|mixed]
//                        [--liftover-chain chain.gz]
//
// Not wired into llmap_main.cpp yet — this is a scaffold for the upcoming
// variant-priors feature. Run-functions print "not yet implemented".
//
// TODO(layer4):
//   - VCF/BCF reader (consider htslib dependency or our own minimal parser)
//   - SV-type aware bucketiser (overlap test for SV intervals, point
//     accumulate for SNV)
//   - bgzip + index sidecar writer
//   - liftover support via UCSC chain file
//   - multi-source merge (called by `variant-merge` later)

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

namespace llmap::cli {

namespace {

struct VariantIngestArgs {
    std::string source;          // catalog id (e.g. "dbsnp", "gnomad_sv")
    std::string vcf_path;        // input VCF/BCF/BED/TSV
    std::string output_path;     // output .priors
    std::string kind = "auto";   // snv | sv | mixed | auto
    std::string liftover_chain;  // optional
    uint32_t bucket_bp = 100;
    float min_af = 0.0f;
    bool tag_pathogenic = false; // for clinvar source
    bool no_af = false;          // for presence-only catalogs (e.g. DGV)
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap variant-ingest --source SRC_ID --vcf INPUT --output OUT.priors\n"
        "                            [--bucket-bp 100] [--min-af 0.0]\n"
        "                            [--kind snv|sv|mixed]\n"
        "                            [--liftover-chain CHAIN]\n"
        "                            [--no-af] [--tag-pathogenic]\n"
        "\n"
        "Ingests a population variant catalog and writes a per-position prior\n"
        "table consumed by `llmap align --variant-priors`.\n"
        "\n"
        "Supported sources: dbsnp, gnomad_genome_snv, gnomad_sv, dgv, dbvar,\n"
        "                   1kg_sv_phase3, hprc_panvariants, clinvar,\n"
        "                   audano_sv_2019, beyter_iceland_2021,\n"
        "                   mouse_genomes_project_v8\n"
        "\n"
        "See knowledge/variants/PRIOR_FORMAT.md for the on-disk layout.\n");
}

bool ParseArgs(int argc, char** argv, VariantIngestArgs& a) {
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--source"          && i + 1 < argc) a.source = argv[++i];
        else if (arg == "--vcf"             && i + 1 < argc) a.vcf_path = argv[++i];
        else if (arg == "--output"          && i + 1 < argc) a.output_path = argv[++i];
        else if (arg == "--kind"            && i + 1 < argc) a.kind = argv[++i];
        else if (arg == "--liftover-chain"  && i + 1 < argc) a.liftover_chain = argv[++i];
        else if (arg == "--bucket-bp"       && i + 1 < argc) a.bucket_bp = static_cast<uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--min-af"          && i + 1 < argc) a.min_af = std::stof(argv[++i]);
        else if (arg == "--no-af")                            a.no_af = true;
        else if (arg == "--tag-pathogenic")                   a.tag_pathogenic = true;
        else if (arg == "--help" || arg == "-h")              { a.help = true; return true; }
        else {
            std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_variant_ingest(int argc, char** argv) {
    VariantIngestArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 1; }
    if (args.help) { PrintUsage(); return 0; }

    if (args.source.empty() || args.vcf_path.empty() || args.output_path.empty()) {
        std::fprintf(stderr, "variant-ingest: --source, --vcf, --output are required.\n");
        PrintUsage();
        return 1;
    }

    std::fprintf(stderr,
        "[variant-ingest] not yet implemented.\n"
        "  source        = %s\n"
        "  input vcf     = %s\n"
        "  output priors = %s\n"
        "  bucket_bp     = %u\n"
        "  min_af        = %.4f\n"
        "  kind          = %s\n"
        "  liftover      = %s\n"
        "  no_af         = %s\n"
        "  tag_pathogenic= %s\n"
        "\n"
        "Scaffold only; see TODO(layer4) markers in src/cli/cmd_variant_ingest.cpp\n",
        args.source.c_str(), args.vcf_path.c_str(), args.output_path.c_str(),
        args.bucket_bp, args.min_af, args.kind.c_str(),
        args.liftover_chain.empty() ? "(none)" : args.liftover_chain.c_str(),
        args.no_af ? "yes" : "no",
        args.tag_pathogenic ? "yes" : "no");
    return 2;  // ENOSYS-style: command recognised, not yet implemented.
}

}  // namespace llmap::cli
