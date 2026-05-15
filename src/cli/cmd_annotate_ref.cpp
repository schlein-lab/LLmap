// LLmap -- `llmap annotate-ref` CLI command.
//
// Reads a reference FASTA, extracts local features per window, applies the
// rule-based classifier, and writes a .regann file to disk.
//
// Usage:
//   llmap annotate-ref --reference REF.fa \
//                      --rules organisms/human/classifier_rules.json \
//                      --output REF.regann \
//                      [--window-bp 1000]

#include "cli/commands.h"

#include "annot/annotation_store.h"
#include "annot/classifier.h"
#include "annot/feature_extractor.h"
#include "io/fasta_reader.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>

namespace llmap::cli {

namespace {

struct AnnotateArgs {
    std::string reference;
    std::string rules;
    std::string output;
    uint32_t window_bp = 1000;
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap annotate-ref --reference REF.fa --rules RULES.json \\\n"
        "                          --output REF.regann [--window-bp 1000]\n"
        "\n"
        "Computes per-window features over the reference, classifies each\n"
        "window via the organism rules file, and writes a .regann annotation\n"
        "file used by `llmap align` for region-aware mapping.\n");
}

bool ParseArgs(int argc, char** argv, AnnotateArgs& a) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--reference" && i + 1 < argc)       a.reference = argv[++i];
        else if (arg == "--rules" && i + 1 < argc)      a.rules = argv[++i];
        else if (arg == "--output" && i + 1 < argc)     a.output = argv[++i];
        else if (arg == "--window-bp" && i + 1 < argc)  a.window_bp = std::stoul(argv[++i]);
        else if (arg == "--help" || arg == "-h")        { a.help = true; return true; }
        else {
            std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int run_annotate_ref(int argc, char** argv) {
    AnnotateArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 1; }
    if (args.help) { PrintUsage(); return 0; }

    if (args.reference.empty() || args.rules.empty() || args.output.empty()) {
        std::fprintf(stderr, "Missing required arg.\n");
        PrintUsage();
        return 1;
    }
    if (!std::filesystem::exists(args.reference)) {
        std::fprintf(stderr, "reference not found: %s\n", args.reference.c_str());
        return 1;
    }
    if (!std::filesystem::exists(args.rules)) {
        std::fprintf(stderr, "rules file not found: %s\n", args.rules.c_str());
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    std::printf("Loading reference %s\n", args.reference.c_str());
    io::FastaReader rdr(args.reference);
    std::vector<std::string> contig_names;
    std::vector<std::string> contig_seqs;
    while (rdr.HasMore()) {
        auto rec = rdr.Next();
        if (rec.IsValid()) {
            contig_names.push_back(rec.name);
            contig_seqs.push_back(std::move(rec.sequence));
        }
    }
    std::printf("  %zu contigs\n", contig_seqs.size());

    std::printf("Loading classifier rules %s\n", args.rules.c_str());
    auto classifier = annot::Classifier::Load(args.rules);
    if (!classifier) return 1;
    std::printf("  %zu rules\n", classifier->NumRules());

    std::printf("Extracting features (window=%u bp)...\n", args.window_bp);
    annot::FeatureExtractorConfig cfg;
    cfg.window_bp = args.window_bp;
    cfg.step_bp = args.window_bp;

    auto features = annot::ExtractFeaturesAll(contig_seqs, cfg);
    std::printf("  %zu windows\n", features.size());

    std::printf("Classifying...\n");
    auto intervals = classifier->Classify(features);
    std::printf("  %zu intervals after merge\n", intervals.size());

    auto store = annot::AnnotationStore::Create(std::move(intervals),
                                                contig_names);
    if (!store->Save(args.output)) {
        std::fprintf(stderr, "failed to write %s\n", args.output.c_str());
        return 1;
    }

    auto t1 = std::chrono::steady_clock::now();
    float secs = std::chrono::duration<float>(t1 - t0).count();
    std::printf("Wrote %s (%zu intervals)\n", args.output.c_str(), store->Size());
    std::printf("Total time: %.1f s\n", secs);
    return 0;
}

}  // namespace llmap::cli
