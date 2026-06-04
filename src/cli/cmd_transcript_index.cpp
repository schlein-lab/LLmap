// LLmap — `llmap transcript-index` CLI command.
//
// Builds a TranscriptKmerIndex from a GENCODE GFF3 anchor source plus
// the built-in splice-site PWMs, and serialises it to disk. The index
// can later be loaded by the upcoming junction-hunter / chimera-aware
// alignment paths (Mode-5 and beyond).
//
// Wired in 2026-06-03 — the engine itself (anchor + index + annot
// transcript_kind_classifier + splice_site_db) was already built but
// unreachable from the CLI surface.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "anchor/anchor_store.h"
#include "annot/splice_site_db.h"
#include "index/transcript_kmer_index.h"

namespace llmap::cli {

namespace {

struct TxIndexArgs {
    std::string gencode_gff;
    std::string reference_fa;
    std::string mane_summary;
    std::string output;
    std::uint8_t k_intra{51};
    std::uint8_t k_junction{31};
    std::uint8_t k_alt{21};
    std::uint32_t max_occ{200};
    bool include_premrna{false};
    bool verbose{false};
    bool help{false};
};

void PrintUsage() {
    std::puts(
        "Usage: llmap transcript-index [options]\n"
        "\n"
        "Build an exon-boundary-aware reverse k-mer index from a GENCODE GFF3\n"
        "anchor source. The index distinguishes IntraExon / JunctionSpanning /\n"
        "BackSpliceSpanning / SterileIntronic / PreMrnaIntronic / ShortRna\n"
        "origins so junction k-mers are not falsely matched against genomic DNA.\n"
        "\n"
        "Required:\n"
        "  --gencode FILE          GENCODE annotation GFF3(.gz)\n"
        "  -o, --output FILE       Output index file (.llmtx)\n"
        "\n"
        "Optional anchor sources:\n"
        "  --reference FILE        Reference FASTA (extracts anchor sequences;\n"
        "                           required for k-mer build — omit only for\n"
        "                           coordinate-only metadata dumps)\n"
        "  --mane FILE             MANE summary TSV\n"
        "\n"
        "Index parameters:\n"
        "  --k-intra INT           k for intra-exon k-mers [51]\n"
        "  --k-junction INT        k for junction-spanning k-mers [31]\n"
        "  --k-alt INT             k for short-RNA fallback [21]\n"
        "  --max-occ INT           Cap occurrences per hash [200]\n"
        "  --include-premrna       Index PreMrnaIntronic k-mers (large)\n"
        "\n"
        "Other:\n"
        "  -v, --verbose           Verbose output\n"
        "  -h, --help              Show this help\n"
        "\n"
        "Example:\n"
        "  llmap transcript-index --gencode gencode.v46.annotation.gff3.gz \\\n"
        "                         --mane MANE.GRCh38.summary.txt \\\n"
        "                         -o gencode.v46.llmtx\n"
    );
}

bool ParseArgs(int argc, char** argv, TxIndexArgs& args) {
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto take_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires a value\n", name);
                return nullptr;
            }
            return argv[++i];
        };
        if (a == "-h" || a == "--help") { args.help = true; return true; }
        else if (a == "-v" || a == "--verbose") args.verbose = true;
        else if (a == "--gencode")     { auto v = take_value("--gencode");     if (!v) return false; args.gencode_gff = v; }
        else if (a == "--reference" || a == "-r") { auto v = take_value("--reference"); if (!v) return false; args.reference_fa = v; }
        else if (a == "--mane")        { auto v = take_value("--mane");        if (!v) return false; args.mane_summary = v; }
        else if (a == "-o" || a == "--output") { auto v = take_value("--output"); if (!v) return false; args.output = v; }
        else if (a == "--k-intra")     { auto v = take_value("--k-intra");     if (!v) return false; args.k_intra = static_cast<std::uint8_t>(std::stoi(v)); }
        else if (a == "--k-junction")  { auto v = take_value("--k-junction");  if (!v) return false; args.k_junction = static_cast<std::uint8_t>(std::stoi(v)); }
        else if (a == "--k-alt")       { auto v = take_value("--k-alt");       if (!v) return false; args.k_alt = static_cast<std::uint8_t>(std::stoi(v)); }
        else if (a == "--max-occ")     { auto v = take_value("--max-occ");     if (!v) return false; args.max_occ = static_cast<std::uint32_t>(std::stoul(v)); }
        else if (a == "--include-premrna") args.include_premrna = true;
        else { std::fprintf(stderr, "error: unknown arg %s\n", a.c_str()); return false; }
    }
    if (args.gencode_gff.empty() || args.output.empty()) {
        std::fprintf(stderr, "error: --gencode and --output are required\n");
        return false;
    }
    return true;
}

}  // namespace

int run_transcript_index(int argc, char** argv) {
    TxIndexArgs args;
    if (!ParseArgs(argc, argv, args)) {
        PrintUsage();
        return 2;
    }
    if (args.help) { PrintUsage(); return 0; }

    auto t0 = std::chrono::steady_clock::now();

    // 1. AnchorStore — load GENCODE annotations (+ optional MANE).
    anchor::AnchorStore store;
    if (args.verbose) std::fprintf(stderr, "[transcript-index] loading GENCODE %s\n", args.gencode_gff.c_str());
    const bool with_seq = !args.reference_fa.empty();
    auto status = store.LoadGencodeGff(args.gencode_gff, args.reference_fa, with_seq);
    if (!status.ok) {
        std::fprintf(stderr, "error: GENCODE load failed: %s (loaded=%zu skipped=%zu)\n",
                      status.error.c_str(), status.records_loaded, status.records_skipped);
        return 1;
    }
    if (args.verbose) std::fprintf(stderr, "[transcript-index]   loaded %zu anchors (skipped %zu)\n",
                                     status.records_loaded, status.records_skipped);
    if (!with_seq) {
        std::fprintf(stderr, "warning: no --reference given; anchor sequences are unavailable so the\n"
                              "         resulting index will be empty (k-mer build needs sequences).\n");
    }

    if (!args.mane_summary.empty()) {
        if (args.verbose) std::fprintf(stderr, "[transcript-index] loading MANE %s\n", args.mane_summary.c_str());
        auto mst = store.LoadMane(args.mane_summary);
        if (!mst.ok) {
            std::fprintf(stderr, "warning: MANE load failed: %s\n", mst.error.c_str());
        }
    }
    store.Reindex();

    // 2. SpliceSiteDb — built-in PWMs only (Burge & Sharp, Sheth, Roca).
    annot::SpliceSiteDb splice;
    splice.LoadDefaults();

    // 3. TranscriptKmerIndex — segregated k-mer hash over the anchor store.
    index::TranscriptKmerIndexConfig cfg;
    cfg.k_intra = args.k_intra;
    cfg.k_junction = args.k_junction;
    cfg.alt_k = args.k_alt;
    cfg.max_occ = args.max_occ;
    cfg.include_premrna_intronic = args.include_premrna;

    index::TranscriptKmerIndex idx;
    if (args.verbose) std::fprintf(stderr, "[transcript-index] building (k_intra=%u, k_junction=%u)\n", cfg.k_intra, cfg.k_junction);
    idx.BuildFromAnchorStore(store, splice, cfg);

    if (args.verbose) {
        std::fprintf(stderr,
            "[transcript-index]   IntraExon=%zu Junction=%zu BackSplice=%zu Sterile=%zu PreMrna=%zu ShortRna=%zu  total=%zu\n",
            idx.TableSize(index::KmerOrigin::IntraExon),
            idx.TableSize(index::KmerOrigin::JunctionSpanning),
            idx.TableSize(index::KmerOrigin::BackSpliceSpanning),
            idx.TableSize(index::KmerOrigin::SterileIntronic),
            idx.TableSize(index::KmerOrigin::PreMrnaIntronic),
            idx.TableSize(index::KmerOrigin::ShortRna),
            idx.TotalKmers());
    }

    // 4. Serialise.
    if (args.verbose) std::fprintf(stderr, "[transcript-index] saving %s\n", args.output.c_str());
    if (!idx.Save(args.output)) {
        std::fprintf(stderr, "error: failed to save %s\n", args.output.c_str());
        return 1;
    }

    auto dt = std::chrono::steady_clock::now() - t0;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
    std::fprintf(stderr, "[transcript-index] done in %lds — %zu k-mers across 6 origin tables\n",
                  static_cast<long>(secs), idx.TotalKmers());
    return 0;
}

}  // namespace llmap::cli
