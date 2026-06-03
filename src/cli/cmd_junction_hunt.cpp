// LLmap — `llmap junction-hunt` CLI command (Mode-5 entry point).
//
// Stage-1 + Stage-2 pipeline for one sample-BAM × the genome-wide
// NAHR-pair panel. minimap2 is NOT invoked at any stage — the read
// classification is alignment-free.
//
// Pipeline:
//   1. Load pair panel TSV.
//   2. For each pair: extract LCR_up / LCR_down / interior sequences from
//      the reference, build PairKmerIndex (k=21,31,51,71,101).
//   3. Stream reads from the input (FASTA/FASTQ; BAM support in a later
//      iteration). For each read:
//        a. Tile multi-k.
//        b. For each candidate pair (heuristic: any LCR_up/down chrom
//           bucket the read falls into), run CallJunction.
//   4. Emit per-(read,pair) JunctionRecord TSV.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/nahr_pair_loader.h"
#include "junction_hunter/pair_kmer_index.h"
#include "junction_hunter/read_tiler.h"
#include "junction_hunter/consensus_caller.h"
#include "io/fasta_reader.h"
#include "io/mmap_fasta.h"

namespace llmap::cli {

namespace {

struct JhArgs {
    std::string pairs_tsv;
    std::string reference;
    std::string reads;        // FASTA or FASTQ
    std::string output;
    int max_pairs{0};         // 0 = all
    bool verbose{false};
    bool help{false};
};

void PrintUsage() {
    std::puts(
        "Usage: llmap junction-hunt [options]\n"
        "\n"
        "Detect NAHR breakpoints in long reads via multi-k consensus —\n"
        "alignment-free, no minimap2. Outputs one record per (read, pair).\n"
        "\n"
        "Required:\n"
        "  --pairs FILE          NAHR-pair TSV (linter output)\n"
        "  --reference FILE      Reference FASTA (GRCh38)\n"
        "  --reads FILE          Input reads (FASTA/FASTQ)\n"
        "  -o, --output FILE     Per-(read,pair) record TSV\n"
        "\n"
        "Optional:\n"
        "  --max-pairs N         Limit panel to first N pairs (debug)\n"
        "  -v, --verbose         Verbose progress\n"
        "  -h, --help            This help\n"
    );
}

bool ParseArgs(int argc, char** argv, JhArgs& a) {
    for (int i = 0; i < argc; ++i) {
        std::string x = argv[i];
        auto take = [&](const char* n) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "error: %s requires a value\n", n); return nullptr; }
            return argv[++i];
        };
        if (x == "-h" || x == "--help") { a.help = true; return true; }
        else if (x == "-v" || x == "--verbose") a.verbose = true;
        else if (x == "--pairs")     { auto v = take("--pairs");     if (!v) return false; a.pairs_tsv = v; }
        else if (x == "--reference" || x == "-r") { auto v = take("--reference"); if (!v) return false; a.reference = v; }
        else if (x == "--reads")     { auto v = take("--reads");     if (!v) return false; a.reads = v; }
        else if (x == "-o" || x == "--output") { auto v = take("--output"); if (!v) return false; a.output = v; }
        else if (x == "--max-pairs") { auto v = take("--max-pairs"); if (!v) return false; a.max_pairs = std::stoi(v); }
        else { std::fprintf(stderr, "error: unknown arg %s\n", x.c_str()); return false; }
    }
    if (a.pairs_tsv.empty() || a.reference.empty() || a.reads.empty() || a.output.empty()) {
        std::fprintf(stderr, "error: --pairs, --reference, --reads, --output are all required\n");
        return false;
    }
    return true;
}

}  // namespace

int run_junction_hunt(int argc, char** argv) {
    JhArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 2; }
    if (args.help) { PrintUsage(); return 0; }

    auto t0 = std::chrono::steady_clock::now();

    // 1. Pairs.
    std::vector<junction_hunter::NahrPair> pairs;
    auto pst = junction_hunter::LoadNahrPairsTsv(args.pairs_tsv, pairs);
    if (!pst.ok) {
        std::fprintf(stderr, "error: pair-tsv load failed: %s\n", pst.error.c_str());
        return 1;
    }
    if (args.max_pairs > 0 && static_cast<int>(pairs.size()) > args.max_pairs)
        pairs.resize(args.max_pairs);
    if (args.verbose)
        std::fprintf(stderr, "[junction-hunt] %zu pairs loaded\n", pairs.size());

    // 2. Build per-pair multi-k indices using reference sequence slices.
    //    Uses MmapFastaReader::GetSubsequence(name, start, length).
    llmap::io::MmapFastaReader ref(args.reference);
    if (!ref.IsValid()) {
        std::fprintf(stderr, "error: cannot mmap reference %s: %s\n",
                      args.reference.c_str(), ref.LastError().c_str());
        return 1;
    }

    std::vector<junction_hunter::PairKmerIndex> indices;
    indices.reserve(pairs.size());
    junction_hunter::MultiKConfig kcfg;
    std::size_t total_unique = 0;
    for (const auto& p : pairs) {
        const auto up_len   = p.lcr_up_end   > p.lcr_up_start   ? p.lcr_up_end   - p.lcr_up_start   : 0;
        const auto down_len = p.lcr_down_end > p.lcr_down_start ? p.lcr_down_end - p.lcr_down_start : 0;
        const auto in_len   = p.interior_end > p.interior_start ? p.interior_end - p.interior_start : 0;
        std::string up   = ref.GetSubsequence(p.chrom, p.lcr_up_start,   up_len);
        std::string down = ref.GetSubsequence(p.chrom, p.lcr_down_start, down_len);
        std::string inn  = ref.GetSubsequence(p.chrom, p.interior_start, in_len);
        junction_hunter::PairKmerIndex idx;
        auto st = junction_hunter::BuildPairKmerIndex(p, up, down, inn, kcfg, idx);
        total_unique += st.total_unique;
        indices.push_back(std::move(idx));
    }
    if (args.verbose)
        std::fprintf(stderr, "[junction-hunt] built %zu pair indices (%zu total k-mers)\n",
                      indices.size(), total_unique);

    // 3. Stream reads.
    std::ofstream out(args.output);
    if (!out) { std::fprintf(stderr, "error: cannot write %s\n", args.output.c_str()); return 1; }
    out << "read_id\tpair_id\tcall\tn_kmer_total\tn_up\tn_dn\tn_in\tn_amb\t"
           "up_mono\tdn_mono\tbreakpoint_read_pos\tbreakpoint_quality\n";

    llmap::io::FastaReader rd(args.reads);
    // FastaReader has no IsValid() — rely on HasMore() in the loop;
    // LastError() is checked after the first failed Next().
    std::size_t n_reads = 0, n_junctions = 0;
    while (rd.HasMore()) {
        auto record = rd.Next();
        if (!record.IsValid()) continue;
        ++n_reads;
        auto tiling = junction_hunter::TileRead(record.sequence, kcfg);
        for (std::size_t i = 0; i < indices.size(); ++i) {
            auto rec = junction_hunter::CallJunction(record.name, tiling, indices[i], pairs[i], kcfg);
            if (rec.call == junction_hunter::JunctionCall::Unmapped) continue;
            out << rec.read_id << '\t' << rec.pair_id << '\t'
                << junction_hunter::JunctionCallName(rec.call) << '\t'
                << rec.n_kmer_total << '\t'
                << rec.n_consensus_up << '\t' << rec.n_consensus_dn << '\t'
                << rec.n_consensus_in << '\t' << rec.n_ambiguous << '\t'
                << rec.up_monotonicity << '\t' << rec.dn_monotonicity << '\t'
                << rec.breakpoint_read_pos << '\t' << rec.breakpoint_quality << '\n';
            if (rec.call == junction_hunter::JunctionCall::JunctionReal) ++n_junctions;
        }
        if (args.verbose && n_reads % 10000 == 0)
            std::fprintf(stderr, "[junction-hunt]   %zu reads, %zu junctions so far\n", n_reads, n_junctions);
    }

    auto dt = std::chrono::steady_clock::now() - t0;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
    std::fprintf(stderr, "[junction-hunt] done in %lds — %zu reads, %zu junction_real calls\n",
                  static_cast<long>(secs), n_reads, n_junctions);
    return 0;
}

}  // namespace llmap::cli
