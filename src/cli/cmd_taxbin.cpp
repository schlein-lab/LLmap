// LLmap — `llmap taxbin` CLI command (Mode-6 entry point).
//
// Mode-6 classifies un-human reads LLmap-natively (no Kraken/BLAST). The full
// pipeline is:
//   Stage 1  self_interference::AllpairPipeline  -> read-to-read clusters
//   Cascade  classical::ClassicalPipeline per species panel entry
//            (specific -> broad; bound reads short-circuit)            -> conf matrix
//   Engine   classify::Taxbin -> likelihood vectors + cluster-consensus collapse
//
// This first version (v0) consumes a precomputed per-read x species confidence
// matrix (--conf) plus optional Stage-1 clusters (--clusters) and runs the
// engine end-to-end, emitting the full likelihood distribution + hard calls.
// The in-process integration of AllpairPipeline (Stage 1) and per-species
// ClassicalPipeline (cascade) is the next step (see LLMAP_MODE6_SPEC.md) and
// will populate TaxbinInput directly instead of reading it from a TSV.
#include "cli/commands.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "classify/taxbin.h"

namespace llmap::cli {
namespace {

struct TaxbinArgs {
    std::string conf;        // per-read x species confidence matrix (TSV)
    std::string clusters;    // optional read_id -> cluster_id (TSV)
    std::string output;
    std::string clusters_out;
    float bind_threshold = 0.80f;
    float cluster_weight = 0.50f;
    float min_margin = 0.10f;
    bool verbose = false;
    bool help = false;
};

void PrintUsage() {
    std::puts(
        "Usage: llmap taxbin [options]\n"
        "\n"
        "Mode-6: LLmap-native taxonomic binding of un-human reads.\n"
        "Turns a per-read x species confidence matrix into an explicit\n"
        "likelihood distribution (incl. a NOVEL/dark class) and collapses\n"
        "ambiguous reads onto their read-to-read cluster's consensus species.\n"
        "\n"
        "Required:\n"
        "  --conf FILE         Per-read x species confidence matrix (TSV;\n"
        "                      header: read_id<TAB>sp1<TAB>sp2..., values in [0,1])\n"
        "  -o, --output FILE   Per-read calls + likelihood TSV\n"
        "\n"
        "Optional:\n"
        "  --clusters FILE     read_id<TAB>cluster_id from Stage 1 (allpair)\n"
        "  --clusters-out FILE Per-cluster consensus summary TSV\n"
        "  --bind-threshold F  Species binding threshold (default 0.80)\n"
        "  --cluster-weight F  Cluster-consensus weight in posterior (default 0.50)\n"
        "  --min-margin F      Min top1-top2 margin for a confident call (default 0.10)\n"
        "  -v, --verbose       Verbose progress\n"
        "  -h, --help          This help\n");
}

bool ParseArgs(int argc, char** argv, TaxbinArgs& a) {
    for (int i = 0; i < argc; ++i) {
        std::string s = argv[i];
        auto next = [&](const char* name) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "taxbin: %s needs a value\n", name); return nullptr; }
            return argv[++i];
        };
        if (s == "-h" || s == "--help") { a.help = true; return true; }
        else if (s == "-v" || s == "--verbose") a.verbose = true;
        else if (s == "--conf") { auto v = next("--conf"); if (!v) return false; a.conf = v; }
        else if (s == "-o" || s == "--output") { auto v = next("--output"); if (!v) return false; a.output = v; }
        else if (s == "--clusters") { auto v = next("--clusters"); if (!v) return false; a.clusters = v; }
        else if (s == "--clusters-out") { auto v = next("--clusters-out"); if (!v) return false; a.clusters_out = v; }
        else if (s == "--bind-threshold") { auto v = next("--bind-threshold"); if (!v) return false; a.bind_threshold = std::atof(v); }
        else if (s == "--cluster-weight") { auto v = next("--cluster-weight"); if (!v) return false; a.cluster_weight = std::atof(v); }
        else if (s == "--min-margin") { auto v = next("--min-margin"); if (!v) return false; a.min_margin = std::atof(v); }
        else { std::fprintf(stderr, "taxbin: unknown option '%s'\n", s.c_str()); return false; }
    }
    return true;
}

// Parse the confidence matrix. Returns false on error.
bool LoadConfMatrix(const std::string& path, classify::TaxbinInput& in) {
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "taxbin: cannot open --conf %s\n", path.c_str()); return false; }
    std::string line;
    if (!std::getline(f, line)) { std::fprintf(stderr, "taxbin: empty --conf\n"); return false; }
    {  // header
        std::istringstream ss(line);
        std::string col; bool first = true;
        while (std::getline(ss, col, '\t')) {
            if (first) { first = false; continue; }  // skip read_id column header
            in.species_labels.push_back(col);
        }
    }
    const std::size_t S = in.species_labels.size();
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok; bool first = true;
        std::vector<float> row; row.reserve(S);
        std::string rid;
        while (std::getline(ss, tok, '\t')) {
            if (first) { rid = tok; first = false; continue; }
            row.push_back(static_cast<float>(std::atof(tok.c_str())));
        }
        row.resize(S, 0.0f);
        in.read_ids.push_back(std::move(rid));
        in.per_read_species_conf.push_back(std::move(row));
    }
    return true;
}

bool LoadClusters(const std::string& path, classify::TaxbinInput& in) {
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "taxbin: cannot open --clusters %s\n", path.c_str()); return false; }
    std::unordered_map<std::string, std::uint32_t> m;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string rid, cid;
        if (!std::getline(ss, rid, '\t')) continue;
        if (!std::getline(ss, cid, '\t')) continue;
        m[rid] = static_cast<std::uint32_t>(std::strtoul(cid.c_str(), nullptr, 10));
    }
    in.cluster_ids.resize(in.read_ids.size(), 0u);
    for (std::size_t r = 0; r < in.read_ids.size(); ++r) {
        auto it = m.find(in.read_ids[r]);
        in.cluster_ids[r] = (it != m.end()) ? it->second : 0u;
    }
    return true;
}

void WriteReads(const std::string& path, const classify::TaxbinResult& res,
                const classify::TaxbinInput& in) {
    std::ofstream o(path);
    o << "read_id\tcluster_id\tcall\ttop_prob\tmargin\tnovel\tby_cluster";
    for (const auto& sp : in.species_labels) o << "\tL_" << sp;
    o << "\tL_NOVEL\n";
    for (const auto& r : res.reads) {
        const char* call = r.novel ? "NOVEL"
                         : (r.top_species >= 0 ? in.species_labels[r.top_species].c_str() : "NOVEL");
        o << r.read_id << '\t' << r.cluster_id << '\t' << call << '\t'
          << r.top_prob << '\t' << r.margin << '\t' << (r.novel ? 1 : 0) << '\t'
          << (r.by_cluster ? 1 : 0);
        for (float v : r.likelihood) o << '\t' << v;
        o << '\n';
    }
}

void WriteClusters(const std::string& path, const classify::TaxbinResult& res,
                   const classify::TaxbinInput& in) {
    std::ofstream o(path);
    o << "cluster_id\tsize\tconsensus\tpurity\n";
    for (const auto& c : res.clusters) {
        const char* cons = (c.consensus_species >= 0)
            ? in.species_labels[c.consensus_species].c_str() : "NOVEL";
        o << c.cluster_id << '\t' << c.size << '\t' << cons << '\t' << c.purity << '\n';
    }
}

}  // namespace

int run_taxbin(int argc, char** argv) {
    TaxbinArgs a;
    if (!ParseArgs(argc, argv, a)) { PrintUsage(); return 64; }
    if (a.help) { PrintUsage(); return 0; }
    if (a.conf.empty() || a.output.empty()) {
        std::fprintf(stderr, "taxbin: --conf and --output are required\n");
        PrintUsage();
        return 64;
    }

    classify::TaxbinInput in;
    if (!LoadConfMatrix(a.conf, in)) return 1;
    if (!a.clusters.empty() && !LoadClusters(a.clusters, in)) return 1;

    classify::TaxbinConfig cfg;
    cfg.bind_threshold = a.bind_threshold;
    cfg.cluster_weight = a.cluster_weight;
    cfg.min_margin = a.min_margin;
    cfg.enable_cluster_collapse = !a.clusters.empty();

    classify::Taxbin engine(cfg);
    classify::TaxbinResult res = engine.Run(in);

    WriteReads(a.output, res, in);
    if (!a.clusters_out.empty()) WriteClusters(a.clusters_out, res, in);

    std::fprintf(stderr,
        "[taxbin] %zu reads, %zu species, %zu novel/dark, %zu collapsed-by-cluster, %zu clusters\n",
        res.reads.size(), res.num_species, res.num_novel,
        res.num_collapsed_by_cluster, res.clusters.size());
    return 0;
}

}  // namespace llmap::cli
