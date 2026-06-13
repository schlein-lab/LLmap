// LLmap — `llmap splice-determinism`: per-position splice determinism D(pos)
// + junction usage from any spliced SAM (llmap OR minimap2 -ax splice).
//
// Reads a SAM (stdin or --in), accumulates, and writes:
//   <prefix>.determinism.bedgraph  — per-base D(pos) in [0,100] (runs merged)
//   <prefix>.junctions.tsv         — ref, donor, acceptor, n_reads
//
// Tool-agnostic by design: the same metric quantifies isoform blurriness for
// llmap and minimap2 output, and (with GENCODE) feeds the boundary-concordance
// test. D(pos) = 100 * modal-state-fraction = 1 - normalised splice entropy.

#include "cli/commands.h"

#include "mapping/splice_determinism.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

struct Args {
    std::string in = "-";
    std::string out_prefix = "splice_determinism";
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap splice-determinism [options]\n\n"
        "Per-position splice determinism D(pos) + junction usage from a spliced\n"
        "SAM (llmap or minimap2 -ax splice). D(pos) = 100 * modal-state fraction\n"
        "= 1 - normalised splice entropy; 100 = constitutive, lower = alternative\n"
        "splicing / blurry boundaries.\n\n"
        "Options:\n"
        "  --in FILE        Input SAM (default: stdin)\n"
        "  --out-prefix P   Output prefix [splice_determinism]\n"
        "                   writes P.determinism.bedgraph + P.junctions.tsv\n"
        "  -h, --help       Show this help\n");
}

bool ParseArgs(int argc, char** argv, Args& a) {
    for (int i = 0; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            a.help = true;
            return true;
        } else if (arg == "--in" && i + 1 < argc) {
            a.in = argv[++i];
        } else if (arg == "--out-prefix" && i + 1 < argc) {
            a.out_prefix = argv[++i];
        } else if (!arg.empty() && arg[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

std::vector<std::string> SplitTab(const std::string& s) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t tab = s.find('\t', start);
        if (tab == std::string::npos) {
            out.push_back(s.substr(start));
            break;
        }
        out.push_back(s.substr(start, tab - start));
        start = tab + 1;
    }
    return out;
}

}  // namespace

int run_splice_determinism(int argc, char** argv) {
    Args a;
    if (!ParseArgs(argc, argv, a)) {
        PrintUsage();
        return 64;
    }
    if (a.help) {
        PrintUsage();
        return 0;
    }

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (a.in != "-") {
        fin.open(a.in);
        if (!fin) {
            std::fprintf(stderr, "error: cannot open %s\n", a.in.c_str());
            return 66;
        }
        in = &fin;
    }

    // Buffer the spliced alignments (primary + supplementary; not secondary /
    // unmapped) — two passes are needed: mark exon-body positions included,
    // then mark intron-skipped positions excluded over the exonic set.
    struct Rec {
        std::string ref;
        std::uint64_t pos0;
        std::string cigar;
    };
    std::vector<Rec> recs;
    std::size_t n_aln = 0;
    std::string line;
    while (std::getline(*in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '@') continue;
        const std::vector<std::string> f = SplitTab(line);
        if (f.size() < 6) continue;
        const int flag = std::atoi(f[1].c_str());
        if ((flag & 0x4) || (flag & 0x100)) continue;  // unmapped / secondary
        if (f[2] == "*" || f[5] == "*") continue;
        const long long pos1 = std::atoll(f[3].c_str());
        if (pos1 <= 0) continue;
        recs.push_back(Rec{f[2], static_cast<std::uint64_t>(pos1 - 1), f[5]});
        ++n_aln;
    }

    mapping::DeterminismAccumulator acc;
    for (const auto& r : recs) acc.MarkIncluded(r.ref, r.pos0, r.cigar);
    for (const auto& r : recs) acc.MarkExcluded(r.ref, r.pos0, r.cigar);

    // BedGraph: merge runs of contiguous positions with equal D.
    const std::string bg_path = a.out_prefix + ".determinism.bedgraph";
    std::ofstream bg(bg_path);
    if (!bg) {
        std::fprintf(stderr, "error: cannot write %s\n", bg_path.c_str());
        return 66;
    }
    std::size_t n_positions = 0;
    for (const auto& [ref, pmap] : acc.positions) {
        bool open = false;
        std::uint64_t run_start = 0, prev = 0;
        double run_d = -1.0;
        for (const auto& [pos, counts] : pmap) {
            const double d = mapping::Determinism(counts);
            ++n_positions;
            const bool contiguous = open && pos == prev + 1 && d == run_d;
            if (!contiguous) {
                if (open) {
                    bg << ref << '\t' << run_start << '\t' << (prev + 1) << '\t'
                       << run_d << '\n';
                }
                run_start = pos;
                run_d = d;
                open = true;
            }
            prev = pos;
        }
        if (open) {
            bg << ref << '\t' << run_start << '\t' << (prev + 1) << '\t' << run_d
               << '\n';
        }
    }
    bg.close();

    // Junctions TSV.
    const std::string jp = a.out_prefix + ".junctions.tsv";
    std::ofstream jt(jp);
    if (!jt) {
        std::fprintf(stderr, "error: cannot write %s\n", jp.c_str());
        return 66;
    }
    jt << "ref\tdonor\tacceptor\tn_reads\n";
    std::size_t n_junctions = 0;
    for (const auto& [ref, jmap] : acc.junctions) {
        for (const auto& usage : mapping::SortedJunctions(acc, ref)) {
            jt << ref << '\t' << usage.donor << '\t' << usage.acceptor << '\t'
               << usage.n_reads << '\n';
            ++n_junctions;
        }
    }
    jt.close();

    std::fprintf(stderr,
        "[splice-determinism] alignments=%zu positions=%zu junctions=%zu\n"
        "  wrote %s, %s\n",
        n_aln, n_positions, n_junctions, bg_path.c_str(), jp.c_str());
    return 0;
}

}  // namespace llmap::cli
