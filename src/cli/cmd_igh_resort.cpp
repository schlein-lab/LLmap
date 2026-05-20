// LLmap — standalone `igh-resort` command.
//
// Re-sorts an EXISTING alignment (from minimap2 or any mapper) within the IGH
// constant region, using paralog-specific exon anchors. Works on the SAM stream
// so it composes with samtools and needs no htslib:
//
//   samtools view -h in.bam
//     | llmap igh-resort --anchors igh_anchors.fa
//     | samtools sort -o resorted.bam
//
// For each primary alignment whose target is an IGH copy (transcriptomic) or
// whose position is inside the IGHC locus (genomic), the read's SEQ is scanned
// against the anchors. When the exact-anchor evidence names a different copy,
// RNAME/POS are re-pointed and provenance tags are added:
//   ir:Z:<original RNAME>   ip:Z:<paralog gene>   ie:i:<distinct exon count>
//
// Reads minimap2 already dropped (absent from the SAM) cannot be recovered, but
// re-sorting the survivors is still better than the collapsed call.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cli/commands.h"
#include "igh_locus/igh_anchor_catalog.h"
#include "igh_locus/igh_region.h"

namespace llmap::cli {

namespace {

struct Opts {
    std::string anchors;
    std::string in_path = "-";
    std::string out_path = "-";
    int max_mismatch = 0;
    int min_exons = 2;
    bool verbose = false;
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap igh-resort --anchors FILE [options] < in.sam > out.sam\n\n"
        "Re-sort IGH-locus reads of an existing SAM/BAM to their true paralog\n"
        "copy using exact exon anchors. Reads SAM on stdin, writes SAM on stdout.\n\n"
        "Required:\n"
        "  --anchors FILE        FASTA of paralog-specific CH exon anchors\n"
        "                        headers: GENE_COPY_HAP_EXON [ loc=chr:start-end ]\n\n"
        "Options:\n"
        "  --in FILE             Input SAM (default: stdin)\n"
        "  --out FILE            Output SAM (default: stdout)\n"
        "  --max-mismatch INT    Tolerated mismatches per exon anchor [0 = exact]\n"
        "  --min-exons INT       Distinct exons required to re-sort [2]\n"
        "  -v, --verbose         Per-read re-sort log to stderr\n"
        "  -h, --help            Show this help\n\n"
        "Example:\n"
        "  samtools view -h in.bam | llmap igh-resort --anchors igh.fa \\\n"
        "    | samtools sort -o resorted.bam\n");
}

bool ParseOpts(int argc, char** argv, Opts& o) {
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { o.help = true; return true; }
        else if (a == "--anchors" && i + 1 < argc) o.anchors = argv[++i];
        else if (a == "--in" && i + 1 < argc) o.in_path = argv[++i];
        else if (a == "--out" && i + 1 < argc) o.out_path = argv[++i];
        else if (a == "--max-mismatch" && i + 1 < argc) o.max_mismatch = std::atoi(argv[++i]);
        else if (a == "--min-exons" && i + 1 < argc) o.min_exons = std::atoi(argv[++i]);
        else if (a == "-v" || a == "--verbose") o.verbose = true;
        else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "Unknown option: %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

std::vector<std::string> SplitTab(const std::string& line) {
    std::vector<std::string> f;
    std::size_t start = 0;
    while (true) {
        std::size_t tab = line.find('\t', start);
        if (tab == std::string::npos) { f.push_back(line.substr(start)); break; }
        f.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return f;
}

}  // namespace

int run_igh_resort(int argc, char** argv) {
    Opts o;
    if (!ParseOpts(argc, argv, o)) return 64;
    if (o.help) { PrintUsage(); return 0; }
    if (o.anchors.empty()) {
        std::fprintf(stderr, "error: --anchors is required\n\n");
        PrintUsage();
        return 64;
    }

    auto catalog = igh_locus::IghAnchorCatalog::LoadFasta(o.anchors, o.max_mismatch);
    if (!catalog) {
        std::fprintf(stderr, "error: could not load IGH anchors from %s\n",
                     o.anchors.c_str());
        return 65;
    }
    const igh_locus::IghRegion region = igh_locus::IghRegion::Default();
    std::fprintf(stderr, "[igh-resort] %zu anchors loaded (max_mismatch=%d)\n",
                 catalog->size(), o.max_mismatch);

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (o.in_path != "-") {
        fin.open(o.in_path);
        if (!fin) { std::fprintf(stderr, "error: cannot open %s\n", o.in_path.c_str()); return 66; }
        in = &fin;
    }
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (o.out_path != "-") {
        fout.open(o.out_path);
        if (!fout) { std::fprintf(stderr, "error: cannot open %s\n", o.out_path.c_str()); return 66; }
        out = &fout;
    }

    std::size_t n_aln = 0, n_in_locus = 0, n_resorted = 0;
    bool pg_emitted = false;
    std::string line;
    while (std::getline(*in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (!line.empty() && line[0] == '@') {
            (*out) << line << '\n';
            if (line.rfind("@PG", 0) == 0) {
                // remember we saw a PG; we append ours after the last header line
            }
            continue;
        }
        if (!pg_emitted) {
            (*out) << "@PG\tID:llmap-igh-resort\tPN:llmap\tVN:1.0.0\n";
            pg_emitted = true;
        }
        if (line.empty()) continue;

        std::vector<std::string> f = SplitTab(line);
        if (f.size() < 11) { (*out) << line << '\n'; continue; }
        ++n_aln;

        const int flag = std::atoi(f[1].c_str());
        const bool unmapped = (flag & 0x4) != 0;
        const bool secondary = (flag & 0x100) != 0;
        const bool supplementary = (flag & 0x800) != 0;
        std::string& rname = f[2];
        std::string& pos = f[3];
        const std::string& seq = f[9];

        if (unmapped || secondary || supplementary || rname == "*" || seq == "*") {
            (*out) << line << '\n';
            continue;
        }

        const auto cur_gene = catalog->GeneOfTarget(rname);
        const std::uint64_t pos0 =
            pos.empty() ? 0 : static_cast<std::uint64_t>(std::strtoull(pos.c_str(), nullptr, 10));
        const bool genomic = !cur_gene.has_value() &&
                             region.Contains(rname, pos0 == 0 ? 0 : pos0 - 1);
        if (!cur_gene.has_value() && !genomic) { (*out) << line << '\n'; continue; }
        ++n_in_locus;

        igh_locus::IghMatch m = catalog->Match(seq);
        if (!m.matched || m.ambiguous_gene || m.n_distinct_exons < o.min_exons) {
            (*out) << line << '\n';
            continue;
        }

        const std::string original_rname = rname;
        bool changed = false;
        if (cur_gene.has_value()) {
            if (rname != m.copy_label) { rname = m.copy_label; changed = true; }
        } else if (m.contig.has_value() && m.start.has_value()) {
            const std::uint64_t cs = *m.start;          // 1-based SAM POS
            const std::uint64_t ce = m.end.value_or(cs);
            const bool inside = (rname == *m.contig) && pos0 >= cs && pos0 < ce;
            if (!inside) {
                rname = *m.contig;
                pos = std::to_string(cs);
                changed = true;
            }
        }

        // Rebuild line with provenance tags appended.
        std::ostringstream os;
        for (std::size_t i = 0; i < f.size(); ++i) {
            if (i) os << '\t';
            os << f[i];
        }
        os << "\tip:Z:" << m.gene << "\tie:i:" << m.n_distinct_exons;
        if (changed) {
            os << "\tir:Z:" << original_rname;
            ++n_resorted;
            if (o.verbose) {
                std::fprintf(stderr, "[igh-resort] %s : %s -> %s (%s, %d exons)\n",
                             f[0].c_str(), original_rname.c_str(),
                             rname.c_str(), m.gene.c_str(), m.n_distinct_exons);
            }
        }
        (*out) << os.str() << '\n';
    }

    out->flush();
    std::fprintf(stderr,
        "[igh-resort] alignments=%zu in_locus=%zu resorted=%zu\n",
        n_aln, n_in_locus, n_resorted);
    return 0;
}

}  // namespace llmap::cli
