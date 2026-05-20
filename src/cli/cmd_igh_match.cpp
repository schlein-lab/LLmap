// LLmap — standalone `igh-match` command.
//
// Alignment-free IGH paralog assignment straight from reads (FASTQ/FASTA), with
// no prior mapping. Each read is scanned against paralog-specific CH exon anchors
// and assigned to the copy with the most distinct exact-exon matches. This is the
// alignment-free principle used to recover IGH paralog reads that minimap2
// mis-assigns, exposed as a first-class command so LLmap can do the per-read copy
// call directly (the `igh-resort` command does the same Match() but needs an
// already-aligned SAM):
//
//   samtools fastq sample.flnc.bam | llmap igh-match --anchors per_genome.fa
//
// Per-read TSV on stdout: read_id, gene, copy, n_exons, strand, ambiguous.
// Per-gene and per-copy tallies are written to stderr.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "cli/commands.h"
#include "core/remote_fetch.h"
#include "igh_locus/igh_anchor_catalog.h"

namespace llmap::cli {

namespace {

struct Opts {
    std::string anchors;
    std::string in_path = "-";
    std::string out_path = "-";
    int max_mismatch = 0;
    int min_exons = 1;
    bool verbose = false;
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap igh-match --anchors FILE [options] < reads.fastq > calls.tsv\n\n"
        "Alignment-free IGH paralog assignment from reads (no prior mapping).\n"
        "Reads FASTQ/FASTA on stdin, writes a per-read paralog-call TSV on stdout.\n\n"
        "Required:\n"
        "  --anchors FILE|URL    FASTA of paralog-specific CH exon anchors\n"
        "                        headers: GENE_COPY_HAP_EXON [ loc=chr:start-end ]\n"
        "                        accepts local path, http(s):// or s3:// URL\n\n"
        "Options:\n"
        "  --in FILE|URL         Input reads (default: stdin; accepts URL)\n"
        "  --out FILE            Output TSV (default: stdout)\n"
        "  --max-mismatch INT    Tolerated mismatches per exon anchor [0 = exact]\n"
        "  --min-exons INT       Distinct exons required for a call [1]\n"
        "  -v, --verbose         Progress to stderr every 1M reads\n"
        "  -h, --help            Show this help\n\n"
        "Example:\n"
        "  samtools fastq -n in.flnc.bam | llmap igh-match --anchors per_genome.fa \\\n"
        "    > calls.tsv\n");
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

// Trim a read header to the first whitespace-delimited token (drops the '@'/'>').
std::string ReadId(const std::string& header) {
    std::size_t b = (!header.empty() && (header[0] == '@' || header[0] == '>')) ? 1 : 0;
    std::size_t e = header.find_first_of(" \t", b);
    return header.substr(b, e == std::string::npos ? std::string::npos : e - b);
}

}  // namespace

int run_igh_match(int argc, char** argv) {
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
    std::fprintf(stderr, "[igh-match] %zu anchors loaded (max_mismatch=%d)\n",
                 catalog->size(), o.max_mismatch);

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (o.in_path != "-") {
        std::string in_local = o.in_path;
        if (core::IsRemotePath(o.in_path)) {
            auto fetched = core::FetchToCache(o.in_path);
            if (!fetched) {
                std::fprintf(stderr, "error: could not fetch %s\n", o.in_path.c_str());
                return 66;
            }
            in_local = *fetched;
        }
        fin.open(in_local);
        if (!fin) { std::fprintf(stderr, "error: cannot open %s\n", in_local.c_str()); return 66; }
        in = &fin;
    }
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (o.out_path != "-") {
        fout.open(o.out_path);
        if (!fout) { std::fprintf(stderr, "error: cannot open %s\n", o.out_path.c_str()); return 66; }
        out = &fout;
    }

    // Detect FASTQ ('@') vs FASTA ('>') from the first non-empty line.
    std::string first;
    while (std::getline(*in, first)) {
        if (!first.empty() && first.back() == '\r') first.pop_back();
        if (!first.empty()) break;
    }
    const bool fastq = !first.empty() && first[0] == '@';
    const bool fasta = !first.empty() && first[0] == '>';
    if (!fastq && !fasta) {
        std::fprintf(stderr, "error: input is neither FASTQ nor FASTA\n");
        return 65;
    }

    (*out) << "read_id\tgene\tcopy\tn_exons\tstrand\tambiguous_gene\t"
              "ambiguous_copy\ttied_copies\n";

    std::size_t n_reads = 0, n_called = 0, n_ambiguous = 0, n_ambiguous_copy = 0;
    std::map<std::string, std::size_t> per_gene, per_copy;

    auto handle = [&](const std::string& header, const std::string& seq) {
        ++n_reads;
        if (o.verbose && n_reads % 1000000 == 0)
            std::fprintf(stderr, "[igh-match] %zu reads, %zu calls\n", n_reads, n_called);
        igh_locus::IghMatch m = catalog->Match(seq);
        if (!m.matched || m.n_distinct_exons < o.min_exons) return;
        ++n_called;
        if (m.ambiguous_gene) ++n_ambiguous;
        if (m.ambiguous_copy) ++n_ambiguous_copy;
        per_gene[m.gene]++;
        per_copy[m.copy_label]++;
        std::string tied;
        for (const auto& c : m.tied_copies) {
            if (!tied.empty()) tied += ",";
            tied += c;
        }
        (*out) << ReadId(header) << '\t' << m.gene << '\t' << m.copy_label << '\t'
               << m.n_distinct_exons << '\t' << (m.is_reverse ? '-' : '+') << '\t'
               << (m.ambiguous_gene ? "1" : "0") << '\t'
               << (m.ambiguous_copy ? "1" : "0") << '\t' << tied << '\n';
    };

    std::string line, seq;
    if (fastq) {
        std::string header = first;  // first line is the first record header
        while (true) {
            if (!std::getline(*in, seq)) break;
            if (!seq.empty() && seq.back() == '\r') seq.pop_back();
            if (!std::getline(*in, line)) break;  // '+'
            if (!std::getline(*in, line)) break;  // qual
            handle(header, seq);
            if (!std::getline(*in, header)) break;
            if (!header.empty() && header.back() == '\r') header.pop_back();
            if (header.empty()) break;
        }
    } else {
        std::string header = first;
        seq.clear();
        while (std::getline(*in, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty() && line[0] == '>') {
                handle(header, seq);
                header = line;
                seq.clear();
            } else {
                seq += line;
            }
        }
        handle(header, seq);  // last record
    }

    out->flush();
    std::fprintf(stderr,
        "[igh-match] reads=%zu calls=%zu ambiguous_gene=%zu ambiguous_copy=%zu\n",
        n_reads, n_called, n_ambiguous, n_ambiguous_copy);
    std::fprintf(stderr, "[igh-match] per-gene:\n");
    for (const auto& [g, c] : per_gene)
        std::fprintf(stderr, "    %s\t%zu\n", g.c_str(), c);
    std::fprintf(stderr, "[igh-match] per-copy:\n");
    for (const auto& [cp, c] : per_copy)
        std::fprintf(stderr, "    %s\t%zu\n", cp.c_str(), c);
    return 0;
}

}  // namespace llmap::cli
