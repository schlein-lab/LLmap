// LLmap -- `llmap annotate-specific-loci`.
//
// Walks knowledge/organisms/<org>/specific_loci/**.json and emits a
// .regann fragment containing Layer-2 intervals for every locus whose
// coordinates fall on a contig present in the supplied reference.
//
// The resulting file can be concatenated with the Layer-1 output of
// `llmap annotate-ref` to produce a combined Layer-1+2 annotation.

#include "annot/annotation_store.h"
#include "annot/annot_types.h"
#include "cli/commands.h"
#include "io/fasta_reader.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace llmap::cli {

namespace {

struct LociArgs {
    std::string reference;
    std::string loci_dir;
    std::string output;
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap annotate-specific-loci \\\n"
        "         --reference REF.fa \\\n"
        "         --loci-dir knowledge/organisms/human/specific_loci \\\n"
        "         --output REF.specific_loci.regann\n"
        "\n"
        "Walks the loci-dir recursively, parses every *.json with a\n"
        "coordinates.grch38 entry, and emits a Layer-2 (.regann) fragment\n"
        "containing one interval per locus whose chr is present in the\n"
        "reference. Output is the standard .regann text format.\n");
}

bool ParseArgs(int argc, char** argv, LociArgs& a) {
    for (int i = 0; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--reference" && i + 1 < argc) a.reference = argv[++i];
        else if (arg == "--loci-dir"  && i + 1 < argc) a.loci_dir = argv[++i];
        else if (arg == "--output"    && i + 1 < argc) a.output = argv[++i];
        else if (arg == "--help" || arg == "-h") { a.help = true; return true; }
        else {
            std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str());
            return false;
        }
    }
    return true;
}

// Trivial JSON field extractor -- finds "key" : <number-or-string-value>
// and returns it as a string. Returns empty if not found.
std::string ExtractField(const std::string& s, const std::string& key) {
    auto k = "\"" + key + "\"";
    auto p = s.find(k);
    if (p == std::string::npos) return {};
    p = s.find(':', p);
    if (p == std::string::npos) return {};
    ++p;
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p;
    if (p >= s.size()) return {};
    if (s[p] == '"') {
        ++p;
        auto end = s.find('"', p);
        if (end == std::string::npos) return {};
        return s.substr(p, end - p);
    } else {
        auto end = p;
        while (end < s.size() && (std::isdigit(static_cast<unsigned char>(s[end])) ||
                                   s[end] == '.' || s[end] == '-' ||
                                   s[end] == 'e' || s[end] == 'E')) ++end;
        return s.substr(p, end - p);
    }
}

// Find a {chr, start, end} block under "coordinates" -> "grch38"
struct LocusCoords {
    std::string chr;
    uint64_t start = 0;
    uint64_t end = 0;
    bool valid = false;
};

LocusCoords ExtractGrch38Coords(const std::string& s) {
    LocusCoords c;
    auto pos = s.find("\"coordinates\"");
    if (pos == std::string::npos) return c;
    auto g = s.find("\"grch38\"", pos);
    if (g == std::string::npos) return c;
    auto block_end = s.find('}', g);
    if (block_end == std::string::npos) return c;
    auto sub = s.substr(g, block_end - g);
    c.chr = ExtractField(sub, "chr");
    auto s_str = ExtractField(sub, "start");
    auto e_str = ExtractField(sub, "end");
    if (c.chr.empty() || s_str.empty() || e_str.empty()) return c;
    try {
        c.start = std::stoull(s_str);
        c.end   = std::stoull(e_str);
        c.valid = c.end > c.start;
    } catch (...) {}
    return c;
}

}  // namespace

int run_annotate_specific_loci(int argc, char** argv) {
    LociArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 1; }
    if (args.help) { PrintUsage(); return 0; }
    if (args.reference.empty() || args.loci_dir.empty() || args.output.empty()) {
        std::fprintf(stderr, "Missing required arg.\n");
        PrintUsage();
        return 1;
    }

    // Load contig names from reference (we only need names, not full seqs)
    io::FastaReader rdr(args.reference);
    std::unordered_map<std::string, uint32_t> contig_to_id;
    std::vector<std::string> contig_names;
    while (rdr.HasMore()) {
        auto rec = rdr.Next();
        if (rec.IsValid()) {
            contig_to_id[rec.name] = static_cast<uint32_t>(contig_names.size());
            contig_names.push_back(rec.name);
        }
    }

    // Also accept "chr14" -> "chr14:105500000-107300000"-style names by
    // prefix-matching: if the reference contig name starts with "chrX",
    // any locus chr "chrX" matches.
    auto match_contig = [&](const std::string& chr) -> int32_t {
        if (auto it = contig_to_id.find(chr); it != contig_to_id.end())
            return static_cast<int32_t>(it->second);
        // Prefix match: reference contig might be "chr14:..."
        for (size_t i = 0; i < contig_names.size(); ++i) {
            if (contig_names[i].rfind(chr + ":", 0) == 0 ||
                contig_names[i] == chr) {
                return static_cast<int32_t>(i);
            }
        }
        return -1;
    };

    std::ofstream out(args.output);
    if (!out) {
        std::fprintf(stderr, "failed to open %s\n", args.output.c_str());
        return 1;
    }
    out << "# LLmap region annotation -- Layer 2 specific_loci fragment\n";
    out << "# ref_name\tstart\tend\tregion_name\tsource\tlayer\tparams\n";

    size_t scanned = 0, matched = 0, skipped_no_match = 0, skipped_invalid = 0;

    for (auto& p : std::filesystem::recursive_directory_iterator(args.loci_dir)) {
        if (!p.is_regular_file()) continue;
        if (p.path().extension() != ".json") continue;
        ++scanned;

        std::ifstream f(p.path());
        std::stringstream ss;
        ss << f.rdbuf();
        std::string s = ss.str();

        auto coords = ExtractGrch38Coords(s);
        if (!coords.valid) { ++skipped_invalid; continue; }

        int32_t contig_idx = match_contig(coords.chr);
        if (contig_idx < 0) { ++skipped_no_match; continue; }

        std::string name = ExtractField(s, "name");
        std::string tax = ExtractField(s, "taxonomy_id");
        if (name.empty()) name = p.path().stem().string();

        // Build a minimal ParamOverride string with sensible defaults per
        // taxonomy_id. Reads mapping_hints fields directly via ExtractField.
        std::string params;
        auto add = [&](const std::string& k, const std::string& v) {
            if (v.empty()) return;
            if (!params.empty()) params += ',';
            params += k + "=" + v;
        };
        add("lambda_scale", ExtractField(s, "lambda_scale"));
        add("anchor_weight_scale", ExtractField(s, "anchor_weight_scale"));
        add("max_occ", ExtractField(s, "max_occ"));
        std::string rmp = ExtractField(s, "report_multi_position");
        if (!rmp.empty()) add("report_multi_position",
                              (rmp == "true" || rmp == "1") ? "1" : "0");
        std::string rps = ExtractField(s, "require_psv_disambiguation");
        if (!rps.empty()) add("require_psv_disambig",
                              (rps == "true" || rps == "1") ? "1" : "0");
        std::string ahm = ExtractField(s, "allow_high_mismatch");
        if (!ahm.empty()) add("allow_high_mismatch",
                              (ahm == "true" || ahm == "1") ? "1" : "0");
        if (params.empty()) params = "-";

        out << contig_names[contig_idx] << '\t'
            << coords.start << '\t'
            << coords.end << '\t'
            << name << '\t'
            << "specific_locus:" << name << '\t'
            << "2" << '\t'
            << params << '\n';
        ++matched;
    }

    std::printf("Scanned: %zu files\n", scanned);
    std::printf("Matched to reference contigs: %zu\n", matched);
    std::printf("Skipped (no contig match): %zu\n", skipped_no_match);
    std::printf("Skipped (no grch38 coords): %zu\n", skipped_invalid);
    std::printf("Output: %s\n", args.output.c_str());
    return 0;
}

}  // namespace llmap::cli
