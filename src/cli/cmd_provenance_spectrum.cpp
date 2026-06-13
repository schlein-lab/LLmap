// LLmap — `llmap provenance-spectrum`: per-read provenance → contamination
// spectrum from any SAM carrying the provenance tags (XB:Z:<class>, XQ:f:<post>,
// XO:i:<bioconfounder-bitmask>). Reads without an XB tag default to Host, so the
// tool is tool-agnostic and the Σ-invariant (every read in exactly one Layer-1
// class) is always checkable.
//
//   <prefix>.contamination_spectrum.tsv  — Layer-1 partition (Σ==N) + Layer-3
//                                          bioconfounder overlay
// Drives the long-/short-read genome acceptance test (which special buckets are
// real in a sample, and does Σ hold).

#include "cli/commands.h"

#include "provenance/contamination_spectrum.h"
#include "provenance/mapping_confusion.h"
#include "provenance/provenance_class.h"
#include "provenance/provenance_resolver.h"
#include "provenance/pseudogene_catalog.h"

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace llmap::cli {

namespace {

struct Args {
    std::string in = "-";
    std::string out_prefix = "provenance";
    bool dup_pos_heuristic = false;
    bool help = false;
};

void PrintUsage() {
    std::fprintf(stderr,
        "Usage: llmap provenance-spectrum [options]\n\n"
        "Aggregate per-read provenance tags from a SAM into the sample's\n"
        "contamination/provenance spectrum (Layer-1 origin partition, Σ==N, plus\n"
        "the Layer-3 bioconfounder overlay). Reads without an XB tag count as\n"
        "Host, so the tool runs on any BAM/SAM.\n\n"
        "Tags read per primary record:\n"
        "  XB:Z:<class>   provenance class (host,exo,xindiv,para,numt,pseudo,\n"
        "                 rdna,mei,refartefact,chim,xsample,dup) [default host]\n"
        "  XQ:f:<post>    collapse posterior [1.0]\n"
        "  XO:i:<mask>    bioconfounder bitmask (Layer-3 overlay) [0]\n\n"
        "Options:\n"
        "  --in FILE        Input SAM (default: stdin)\n"
        "  --out-prefix P   Output prefix [provenance]\n"
        "  --dup-position-heuristic\n"
        "                   Flag coordinate duplicates in UNMARKED BAMs: reads\n"
        "                   sharing the markdup signature (unclipped 5' position +\n"
        "                   strand) as a prior primary read get the `dup` class\n"
        "                   (OR'd with the 0x400 flag). Off by default — only use\n"
        "                   when upstream MarkDuplicates was NOT run.\n"
        "  -h, --help       Show this help\n");
}

bool ParseArgs(int argc, char** argv, Args& a) {
    for (int i = 0; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { a.help = true; return true; }
        else if (arg == "--in" && i + 1 < argc) a.in = argv[++i];
        else if (arg == "--out-prefix" && i + 1 < argc) a.out_prefix = argv[++i];
        else if (arg == "--dup-position-heuristic") a.dup_pos_heuristic = true;
        else { std::fprintf(stderr, "Unknown arg: %s\n", arg.c_str()); return false; }
    }
    return true;
}

// Extract an optional tag value from SAM optional fields (tokens[11..]).
bool FindTag(const std::vector<std::string>& tok, const char* key, std::string& val) {
    const std::string pre = std::string(key);   // e.g. "XB:Z:"
    for (std::size_t i = 11; i < tok.size(); ++i) {
        if (tok[i].rfind(pre, 0) == 0) { val = tok[i].substr(pre.size()); return true; }
    }
    return false;
}

// Derive a read's provenance from plain BAM fields (no XB tag) via the cheap,
// align-free signals available in any sorted BAM. Catalog-backed classes
// (para via PSV, numt via NUMT-BED+chrM, exo via a contaminant panel, mei via a
// TE catalog) need their reference data and correctly stay Host/None here.
// markdup signature of a primary read: refid + 5' UNCLIPPED coordinate in the
// read's orientation + strand. Two reads with the same signature are coordinate
// duplicates — the same key samtools markdup uses for single reads, so we flag
// PCR/optical duplicates even when MarkDuplicates was never run. Empty for
// unmapped reads (no coordinate). CIGAR ops parsed: leading/trailing S/H are the
// clip; M/D/N/=/X consume the reference.
std::string DupSignature(const std::vector<std::string>& tok) {
    const int flag = std::atoi(tok[1].c_str());
    if (flag & 0x4) return {};                         // unmapped → no coordinate
    const std::string& cig = tok[5];
    if (cig == "*" || cig.empty()) return {};
    const std::int64_t pos = std::atoll(tok[3].c_str());
    const bool reverse = (flag & 0x10) != 0;

    std::int64_t ref_span = 0, lead_clip = 0, trail_clip = 0;
    std::int64_t num = 0;
    bool seen_aligned = false, only_clip_so_far = true;
    std::int64_t last_clip = 0;
    for (char c : cig) {
        if (std::isdigit(static_cast<unsigned char>(c))) { num = num * 10 + (c - '0'); continue; }
        switch (c) {
            case 'S': case 'H':
                if (!seen_aligned && only_clip_so_far) lead_clip = num;
                last_clip = num;                       // trailing clip = last clip op seen
                break;
            case 'M': case 'D': case 'N': case '=': case 'X':
                ref_span += num; seen_aligned = true; only_clip_so_far = false; last_clip = 0;
                break;
            case 'I': seen_aligned = true; only_clip_so_far = false; last_clip = 0; break;
            default: break;
        }
        num = 0;
    }
    trail_clip = last_clip;
    // 5' unclipped coordinate in read orientation.
    const std::int64_t coord = reverse ? (pos + ref_span - 1 + trail_clip)
                                       : (pos - lead_clip);
    return tok[2] + (reverse ? "-" : "+") + std::to_string(coord);
}

provenance::ReadProvenance DeriveFromBam(const std::vector<std::string>& tok,
                                         const provenance::PseudogeneCatalog& cat,
                                         bool position_dup = false) {
    const int flag = std::atoi(tok[1].c_str());
    const int mapq = std::atoi(tok[4].c_str());
    const bool unmapped = (flag & 0x4) != 0;
    std::string dummy;

    provenance::ReadEvidence ev;
    ev.aligned_bases = static_cast<std::uint32_t>(tok[9].size());
    ev.host_posterior = unmapped ? 0.0f : std::min(1.0f, static_cast<float>(mapq) / 60.0f);
    // PCR/optical duplicate: the 0x400 flag (if MarkDuplicates ran) OR the
    // coordinate-signature heuristic (for unmarked BAMs, when enabled).
    ev.is_duplicate = ((flag & 0x400) != 0) || position_dup;
    ev.is_chimera = FindTag(tok, "SA:Z:", dummy);     // supplementary mapping (crude)

    provenance::MappingConfusionEvidence me;
    me.mapq = static_cast<std::uint32_t>(mapq);
    me.in_repeat_array = (!unmapped && mapq < 5);     // low-MAPQ pileup ⇒ repeat/rdna
    me.at_pseudogene_parent_locus =
        cat.Lookup(tok[2], static_cast<std::uint64_t>(std::atoll(tok[3].c_str()))) != nullptr;
    me.read_is_spliced = tok[5].find('N') != std::string::npos;
    ev.mapping = provenance::ClassifyMappingConfusion(me);

    return provenance::ResolveProvenance(ev);
}

}  // namespace

int run_provenance_spectrum(int argc, char** argv) {
    Args a;
    if (!ParseArgs(argc, argv, a)) { PrintUsage(); return 2; }
    if (a.help) { PrintUsage(); return 0; }

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (a.in != "-") {
        fin.open(a.in);
        if (!fin) { std::fprintf(stderr, "cannot open %s\n", a.in.c_str()); return 1; }
        in = &fin;
    }

    // Pseudogene parent/pseudogene pair catalog (built-in starter; the genome
    // test can later --gencode/--bed a curated source).
    provenance::PseudogeneCatalog pseudo_cat;
    pseudo_cat.LoadBuiltinStarter();

    provenance::ContaminationSpectrum spectrum;
    std::uint64_t n_input = 0, n_tagged = 0, n_derived = 0, n_pos_dup = 0;
    std::unordered_set<std::string> seen_sigs;   // markdup signatures (heuristic)
    std::string line;
    while (std::getline(*in, line)) {
        if (line.empty() || line[0] == '@') continue;
        std::vector<std::string> tok;
        std::stringstream ss(line);
        std::string f;
        while (std::getline(ss, f, '\t')) tok.push_back(f);
        if (tok.size() < 11) continue;
        const int flag = std::atoi(tok[1].c_str());
        if (flag & 0x100 || flag & 0x800) continue;   // skip secondary/supplementary
        ++n_input;

        provenance::ReadProvenance rp;
        std::string v;
        if (FindTag(tok, "XB:Z:", v)) {
            // Pre-tagged path (tool-agnostic): trust the upstream provenance tags.
            ++n_tagged;
            if (auto c = provenance::ParseProvenanceClass(v)) rp.origin = *c;
            if (FindTag(tok, "XQ:f:", v)) rp.posterior = std::strtof(v.c_str(), nullptr);
            if (FindTag(tok, "XO:i:", v))
                rp.bioconfounder =
                    static_cast<std::uint16_t>(std::strtoul(v.c_str(), nullptr, 10));
            rp.aligned_bases = static_cast<std::uint32_t>(tok[9].size());
        } else {
            // Derive from plain BAM fields (the genome-test path).
            ++n_derived;
            bool position_dup = false;
            if (a.dup_pos_heuristic) {
                const std::string sig = DupSignature(tok);
                // First read at a signature is the original; later ones are dups.
                if (!sig.empty() && !seen_sigs.insert(sig).second) {
                    position_dup = true;
                    ++n_pos_dup;
                }
            }
            rp = DeriveFromBam(tok, pseudo_cat, position_dup);
        }
        spectrum.Add(rp);
    }
    std::fprintf(stderr, "provenance-spectrum: %llu reads (%llu pre-tagged, "
                 "%llu derived-from-BAM)\n",
                 static_cast<unsigned long long>(n_input),
                 static_cast<unsigned long long>(n_tagged),
                 static_cast<unsigned long long>(n_derived));
    if (a.dup_pos_heuristic)
        std::fprintf(stderr, "provenance-spectrum: %llu coordinate-duplicates flagged "
                     "by position heuristic (unmarked-BAM mode)\n",
                     static_cast<unsigned long long>(n_pos_dup));

    const std::string out = a.out_prefix + ".contamination_spectrum.tsv";
    if (!spectrum.WriteTsv(out)) {
        std::fprintf(stderr, "cannot write %s\n", out.c_str());
        return 1;
    }
    const bool lossless = spectrum.CheckLossless(n_input);
    std::fprintf(stderr, "provenance-spectrum: %llu reads → %s  (Σ-invariant: %s)\n",
                 static_cast<unsigned long long>(n_input), out.c_str(),
                 lossless ? "OK" : "VIOLATED");
    std::fprintf(stderr, "%s", spectrum.ToString().c_str());
    return lossless ? 0 : 3;
}

}  // namespace llmap::cli
