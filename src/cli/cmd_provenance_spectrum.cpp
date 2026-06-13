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

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

struct Args {
    std::string in = "-";
    std::string out_prefix = "provenance";
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
        "  -h, --help       Show this help\n");
}

bool ParseArgs(int argc, char** argv, Args& a) {
    for (int i = 0; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { a.help = true; return true; }
        else if (arg == "--in" && i + 1 < argc) a.in = argv[++i];
        else if (arg == "--out-prefix" && i + 1 < argc) a.out_prefix = argv[++i];
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
provenance::ReadProvenance DeriveFromBam(const std::vector<std::string>& tok,
                                         const provenance::PseudogeneCatalog& cat) {
    const int flag = std::atoi(tok[1].c_str());
    const int mapq = std::atoi(tok[4].c_str());
    const bool unmapped = (flag & 0x4) != 0;
    std::string dummy;

    provenance::ReadEvidence ev;
    ev.aligned_bases = static_cast<std::uint32_t>(tok[9].size());
    ev.host_posterior = unmapped ? 0.0f : std::min(1.0f, static_cast<float>(mapq) / 60.0f);
    ev.is_duplicate = (flag & 0x400) != 0;            // PCR/optical duplicate
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
    std::uint64_t n_input = 0, n_tagged = 0, n_derived = 0;
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
            rp = DeriveFromBam(tok, pseudo_cat);
        }
        spectrum.Add(rp);
    }
    std::fprintf(stderr, "provenance-spectrum: %llu reads (%llu pre-tagged, "
                 "%llu derived-from-BAM)\n",
                 static_cast<unsigned long long>(n_input),
                 static_cast<unsigned long long>(n_tagged),
                 static_cast<unsigned long long>(n_derived));

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
