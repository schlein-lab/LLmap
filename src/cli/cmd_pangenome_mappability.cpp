// LLmap — `llmap pangenome-mappability`: build the cross-sample M(pos) mapping-
// determinism track from raw pangenome BAMs (the Operator's population-prior
// directive — a placement is more probable when reproducibly mappable across many
// pangenome samples). The second free aggregate of the same HPRC run that produced
// the provenance baseline.
//
//   persample:  samtools view BAM region | llmap pangenome-mappability --persample
//               → per-window mean mapping confidence (mean MAPQ→[0,1]) for ONE
//                 sample, as `ref<TAB>window_start<TAB>conf` lines.
//   build:      llmap pangenome-mappability --build OUT.track  S1.win S2.win ...
//               → AddSample each sample's per-window confidence into M(pos) and
//                 Save the track. M = cross-sample mean; n = contributing samples.

#include "cli/commands.h"

#include "core/pangenome_mappability.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

void PrintUsage() {
    std::fprintf(stderr,
        "Usage:\n"
        "  samtools view BAM region | llmap pangenome-mappability --persample [--window 100]\n"
        "      Per-window mean mapping confidence (mean MAPQ/60) for ONE sample's\n"
        "      SAM stream → ref<TAB>window_start<TAB>conf on stdout.\n"
        "  llmap pangenome-mappability --build OUT.track [--window 100] S1.win S2.win ...\n"
        "      Aggregate N per-sample window files into the M(pos) track\n"
        "      (cross-sample mean + n_samples).\n");
}

// One sample's SAM on stdin → per-window mean confidence. confidence = MAPQ/60
// clamped to [0,1] (a per-read mapping-certainty proxy); the window mean is the
// per-sample value the track then averages across samples.
int PerSample(std::uint32_t window) {
    struct Acc { double sum = 0.0; std::uint64_t n = 0; };
    std::map<std::string, std::map<std::uint64_t, Acc>> bins;

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty() || line[0] == '@') continue;
        std::stringstream ss(line);
        std::vector<std::string> t;
        std::string c;
        for (int i = 0; i < 5 && std::getline(ss, c, '\t'); ++i) t.push_back(c);
        if (t.size() < 5) continue;
        const std::uint32_t flag = static_cast<std::uint32_t>(std::strtoul(t[1].c_str(), nullptr, 10));
        if (flag & 0x4) continue;                 // unmapped
        if (t[2] == "*") continue;
        const std::uint64_t pos = std::strtoull(t[3].c_str(), nullptr, 10);
        if (pos == 0) continue;
        const long mapq = std::strtol(t[4].c_str(), nullptr, 10);
        const float conf = std::min(1.0f, std::max(0.0f, static_cast<float>(mapq) / 60.0f));
        const std::uint64_t w = (pos / window) * window;
        auto& a = bins[t[2]][w];
        a.sum += conf;
        ++a.n;
    }

    std::uint64_t emitted = 0;
    for (const auto& [ref, wins] : bins)
        for (const auto& [w, a] : wins) {
            std::printf("%s\t%llu\t%.6f\n", ref.c_str(),
                        static_cast<unsigned long long>(w), a.sum / static_cast<double>(a.n));
            ++emitted;
        }
    std::fprintf(stderr, "persample: %llu windows over %zu refs\n",
                 static_cast<unsigned long long>(emitted), bins.size());
    return 0;
}

int Build(const std::string& out, std::uint32_t window,
          const std::vector<std::string>& files) {
    core::PangenomeMappability track(window);
    std::size_t lines = 0;
    for (const auto& path : files) {
        std::ifstream f(path);
        if (!f) { std::fprintf(stderr, "cannot read %s\n", path.c_str()); continue; }
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::vector<std::string> t;
            std::string c;
            while (std::getline(ss, c, '\t')) t.push_back(c);
            if (t.size() < 3) continue;
            const std::uint64_t pos = std::strtoull(t[1].c_str(), nullptr, 10);
            const float conf = static_cast<float>(std::strtod(t[2].c_str(), nullptr));
            track.AddSample(t[0], pos, conf);   // one value per (sample, window)
            ++lines;
        }
    }
    if (!track.Save(out)) {
        std::fprintf(stderr, "cannot write %s\n", out.c_str());
        return 1;
    }
    std::fprintf(stderr, "M(pos) track: %zu windows from %zu samples (%zu rows) → %s\n",
                 track.n_windows(), files.size(), lines, out.c_str());
    return 0;
}

}  // namespace

int run_pangenome_mappability(int argc, char** argv) {
    bool persample = false;
    std::string build_out;
    std::uint32_t window = 100;
    std::vector<std::string> files;
    for (int i = 0; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") { PrintUsage(); return 0; }
        else if (a == "--persample") persample = true;
        else if (a == "--build" && i + 1 < argc) build_out = argv[++i];
        else if (a == "--window" && i + 1 < argc)
            window = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        else files.push_back(a);
    }
    if (window == 0) window = 1;
    if (persample) return PerSample(window);
    if (!build_out.empty()) return Build(build_out, window, files);
    PrintUsage();
    return 2;
}

}  // namespace llmap::cli
