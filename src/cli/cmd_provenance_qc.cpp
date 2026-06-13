// LLmap — `llmap provenance-qc`: population-baseline QC of the provenance
// spectrum (the Operator's directive — quantify the buckets against EXPECTED
// values derived by running LLmap over the whole pangenome).
//
//   build mode:  llmap provenance-qc --build-baseline OUT.tsv  S1.tsv S2.tsv …
//                aggregate N per-sample contamination_spectrum.tsv into the
//                pangenome baseline (per-class expected mean + sd).
//   qc mode:     llmap provenance-qc --baseline B.tsv --sample S.tsv
//                z-score + verdict (PASS/WARN/FLAG) of one sample vs the baseline
//                — e.g. an exo class never seen in the pangenome → FLAG.

#include "cli/commands.h"

#include "provenance/provenance_baseline.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace llmap::cli {

namespace {

// Parse a `contamination_spectrum.tsv` (class<TAB>n_reads<TAB>fraction<TAB>…)
// into class→fraction, skipping `#` comments and the header. The Layer-3 overlay
// rows (after the bioconfounder header) are read too — they QC just as well.
provenance::SampleSpectrum ParseSpectrum(const std::string& path) {
    provenance::SampleSpectrum s;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::vector<std::string> t;
        std::string c;
        while (std::getline(ss, c, '\t')) t.push_back(c);
        if (t.size() < 3) continue;
        if (t[0] == "class" || t[0] == "flag") continue;   // header lines
        s[t[0]] = std::strtod(t[2].c_str(), nullptr);       // fraction (col 2)
    }
    return s;
}

void PrintUsage() {
    std::fprintf(stderr,
        "Usage:\n"
        "  llmap provenance-qc --build-baseline OUT.tsv  SPEC1.tsv SPEC2.tsv ...\n"
        "      Aggregate N per-sample spectra (whole pangenome) into the baseline\n"
        "      (per-class expected mean + sd).\n"
        "  llmap provenance-qc --baseline B.tsv --sample S.tsv [--warn-z 2] [--flag-z 4]\n"
        "      QC one sample's spectrum vs the baseline → z-score + verdict.\n");
}

int BuildBaseline(const std::string& out, const std::vector<std::string>& specs) {
    std::vector<provenance::SampleSpectrum> samples;
    samples.reserve(specs.size());
    for (const auto& p : specs) samples.push_back(ParseSpectrum(p));
    const auto base = provenance::AggregateBaseline(samples);

    std::FILE* f = std::fopen(out.c_str(), "w");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", out.c_str()); return 1; }
    std::fprintf(f, "# provenance baseline over %u pangenome samples\n", base.n_samples);
    std::fprintf(f, "class\tmean\tsd\tn_samples\n");
    for (const auto& [cls, st] : base.per_class)
        std::fprintf(f, "%s\t%.8f\t%.8f\t%u\n", cls.c_str(), st.mean, st.sd, st.n_samples);
    std::fclose(f);
    std::fprintf(stderr, "baseline: %u samples, %zu classes → %s\n",
                 base.n_samples, base.per_class.size(), out.c_str());
    return 0;
}

provenance::ProvenanceBaseline LoadBaseline(const std::string& path) {
    provenance::ProvenanceBaseline b;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::vector<std::string> t;
        std::string c;
        while (std::getline(ss, c, '\t')) t.push_back(c);
        if (t.size() < 4 || t[0] == "class") continue;
        provenance::ClassStat st;
        st.mean = std::strtod(t[1].c_str(), nullptr);
        st.sd = std::strtod(t[2].c_str(), nullptr);
        st.n_samples = static_cast<std::uint32_t>(std::strtoul(t[3].c_str(), nullptr, 10));
        b.per_class[t[0]] = st;
        b.n_samples = std::max(b.n_samples, st.n_samples);
    }
    return b;
}

int QcSample(const std::string& baseline, const std::string& sample,
             double warn_z, double flag_z) {
    const auto base = LoadBaseline(baseline);
    const auto spec = ParseSpectrum(sample);
    const auto qc = provenance::QcAgainstBaseline(spec, base, warn_z, flag_z);

    std::printf("class\tobserved\texpected\tsd\tz\tverdict\n");
    int n_flag = 0;
    for (const auto& e : qc) {
        std::printf("%s\t%.6f\t%.6f\t%.6f\t%+.2f\t%s\n", e.cls.c_str(), e.observed,
                    e.expected, e.sd, e.zscore, provenance::QcVerdictName(e.verdict));
        if (e.verdict == provenance::QcVerdict::Flag) ++n_flag;
    }
    std::fprintf(stderr, "provenance-qc: %zu classes, %d FLAG vs baseline (%u samples)\n",
                 qc.size(), n_flag, base.n_samples);
    return n_flag > 0 ? 10 : 0;   // non-zero exit when a class is flagged
}

}  // namespace

int run_provenance_qc(int argc, char** argv) {
    std::string build_out, baseline, sample;
    double warn_z = 2.0, flag_z = 4.0;
    std::vector<std::string> specs;
    for (int i = 0; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") { PrintUsage(); return 0; }
        else if (a == "--build-baseline" && i + 1 < argc) build_out = argv[++i];
        else if (a == "--baseline" && i + 1 < argc) baseline = argv[++i];
        else if (a == "--sample" && i + 1 < argc) sample = argv[++i];
        else if (a == "--warn-z" && i + 1 < argc) warn_z = std::strtod(argv[++i], nullptr);
        else if (a == "--flag-z" && i + 1 < argc) flag_z = std::strtod(argv[++i], nullptr);
        else specs.push_back(a);   // positional = spectrum files for --build-baseline
    }
    if (!build_out.empty()) return BuildBaseline(build_out, specs);
    if (!baseline.empty() && !sample.empty()) return QcSample(baseline, sample, warn_z, flag_z);
    PrintUsage();
    return 2;
}

}  // namespace llmap::cli
