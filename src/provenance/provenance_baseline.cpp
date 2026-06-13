// LLmap — Provenance population baseline + QC implementation.

#include "provenance/provenance_baseline.h"

#include <cmath>
#include <limits>
#include <set>

namespace llmap::provenance {

const char* QcVerdictName(QcVerdict v) noexcept {
    switch (v) {
        case QcVerdict::Pass: return "PASS";
        case QcVerdict::Warn: return "WARN";
        case QcVerdict::Flag: return "FLAG";
    }
    return "PASS";
}

ProvenanceBaseline AggregateBaseline(const std::vector<SampleSpectrum>& samples) {
    ProvenanceBaseline b;
    b.n_samples = static_cast<std::uint32_t>(samples.size());
    if (samples.empty()) return b;

    std::set<std::string> classes;
    for (const auto& s : samples) {
        for (const auto& [cls, _] : s) classes.insert(cls);
    }

    const double n = static_cast<double>(samples.size());
    for (const auto& cls : classes) {
        double sum = 0.0;
        for (const auto& s : samples) {
            const auto it = s.find(cls);
            sum += (it != s.end()) ? it->second : 0.0;  // absent ⇒ observed at 0
        }
        const double mean = sum / n;

        double ss = 0.0;
        for (const auto& s : samples) {
            const auto it = s.find(cls);
            const double v = (it != s.end()) ? it->second : 0.0;
            ss += (v - mean) * (v - mean);
        }
        ClassStat st;
        st.mean = mean;
        // sample standard deviation (n-1); 0 for a single sample.
        st.sd = (samples.size() > 1) ? std::sqrt(ss / (n - 1.0)) : 0.0;
        st.n_samples = b.n_samples;
        b.per_class[cls] = st;
    }
    return b;
}

std::vector<QcEntry> QcAgainstBaseline(const SampleSpectrum& sample,
                                       const ProvenanceBaseline& baseline,
                                       double warn_z, double flag_z) {
    std::set<std::string> classes;
    for (const auto& [cls, _] : sample) classes.insert(cls);
    for (const auto& [cls, _] : baseline.per_class) classes.insert(cls);

    std::vector<QcEntry> out;
    out.reserve(classes.size());
    for (const auto& cls : classes) {
        QcEntry e;
        e.cls = cls;
        const auto sit = sample.find(cls);
        e.observed = (sit != sample.end()) ? sit->second : 0.0;
        const auto bit = baseline.per_class.find(cls);
        e.expected = (bit != baseline.per_class.end()) ? bit->second.mean : 0.0;
        e.sd = (bit != baseline.per_class.end()) ? bit->second.sd : 0.0;

        if (e.sd > 0.0) {
            e.zscore = (e.observed - e.expected) / e.sd;
        } else {
            // No spread in the baseline: any positive excess is an anomaly
            // (e.g. an exogenous class never seen across the pangenome).
            e.zscore = (e.observed > e.expected)
                           ? std::numeric_limits<double>::infinity()
                           : 0.0;
        }

        const double az = std::abs(e.zscore);
        e.verdict = (az >= flag_z) ? QcVerdict::Flag
                    : (az >= warn_z) ? QcVerdict::Warn
                                     : QcVerdict::Pass;
        out.push_back(std::move(e));
    }
    return out;
}

}  // namespace llmap::provenance
