// LLmap — Provenance population baseline + per-sample QC-vs-expected.
//
// Operator directive (2026-06-13): quantify the per-read provenance buckets and
// output them as a QC statistic AGAINST EXPECTED values — where the expectation
// is derived by applying LLmap across the WHOLE pangenome. So:
//   1. run `provenance-spectrum` on every pangenome sample → per-sample spectrum
//   2. AggregateBaseline() over all of them → per-class expected mean + sd
//   3. QcAgainstBaseline() a sample vs the baseline → per-class z-score + verdict
//      (e.g. exo:ebv at 5σ above the pangenome baseline = a contamination alert;
//       chim/dmg within baseline = OK).
//
// This is the population-QC layer the contamination/provenance mode was built
// for: a sample's spectrum only becomes actionable against an expected
// distribution. Pure stats over the spectra (stdlib only) — the CLI feeds it the
// per-sample fraction maps; running LLmap across the pangenome is the compute
// side.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace llmap::provenance {

// One sample's observed spectrum: provenance class (PV string) → read fraction.
using SampleSpectrum = std::map<std::string, double>;

// Expected distribution of one class across the pangenome.
struct ClassStat {
    double        mean{0.0};       // mean fraction across samples
    double        sd{0.0};         // sample standard deviation
    std::uint32_t n_samples{0};    // samples that contributed (== baseline size)
};

struct ProvenanceBaseline {
    std::map<std::string, ClassStat> per_class;
    std::uint32_t n_samples{0};
};

// Aggregate per-sample spectra into the population baseline. A class absent from
// a sample contributes a 0 fraction (it WAS observed at 0), so the mean/sd are
// over all samples, not just those where the class appeared.
[[nodiscard]] ProvenanceBaseline AggregateBaseline(
    const std::vector<SampleSpectrum>& samples);

enum class QcVerdict : std::uint8_t { Pass, Warn, Flag };
[[nodiscard]] const char* QcVerdictName(QcVerdict v) noexcept;

struct QcEntry {
    std::string cls;
    double      observed{0.0};   // sample fraction
    double      expected{0.0};   // baseline mean
    double      sd{0.0};         // baseline sd
    double      zscore{0.0};     // (observed - expected) / sd
    QcVerdict   verdict{QcVerdict::Pass};
};

// QC a sample against the baseline. Reports every class in the union of the
// sample and the baseline (so a class never seen in the pangenome but present
// in the sample — e.g. an exogenous taxon — flags strongly). When the baseline
// sd is 0 (a class that never varies), any positive observed excess flags.
// Verdicts: |z| >= flag_z ⇒ Flag, >= warn_z ⇒ Warn, else Pass.
[[nodiscard]] std::vector<QcEntry> QcAgainstBaseline(
    const SampleSpectrum& sample, const ProvenanceBaseline& baseline,
    double warn_z = 2.0, double flag_z = 4.0);

}  // namespace llmap::provenance
