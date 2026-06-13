// LLmap — provenance_baseline tests (population baseline + QC-vs-expected).

#include "provenance/provenance_baseline.h"

#include <gtest/gtest.h>

#include <vector>

namespace {

using llmap::provenance::AggregateBaseline;
using llmap::provenance::QcAgainstBaseline;
using llmap::provenance::QcVerdict;
using llmap::provenance::QcVerdictName;
using llmap::provenance::SampleSpectrum;

std::vector<SampleSpectrum> ThreeSamples() {
    return {
        SampleSpectrum{{"host", 0.90}, {"chim", 0.01}, {"exo:ebv", 0.06}},
        SampleSpectrum{{"host", 0.92}, {"chim", 0.02}},  // no exo:ebv → 0
        SampleSpectrum{{"host", 0.94}, {"chim", 0.03}},
    };
}

TEST(ProvenanceBaseline, AggregateMeanAndSd) {
    const auto b = AggregateBaseline(ThreeSamples());
    EXPECT_EQ(b.n_samples, 3u);
    // host: 0.90/0.92/0.94 → mean 0.92, sample-sd 0.02
    EXPECT_NEAR(b.per_class.at("host").mean, 0.92, 1e-9);
    EXPECT_NEAR(b.per_class.at("host").sd, 0.02, 1e-9);
    // chim: 0.01/0.02/0.03 → mean 0.02, sd 0.01
    EXPECT_NEAR(b.per_class.at("chim").mean, 0.02, 1e-9);
    EXPECT_NEAR(b.per_class.at("chim").sd, 0.01, 1e-9);
    // exo:ebv present in only 1/3 samples (0.06, 0, 0) → mean 0.02
    EXPECT_NEAR(b.per_class.at("exo:ebv").mean, 0.02, 1e-9);
}

TEST(ProvenanceBaseline, NormalSamplePasses) {
    const auto b = AggregateBaseline(ThreeSamples());
    const SampleSpectrum s{{"host", 0.92}, {"chim", 0.02}, {"exo:ebv", 0.02}};
    const auto qc = QcAgainstBaseline(s, b);
    for (const auto& e : qc) {
        EXPECT_EQ(e.verdict, QcVerdict::Pass) << e.cls << " z=" << e.zscore;
    }
}

TEST(ProvenanceBaseline, OutlierFlags) {
    const auto b = AggregateBaseline(ThreeSamples());
    // chim wildly above baseline (0.50 vs mean 0.02, sd 0.01) → huge z → Flag.
    const SampleSpectrum s{{"host", 0.50}, {"chim", 0.50}};
    const auto qc = QcAgainstBaseline(s, b, /*warn_z=*/2.0, /*flag_z=*/4.0);
    bool chim_flagged = false;
    for (const auto& e : qc) {
        if (e.cls == "chim") {
            chim_flagged = e.verdict == QcVerdict::Flag;
            EXPECT_GT(e.zscore, 4.0);
        }
    }
    EXPECT_TRUE(chim_flagged);
}

TEST(ProvenanceBaseline, ClassNeverInBaselineFlags) {
    const auto b = AggregateBaseline(ThreeSamples());
    // A class the pangenome never showed (sd 0, expected 0) but present here →
    // infinite z → Flag (e.g. an exogenous taxon contamination).
    const SampleSpectrum s{{"host", 0.80}, {"exo:mycoplasma", 0.20}};
    const auto qc = QcAgainstBaseline(s, b);
    bool myco_flagged = false;
    for (const auto& e : qc) {
        if (e.cls == "exo:mycoplasma") myco_flagged = e.verdict == QcVerdict::Flag;
    }
    EXPECT_TRUE(myco_flagged);
}

TEST(ProvenanceBaseline, EmptyBaseline) {
    const auto b = AggregateBaseline({});
    EXPECT_EQ(b.n_samples, 0u);
    EXPECT_TRUE(b.per_class.empty());
}

TEST(ProvenanceBaseline, VerdictNames) {
    EXPECT_STREQ(QcVerdictName(QcVerdict::Pass), "PASS");
    EXPECT_STREQ(QcVerdictName(QcVerdict::Warn), "WARN");
    EXPECT_STREQ(QcVerdictName(QcVerdict::Flag), "FLAG");
}

}  // namespace
