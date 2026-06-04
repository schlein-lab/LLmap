// LLmap — junction_hunter Mode-5 caller tests.
//
// Synthetic NAHR locus: a 2.4 kb LCR_up, an LCR_down that is LCR_up with
// a PSV every 30 bp (~97 % identity, so every k≥51 window is unique to
// its copy while many k=21 windows are shared), and an independent
// interior. Reads are stitched from these to exercise each call class.
//
// The point of the suite is the discrimination the earlier prototype
// could not make: a chimeric artefact has junction-like k-mer membership
// but scrambled in-region offsets, and must be rejected via the real
// Spearman monotonicity test rather than waved through.

#include "junction_hunter/pair_kmer_index.h"
#include "junction_hunter/read_tiler.h"
#include "junction_hunter/consensus_caller.h"
#include "junction_hunter/junction_hunter_types.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <string>

using namespace llmap::junction_hunter;

namespace {

std::string RandomDna(std::size_t n, std::mt19937_64& rng) {
    static const char b[4] = {'A', 'C', 'G', 'T'};
    std::string s(n, 'A');
    for (std::size_t i = 0; i < n; ++i) s[i] = b[rng() & 3u];
    return s;
}

// Copy `up` and mutate one base every `period` to a different base →
// paralogous-sequence variants spaced so long k are copy-unique.
std::string MakeDownWithPsv(const std::string& up, std::size_t period,
                            std::mt19937_64& rng) {
    static const char b[4] = {'A', 'C', 'G', 'T'};
    std::string dn = up;
    for (std::size_t p = period; p < dn.size(); p += period) {
        char c = dn[p];
        char nc = c;
        while (nc == c) nc = b[rng() & 3u];
        dn[p] = nc;
    }
    return dn;
}

struct Fixture {
    std::string up, down, interior;
    NahrPair pair;
    PairKmerIndex pki;
    MultiKConfig cfg;
};

Fixture MakeFixture(std::size_t psv_period = 30, bool identical = false) {
    std::mt19937_64 rng(0xC0FFEEu + psv_period + (identical ? 7u : 0u));
    Fixture f;
    f.up = RandomDna(2400, rng);
    f.down = identical ? f.up : MakeDownWithPsv(f.up, psv_period, rng);
    f.interior = RandomDna(3000, rng);

    f.pair.pair_id = "NAHR_TEST_01";
    f.pair.chrom = "chrT";
    f.pair.lcr_up_start = 1'000'000;
    f.pair.lcr_up_end = 1'000'000 + f.up.size();
    f.pair.lcr_down_start = 1'200'000;
    f.pair.lcr_down_end = 1'200'000 + f.down.size();
    f.pair.interior_start = 1'100'000;
    f.pair.interior_end = 1'100'000 + f.interior.size();
    f.pair.lcr_identity = identical ? 1.0f : 0.97f;

    f.cfg.k_values = {21, 31, 51, 71, 101};
    f.cfg.consensus_min = 3;
    f.cfg.monotonicity_min = 0.95f;
    f.cfg.min_psv_switches = 3;
    f.cfg.psv_k_min = 51;

    BuildPairKmerIndex(f.pair, f.up, f.down, f.interior, f.cfg, f.pki);
    return f;
}

JunctionRecord Call(const Fixture& f, const std::string& read,
                    const MultiKConfig* override = nullptr) {
    const MultiKConfig& cfg = override ? *override : f.cfg;
    ReadTiling tiling = TileRead(read, cfg);
    return CallJunction("read", tiling, f.pki, f.pair, cfg);
}

}  // namespace

TEST(JunctionHunter, CanonicalUp) {
    auto f = MakeFixture();
    auto rec = Call(f, f.up.substr(600, 400));
    EXPECT_EQ(rec.call, JunctionCall::CanonicalUp);
    EXPECT_GE(rec.n_psv_up, 3u);
    EXPECT_EQ(rec.n_psv_dn, 0u);
    EXPECT_GE(rec.up_monotonicity, 0.95f);
}

TEST(JunctionHunter, CanonicalDown) {
    auto f = MakeFixture();
    auto rec = Call(f, f.down.substr(700, 400));
    EXPECT_EQ(rec.call, JunctionCall::CanonicalDown);
    EXPECT_GE(rec.n_psv_dn, 3u);
    EXPECT_EQ(rec.n_psv_up, 0u);
}

TEST(JunctionHunter, CanonicalInterior) {
    auto f = MakeFixture();
    auto rec = Call(f, f.interior.substr(500, 600));
    EXPECT_EQ(rec.call, JunctionCall::CanonicalInterior);
}

TEST(JunctionHunter, JunctionReal) {
    auto f = MakeFixture();
    // Homologous continuation: up[600..900) then down[900..1300).
    std::string read = f.up.substr(600, 300) + f.down.substr(900, 400);
    auto rec = Call(f, read);
    EXPECT_EQ(rec.call, JunctionCall::JunctionReal);
    EXPECT_GE(rec.up_monotonicity, 0.95f);
    EXPECT_GE(rec.dn_monotonicity, 0.95f);
    EXPECT_GE(rec.n_psv_up, 3u);
    EXPECT_GE(rec.n_psv_dn, 3u);
    // Breakpoint should land near the stitch (read pos ~300) and map into
    // the up/down reference windows.
    EXPECT_GT(rec.breakpoint_read_pos, 150u);
    EXPECT_LT(rec.breakpoint_read_pos, 450u);
    EXPECT_GE(rec.breakpoint_genomic_up, f.pair.lcr_up_start);
    EXPECT_GE(rec.breakpoint_genomic_dn, f.pair.lcr_down_start);
}

TEST(JunctionHunter, ChimeraArtifactRejectedByMonotonicity) {
    auto f = MakeFixture();
    // up part monotonic, down part = 6 blocks from scrambled offsets so
    // the down offsets zig-zag → |Spearman ρ| collapses → ChimeraArtifact.
    std::string read = f.up.substr(600, 300);
    const std::size_t starts[6] = {1500, 300, 1400, 400, 1300, 500};
    for (std::size_t s : starts) read += f.down.substr(s, 130);
    auto rec = Call(f, read);
    EXPECT_EQ(rec.call, JunctionCall::ChimeraArtifact);
    EXPECT_LT(rec.dn_monotonicity, 0.95f);  // the broken half
}

TEST(JunctionHunter, ParalogAmbiguousWhenCopiesIdentical) {
    auto f = MakeFixture(/*psv_period=*/30, /*identical=*/true);
    auto rec = Call(f, f.up.substr(400, 500));
    // Every k-mer is shared between the two identical copies → no
    // unambiguous signal → cannot decide.
    EXPECT_EQ(rec.call, JunctionCall::ParalogAmbiguous);
}

TEST(JunctionHunter, UnmappedWhenUnrelated) {
    auto f = MakeFixture();
    std::mt19937_64 rng(99u);
    auto rec = Call(f, RandomDna(500, rng));
    EXPECT_EQ(rec.call, JunctionCall::Unmapped);
}

TEST(JunctionHunter, MonotonicButThinPsvDemotedToAmbiguous) {
    auto f = MakeFixture();
    // A genuine monotonic junction, but require an impossibly high PSV
    // count → geometry passes, evidence is judged insufficient.
    std::string read = f.up.substr(600, 300) + f.down.substr(900, 400);
    MultiKConfig strict = f.cfg;
    strict.min_psv_switches = 100000;
    auto rec = Call(f, read, &strict);
    EXPECT_EQ(rec.call, JunctionCall::ParalogAmbiguous);
    EXPECT_GE(rec.up_monotonicity, 0.95f);
    EXPECT_GE(rec.dn_monotonicity, 0.95f);
}
