// LLmap — Spliced-chain joiner tests.
//
// Synthesise short sequences of LinearSubChain values that exercise
// the joiner's merge / flush decisions, then check the resulting
// SplicedChainResult + EmitSplicedCigar.

#include "mapping/chain_spliced.h"

#include <gtest/gtest.h>

using namespace llmap::mapping;

namespace {

LinearSubChain MakeSub(std::string ref,
                        std::uint64_t r_start, std::uint64_t r_end,
                        std::uint32_t q_start, std::uint32_t q_end,
                        std::string cigar = "100M",
                        char strand = '+') {
    LinearSubChain s;
    s.ref_id = std::move(ref);
    s.ref_start = r_start;
    s.ref_end = r_end;
    s.query_start = q_start;
    s.query_end = q_end;
    s.cigar = std::move(cigar);
    s.score = static_cast<std::int32_t>(r_end - r_start);
    s.strand = strand;
    return s;
}

}  // namespace

TEST(ChainSpliced, EmptyInputYieldsEmptyOutput) {
    auto out = JoinSplicedChains({}, {}, {});
    EXPECT_TRUE(out.chains.empty());
    EXPECT_EQ(out.n_singletons_kept, 0u);
}

TEST(ChainSpliced, SingletonYieldsSingleSplicedChain) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
    };
    auto out = JoinSplicedChains(in, {}, {});
    ASSERT_EQ(out.chains.size(), 1u);
    EXPECT_EQ(out.chains[0].sub_chains.size(), 1u);
    EXPECT_EQ(out.chains[0].junctions.size(), 0u);
    EXPECT_EQ(out.n_singletons_kept, 1u);
}

TEST(ChainSpliced, MergesIntronGapWhenProbabilityHigh) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
        MakeSub("chr14", 8100, 8200, 100, 200),  // 7000-bp intron
    };
    std::vector<float> p = {0.95f};
    JoinerConfig cfg;
    auto out = JoinSplicedChains(in, p, cfg);
    ASSERT_EQ(out.chains.size(), 1u);
    EXPECT_EQ(out.chains[0].sub_chains.size(), 2u);
    ASSERT_EQ(out.chains[0].junctions.size(), 1u);
    EXPECT_TRUE(out.chains[0].junctions[0].is_confirmed);
    EXPECT_NEAR(out.chains[0].junctions[0].probability, 0.95f, 1e-5f);
}

TEST(ChainSpliced, RefusesMergeWhenProbabilityBelowThreshold) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
        MakeSub("chr14", 8100, 8200, 100, 200),
    };
    std::vector<float> p = {0.10f};
    JoinerConfig cfg;
    cfg.min_junction_probability = 0.30f;
    auto out = JoinSplicedChains(in, p, cfg);
    EXPECT_EQ(out.chains.size(), 2u);
    EXPECT_EQ(out.n_singletons_kept, 2u);
}

TEST(ChainSpliced, RefusesMergeAcrossChroms) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
        MakeSub("chr15", 8100, 8200, 100, 200),
    };
    std::vector<float> probs95 = {0.95f};
    auto out = JoinSplicedChains(in, probs95, {});
    EXPECT_EQ(out.chains.size(), 2u);
    EXPECT_EQ(out.n_singletons_kept, 2u);
}

TEST(ChainSpliced, RefusesMergeOnLongQueryGap) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
        MakeSub("chr14", 8100, 8200, 500, 600),   // query gap 400 → too big
    };
    std::vector<float> probs95 = {0.95f};
    auto out = JoinSplicedChains(in, probs95, {});
    EXPECT_EQ(out.chains.size(), 2u);
}

TEST(ChainSpliced, RefusesMergeOnGapBelowMinIntron) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100, 0, 100),
        MakeSub("chr14", 1110, 1210, 100, 200),  // gap 10 < min 50
    };
    std::vector<float> probs95 = {0.95f};
    auto out = JoinSplicedChains(in, probs95, {});
    EXPECT_EQ(out.chains.size(), 2u);
}

TEST(ChainSpliced, RefusesMergeOnGapAboveMaxIntron) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14",      1000, 1100, 0, 100),
        MakeSub("chr14", 50'000'000, 50'000'100, 100, 200),  // 50 Mb intron
    };
    std::vector<float> probs95 = {0.95f};
    auto out = JoinSplicedChains(in, probs95, {});
    EXPECT_EQ(out.chains.size(), 2u);
}

TEST(ChainSpliced, MultipleJunctionsInOneChain) {
    std::vector<LinearSubChain> in = {
        MakeSub("chr14", 1000, 1100,  0, 100, "100M"),
        MakeSub("chr14", 5000, 5100, 100, 200, "100M"),
        MakeSub("chr14", 9000, 9100, 200, 300, "100M"),
    };
    auto out = JoinSplicedChains(in, std::vector<float>{0.9f, 0.9f}, {});
    ASSERT_EQ(out.chains.size(), 1u);
    EXPECT_EQ(out.chains[0].sub_chains.size(), 3u);
    EXPECT_EQ(out.chains[0].junctions.size(), 2u);
}

TEST(ChainSpliced, EmitCigarContainsNOpForConfirmedJunction) {
    SplicedChain sc;
    sc.sub_chains.push_back(MakeSub("chr14", 1000, 1100, 0, 100, "100M"));
    sc.sub_chains.push_back(MakeSub("chr14", 8100, 8200, 100, 200, "100M"));
    Junction j;
    j.donor_ref_pos = 1099;
    j.acceptor_ref_pos = 8100;
    j.probability = 0.95f;
    j.is_confirmed = true;
    sc.junctions.push_back(j);

    auto cigar = EmitSplicedCigar(sc);
    EXPECT_NE(cigar.find('N'), std::string::npos)
        << "confirmed junction → CIGAR-N; got: " << cigar;
    // 8100 - 1100 = 7000 bp intron
    EXPECT_NE(cigar.find("7000N"), std::string::npos)
        << "expected '7000N' in CIGAR; got: " << cigar;
}

TEST(ChainSpliced, EmitCigarUsesDForUnconfirmedJunction) {
    SplicedChain sc;
    sc.sub_chains.push_back(MakeSub("chr14", 1000, 1100, 0, 100, "100M"));
    sc.sub_chains.push_back(MakeSub("chr14", 8100, 8200, 100, 200, "100M"));
    Junction j;
    j.donor_ref_pos = 1099;
    j.acceptor_ref_pos = 8100;
    j.probability = 0.20f;
    j.is_confirmed = false;
    sc.junctions.push_back(j);

    auto cigar = EmitSplicedCigar(sc);
    EXPECT_NE(cigar.find("7000D"), std::string::npos)
        << "unconfirmed junction → CIGAR-D; got: " << cigar;
    EXPECT_EQ(cigar.find('N'), std::string::npos);
}

TEST(ChainSpliced, EmitCigarSingleSubChainReturnsSubCigar) {
    SplicedChain sc;
    sc.sub_chains.push_back(MakeSub("chr14", 1000, 1100, 0, 100, "100M"));
    EXPECT_EQ(EmitSplicedCigar(sc), "100M");
}
