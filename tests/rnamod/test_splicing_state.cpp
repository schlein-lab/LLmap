// LLmap — SplicingStateClassifier tests.

#include "rnamod/splicing_state.h"

#include <gtest/gtest.h>

#include <vector>

using namespace llmap::splicing;

namespace {

ObservedJunction MakeJunction(std::uint64_t donor,
                                std::uint64_t acceptor,
                                std::uint8_t klass = 0,
                                bool annotated = true,
                                bool retained = false) {
    ObservedJunction j;
    j.donor_genomic_pos = donor;
    j.acceptor_genomic_pos = acceptor;
    j.spliceosome_class = klass;
    j.is_annotated = annotated;
    j.is_retained = retained;
    return j;
}

}  // namespace

TEST(SplicingState, EmptyJunctionsYieldsUnspliced) {
    SplicingStateClassifier c;
    auto inf = c.Classify({});
    EXPECT_EQ(inf.dominant, SplicingState::Unspliced);
    EXPECT_NEAR(inf.confidence, 1.0f, 1e-5f);
}

TEST(SplicingState, AllAnnotatedCanonicalScoresHigh) {
    std::vector<ObservedJunction> j = {
        MakeJunction(100, 500),
        MakeJunction(700, 1000),
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::Canonical);
    EXPECT_GT(inf.confidence, 0.9f);
}

TEST(SplicingState, SingleRetainedIntron) {
    std::vector<ObservedJunction> j = {
        MakeJunction(100, 500, /*klass*/ 0, true, /*retained*/ true),
        MakeJunction(700, 1000),
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::IntronRetained);
    ASSERT_EQ(inf.retained_intron_indices.size(), 1u);
    EXPECT_EQ(inf.retained_intron_indices[0], 0u);
}

TEST(SplicingState, MultipleRetainedYieldsPartiallySpliced) {
    std::vector<ObservedJunction> j = {
        MakeJunction(100, 500, 0, true, /*retained*/ true),
        MakeJunction(700, 1000, 0, true, /*retained*/ true),
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::PartiallySpliced);
    EXPECT_EQ(inf.retained_intron_indices.size(), 2u);
}

TEST(SplicingState, BackSpliceYieldsCircular) {
    std::vector<ObservedJunction> j = {
        MakeJunction(5000, 1000, /*klass*/ 3),  // donor>acceptor + back-class
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::BackSplicedCircular);
}

TEST(SplicingState, LongRefGapUnannotatedYieldsRecursive) {
    std::vector<ObservedJunction> j = {
        MakeJunction(100, 200'000, /*klass*/ 0, /*annotated*/ false),
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::RecursiveSpliced);
    ASSERT_TRUE(inf.recursive_splice_site_pos.has_value());
}

TEST(SplicingState, AllNonCanonicalYieldsNovelSplicing) {
    std::vector<ObservedJunction> j = {
        MakeJunction(100, 500, /*klass*/ 2, /*annotated*/ false),
        MakeJunction(700, 1000, /*klass*/ 2, /*annotated*/ false),
    };
    SplicingStateClassifier c;
    auto inf = c.Classify(j);
    EXPECT_EQ(inf.dominant, SplicingState::NovelSplicingState);
}

TEST(SplicingState, NameTableCoversAllValues) {
    using S = SplicingState;
    constexpr S all[] = {
        S::Unknown, S::Canonical, S::Unspliced, S::PartiallySpliced,
        S::IntronRetained, S::Lariat, S::RecursiveSpliced, S::TransSpliced,
        S::AlternativeCassetteIn, S::AlternativeCassetteOut,
        S::Alt5ss, S::Alt3ss, S::MutuallyExclusive,
        S::BackSplicedCircular, S::HalfSplicedCotrans, S::NovelSplicingState,
    };
    for (auto s : all) {
        const char* n = SplicingStateName(s);
        ASSERT_NE(n, nullptr);
        EXPECT_GT(std::string(n).size(), 0u);
    }
}
