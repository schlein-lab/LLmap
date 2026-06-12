// LLmap — transcript_stage unit tests.
//
// Pure-logic tests with synthetic reference strings + a fake junction scorer.
// No genomic data, no I/O.

#include "mapping/transcript_stage.h"
#include "mapping/chain_spliced.h"

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <vector>

namespace {

using llmap::mapping::ApplyTranscriptStage;
using llmap::mapping::ExtractJunctionMotifs;
using llmap::mapping::JunctionScorer;
using llmap::mapping::LinearSubChain;
using llmap::mapping::RefSeqLookup;
using llmap::mapping::SplicedAlignment;
using llmap::mapping::TranscriptStageConfig;

// A reference of `len` 'A's with optional 2-bp patches at given offsets.
std::string MakeRef(std::size_t len) { return std::string(len, 'A'); }

void Patch(std::string& ref, std::size_t pos, std::string_view motif) {
    for (std::size_t i = 0; i < motif.size(); ++i) ref[pos + i] = motif[i];
}

LinearSubChain Sub(std::string ref_id, std::uint64_t rs, std::uint64_t re,
                   std::uint32_t qs, std::uint32_t qe, char strand,
                   std::string cigar) {
    LinearSubChain s;
    s.ref_id = std::move(ref_id);
    s.ref_start = rs;
    s.ref_end = re;
    s.query_start = qs;
    s.query_end = qe;
    s.strand = strand;
    s.score = 100;
    s.cigar = std::move(cigar);
    return s;
}

// Fake scorer: high when canonical GT..AG sense motifs, low otherwise.
JunctionScorer CanonicalScorer() {
    return [](std::string_view donor2, std::string_view acceptor2,
              std::string_view, std::string_view) -> float {
        return (donor2 == "GT" && acceptor2 == "AG") ? 0.95f : 0.05f;
    };
}

RefSeqLookup LookupFor(const std::string& name, const std::string& seq) {
    return [name, seq](std::string_view id) -> std::string_view {
        return id == name ? std::string_view(seq) : std::string_view{};
    };
}

// ---------------------------------------------------------------------------
// ExtractJunctionMotifs
// ---------------------------------------------------------------------------

TEST(ExtractMotifs, ForwardCanonical) {
    // exon[0,100) intron[100,1100) exon[1100,1200); GT at 100, AG at 1098.
    std::string ref = MakeRef(1200);
    Patch(ref, 100, "GT");
    Patch(ref, 1098, "AG");
    auto a = Sub("chr1", 0, 100, 0, 100, '+', "100M");
    auto b = Sub("chr1", 1100, 1200, 100, 200, '+', "100M");
    std::string d, ac, i3, i5;
    ASSERT_TRUE(ExtractJunctionMotifs(a, b, ref, d, ac, i3, i5));
    EXPECT_EQ(d, "GT");
    EXPECT_EQ(ac, "AG");
}

TEST(ExtractMotifs, ReverseComplemented) {
    // For '-' strand, sense donor is at b side: ref[1098,1100) revcomp == GT
    //   => ref[1098,1100) = "AC". Sense acceptor at a side: ref[100,102)
    //   revcomp == AG => ref[100,102) = "CT".
    std::string ref = MakeRef(1200);
    Patch(ref, 1098, "AC");
    Patch(ref, 100, "CT");
    auto a = Sub("chr1", 0, 100, 0, 100, '-', "100M");
    auto b = Sub("chr1", 1100, 1200, 100, 200, '-', "100M");
    std::string d, ac, i3, i5;
    ASSERT_TRUE(ExtractJunctionMotifs(a, b, ref, d, ac, i3, i5));
    EXPECT_EQ(d, "GT");
    EXPECT_EQ(ac, "AG");
}

TEST(ExtractMotifs, OutOfBoundsRejected) {
    std::string ref = MakeRef(150);  // too short for b at 1100
    auto a = Sub("chr1", 0, 100, 0, 100, '+', "100M");
    auto b = Sub("chr1", 1100, 1200, 100, 200, '+', "100M");
    std::string d, ac, i3, i5;
    EXPECT_FALSE(ExtractJunctionMotifs(a, b, ref, d, ac, i3, i5));
}

// ---------------------------------------------------------------------------
// ApplyTranscriptStage
// ---------------------------------------------------------------------------

TEST(TranscriptStage, MergesCanonicalJunction) {
    std::string ref = MakeRef(1200);
    Patch(ref, 100, "GT");
    Patch(ref, 1098, "AG");
    std::vector<LinearSubChain> subs = {
        Sub("chr1", 0, 100, 0, 100, '+', "100M"),
        Sub("chr1", 1100, 1200, 100, 200, '+', "100M"),
    };
    TranscriptStageConfig cfg;  // default min_junction_probability 0.30
    auto out = ApplyTranscriptStage(subs, LookupFor("chr1", ref),
                                    CanonicalScorer(), cfg);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(out[0].is_spliced);
    EXPECT_EQ(out[0].ref_start, 0u);
    EXPECT_EQ(out[0].ref_end, 1200u);
    EXPECT_EQ(out[0].cigar, "100M1000N100M");
    ASSERT_EQ(out[0].junctions.size(), 1u);
    EXPECT_EQ(out[0].junctions[0].first, 99u);    // donor_ref_pos = a.ref_end-1
    EXPECT_EQ(out[0].junctions[0].second, 1100u);  // acceptor_ref_pos = b.ref_start
}

TEST(TranscriptStage, NonCanonicalStaysSeparate) {
    // No GT/AG patch → scorer returns 0.05 < 0.30 → no merge → 2 singletons.
    std::string ref = MakeRef(1200);
    std::vector<LinearSubChain> subs = {
        Sub("chr1", 0, 100, 0, 100, '+', "100M"),
        Sub("chr1", 1100, 1200, 100, 200, '+', "100M"),
    };
    auto out = ApplyTranscriptStage(subs, LookupFor("chr1", ref),
                                    CanonicalScorer(), {});
    ASSERT_EQ(out.size(), 2u);
    EXPECT_FALSE(out[0].is_spliced);
    EXPECT_FALSE(out[1].is_spliced);
}

TEST(TranscriptStage, SingleSubChainUnchanged) {
    std::string ref = MakeRef(200);
    std::vector<LinearSubChain> subs = {
        Sub("chr1", 0, 100, 0, 100, '+', "100M"),
    };
    auto out = ApplyTranscriptStage(subs, LookupFor("chr1", ref),
                                    CanonicalScorer(), {});
    ASSERT_EQ(out.size(), 1u);
    EXPECT_FALSE(out[0].is_spliced);
    EXPECT_EQ(out[0].cigar, "100M");
}

TEST(TranscriptStage, ReverseStrandMerges) {
    std::string ref = MakeRef(1200);
    Patch(ref, 1098, "AC");  // sense donor GT after revcomp
    Patch(ref, 100, "CT");   // sense acceptor AG after revcomp
    std::vector<LinearSubChain> subs = {
        Sub("chr1", 0, 100, 0, 100, '-', "100M"),
        Sub("chr1", 1100, 1200, 100, 200, '-', "100M"),
    };
    auto out = ApplyTranscriptStage(subs, LookupFor("chr1", ref),
                                    CanonicalScorer(), {});
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(out[0].is_spliced);
    EXPECT_EQ(out[0].strand, '-');
    EXPECT_EQ(out[0].cigar, "100M1000N100M");
}

TEST(TranscriptStage, WeakMotifStillMergesWithNop) {
    // R-A: intron-like geometry but weak/non-canonical motif (scorer 0.05).
    // With the transcript-mode floor threshold the read still merges into ONE
    // spliced alignment, the gap is an N op (intron skip, not a D deletion),
    // and the low motif confidence is preserved for the jM tag.
    std::string ref = MakeRef(1200);  // no GT/AG patch → scorer returns 0.05
    std::vector<LinearSubChain> subs = {
        Sub("chr1", 0, 100, 0, 100, '+', "100M"),
        Sub("chr1", 1100, 1200, 100, 200, '+', "100M"),
    };
    TranscriptStageConfig cfg;
    cfg.joiner.min_junction_probability = 0.05f;  // motif = evidence, not gate
    auto out = ApplyTranscriptStage(subs, LookupFor("chr1", ref),
                                    CanonicalScorer(), cfg);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_TRUE(out[0].is_spliced);
    EXPECT_EQ(out[0].cigar, "100M1000N100M");  // N op despite weak motif
    ASSERT_EQ(out[0].junction_conf.size(), 1u);
    EXPECT_LT(out[0].junction_conf[0], 0.30f);  // low confidence preserved
}

TEST(TranscriptStage, EmptyInput) {
    auto out = ApplyTranscriptStage({}, LookupFor("chr1", ""), CanonicalScorer(),
                                    {});
    EXPECT_TRUE(out.empty());
}

}  // namespace
