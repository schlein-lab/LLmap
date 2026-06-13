// Unit tests for splice-site (GT/AG) boundary snapping.

#include "mapping/splice_snap.h"

#include <gtest/gtest.h>

#include <string>

namespace llmap::mapping {
namespace {

LinearSubChain Sub(std::uint64_t rs, std::uint64_t re,
                   std::uint32_t qs, std::uint32_t qe, char strand = '+') {
    LinearSubChain c;
    c.ref_start = rs; c.ref_end = re;
    c.query_start = qs; c.query_end = qe;
    c.strand = strand;
    return c;
}

TEST(SpliceSnap, ForwardCanonicalAtAbuttingBoundary) {
    // exon-a ref[0,100] q[0,100], exon-b ref[200,300] q[100,200]; intron ref[100,200].
    std::string ref(300, 'A');
    ref[100] = 'G'; ref[101] = 'T';   // donor GT at intron 5'
    ref[198] = 'A'; ref[199] = 'G';   // acceptor AG at intron 3'
    auto r = SnapJunction(Sub(0, 100, 0, 100), Sub(200, 300, 100, 200), ref, '+');
    EXPECT_TRUE(r.snapped);
    EXPECT_EQ(r.donor_ref_pos, 100u);
    EXPECT_EQ(r.acceptor_ref_pos, 200u);
    EXPECT_EQ(r.query_split, 100u);
    EXPECT_FLOAT_EQ(r.motif_score, 1.0f);
}

TEST(SpliceSnap, ClosesSlopGapByExtendingBothExons) {
    // Seed under-aligns: exon-a ends at ref/q 98, exon-b starts at ref/q 202
    // (q_gap = 4). True donor GT@100, acceptor AG@198 → snap extends each exon
    // by 2 and the read split lands at q=100 (gap closed).
    std::string ref(300, 'A');
    ref[100] = 'G'; ref[101] = 'T';
    ref[198] = 'A'; ref[199] = 'G';
    auto r = SnapJunction(Sub(0, 98, 0, 98), Sub(202, 300, 102, 200), ref, '+');
    EXPECT_TRUE(r.snapped);
    EXPECT_EQ(r.donor_ref_pos, 100u);     // exon-a extended +2
    EXPECT_EQ(r.acceptor_ref_pos, 200u);  // exon-b extended back +2
    EXPECT_EQ(r.query_split, 100u);       // single split → q_gap = 0
}

TEST(SpliceSnap, ReverseStrandCtAc) {
    // '-' transcript: forward-ref intron motif is CT..AC.
    std::string ref(300, 'A');
    ref[100] = 'C'; ref[101] = 'T';
    ref[198] = 'A'; ref[199] = 'C';
    auto r = SnapJunction(Sub(0, 100, 0, 100, '-'), Sub(200, 300, 100, 200, '-'),
                          ref, '-');
    EXPECT_TRUE(r.snapped);
    EXPECT_EQ(r.donor_ref_pos, 100u);
    EXPECT_EQ(r.acceptor_ref_pos, 200u);
}

TEST(SpliceSnap, NoCanonicalSiteIsNoOp) {
    std::string ref(300, 'A');         // no GT/AG motif anywhere in the window
    auto r = SnapJunction(Sub(0, 100, 0, 100), Sub(200, 300, 100, 200), ref, '+');
    EXPECT_FALSE(r.snapped);
    EXPECT_EQ(r.donor_ref_pos, 100u);  // boundary unchanged
    EXPECT_EQ(r.acceptor_ref_pos, 200u);
    EXPECT_EQ(r.query_split, 100u);
}

TEST(SpliceSnap, PrefersSmallestDisplacement) {
    // Canonical at the abutting boundary (q=100, displacement 0) AND further out
    // (q=110). The nearer one must win.
    std::string ref(400, 'A');
    ref[100] = 'G'; ref[101] = 'T'; ref[198] = 'A'; ref[199] = 'G';  // disp 0
    ref[110] = 'G'; ref[111] = 'T'; ref[208] = 'A'; ref[209] = 'G';  // disp 20
    auto r = SnapJunction(Sub(0, 100, 0, 100), Sub(200, 320, 100, 220), ref, '+');
    EXPECT_TRUE(r.snapped);
    EXPECT_EQ(r.donor_ref_pos, 100u);
    EXPECT_EQ(r.query_split, 100u);
}

}  // namespace
}  // namespace llmap::mapping
