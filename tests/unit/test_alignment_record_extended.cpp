// LLmap — extended AlignmentRecord tests for the seven new statuses
// landed in Plan-Block 1 (Transcript-Mode lossless schema).
//
// Splits from test_alignment_record.cpp so that file keeps documenting
// the original 3-status contract verbatim. This file documents how the
// new statuses (MappedSterile, MappedPreMrna, MappedCircular,
// PartialMappedExact, ChimericIntraRegion, ChimericInterChrom,
// DarkNovel) preserve the same invariant guarantees.

#include "core/alignment_record.h"

#include <gtest/gtest.h>

using namespace llmap;

namespace {

AlignmentHit hit(std::string target = "chr14", std::uint64_t s = 105'000'000,
                  std::uint64_t e = 105'010'000) {
    AlignmentHit h;
    h.target_id = std::move(target);
    h.start = s;
    h.end = e;
    h.cigar.ops = "10000M";
    h.score = 9500;
    h.nm = 50;
    return h;
}

}  // namespace

// ===========================================================================
// MappedSterile — IGH sterile germline transcript (Iγ→Sγ→Cγ)
// ===========================================================================

TEST(AlignmentRecordExt, MappedSterileRequiresAllThreeIds) {
    SterileDetail d;
    d.i_promoter_id = "I_gamma4_promoter";
    d.s_region_id   = "S_gamma4";
    d.c_gene_id     = "IGHG4";

    auto r = make_mapped_sterile("read_sg_1", 2'500, hit(), d);
    EXPECT_EQ(r.status, AlignmentStatus::MappedSterile);
    EXPECT_EQ(r.transcript_kind, core::TranscriptKind::SterileGermline);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, MappedSterileMissingIPromoterViolatesInvariant) {
    SterileDetail d;
    // d.i_promoter_id intentionally empty
    d.s_region_id = "S_gamma4";
    d.c_gene_id   = "IGHG4";

    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::MappedSterile;
    r.primary = hit();
    r.sterile = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// MappedPreMrna — retained-intron / pre-mRNA
// ===========================================================================

TEST(AlignmentRecordExt, MappedPreMrnaWithRetainedIntrons) {
    PreMrnaDetail d;
    d.retained_introns = {{500, 1500}, {3000, 3800}};
    d.pct_introns_retained = 0.40f;

    auto r = make_mapped_pre_mrna("read_premrna_1", 5'000, hit(), d);
    EXPECT_EQ(r.status, AlignmentStatus::MappedPreMrna);
    EXPECT_TRUE(r.is_lossless_consistent());
    EXPECT_EQ(r.pre_mrna->retained_introns.size(), 2u);
}

TEST(AlignmentRecordExt, MappedPreMrnaPctOutOfRangeFails) {
    PreMrnaDetail d;
    d.pct_introns_retained = 1.5f;  // > 1.0 is bogus
    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::MappedPreMrna;
    r.primary = hit();
    r.pre_mrna = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// MappedCircular — circRNA back-splice
// ===========================================================================

TEST(AlignmentRecordExt, MappedCircularValid) {
    CircularDetail d;
    d.backsplice_acceptor_pos = 100'000;
    d.backsplice_donor_pos    = 105'000;  // donor > acceptor in linear genome
    d.host_gene_id = "IGHG4";

    auto r = make_mapped_circular("read_circ_1", 1'200, hit(), d);
    EXPECT_EQ(r.status, AlignmentStatus::MappedCircular);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, MappedCircularDonorBeforeAcceptorFails) {
    // donor < acceptor means it's a regular forward splice, not a circle
    CircularDetail d;
    d.backsplice_acceptor_pos = 105'000;
    d.backsplice_donor_pos    = 100'000;
    d.host_gene_id = "IGHG4";
    AlignmentRecord r;
    r.read_id = "read_x";
    r.status = AlignmentStatus::MappedCircular;
    r.primary = hit();
    r.circular = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// PartialMappedExact — 100%-match substring with preserved tail
// ===========================================================================

TEST(AlignmentRecordExt, PartialMappedExactValid) {
    PartialMatchDetail d;
    d.read_offset = 200;
    d.match_length = 250;
    d.anchor_id = "GENCODE:ENST00000390557.4:exon1";
    d.unmatched_tail = "AAACGTGTGATCGAT";

    auto r = make_partial_exact("read_partial_1", 800, d);
    EXPECT_EQ(r.status, AlignmentStatus::PartialMappedExact);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, PartialMappedExactZeroLengthFails) {
    PartialMatchDetail d;
    d.match_length = 0;
    d.anchor_id = "x";
    AlignmentRecord r;
    r.read_id = "x";
    r.status = AlignmentStatus::PartialMappedExact;
    r.partial_match = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// Chimeric — intra-region + cross-chrom + VDJ
// ===========================================================================

TEST(AlignmentRecordExt, ChimericIntraRegionVdjClassSwitch) {
    ChimericDetail d;
    d.parts = {hit("chr14", 105'860'000, 105'860'400),    // IGHM-CH1
               hit("chr14", 105'625'700, 105'625'994)};   // IGHG4-CH1
    d.part_probabilities = {0.55f, 0.45f};
    d.kind = 'V';  // VDJ class-switch
    d.vdj_class_switch_detected = true;
    d.genomic_distance_bp = 234'000;

    auto r = make_chimeric("read_csr_1", 2'400, d);
    EXPECT_EQ(r.status, AlignmentStatus::ChimericIntraRegion);
    EXPECT_TRUE(r.is_lossless_consistent());
    EXPECT_TRUE(r.chimeric->vdj_class_switch_detected);
}

TEST(AlignmentRecordExt, ChimericInterChromValid) {
    ChimericDetail d;
    d.parts = {hit("chr1", 1'000, 1'500),
               hit("chr22", 50'000'000, 50'000'500)};
    d.part_probabilities = {0.5f, 0.5f};
    d.kind = 'X';

    auto r = make_chimeric("read_trans_1", 2'000, d);
    EXPECT_EQ(r.status, AlignmentStatus::ChimericInterChrom);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, ChimericLessThanTwoPartsFails) {
    ChimericDetail d;
    d.parts = {hit()};
    d.part_probabilities = {1.0f};
    d.kind = 'I';
    AlignmentRecord r;
    r.read_id = "x";
    r.status = AlignmentStatus::ChimericIntraRegion;
    r.chimeric = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, ChimericPartProbsSizeMismatchFails) {
    ChimericDetail d;
    d.parts = {hit(), hit()};
    d.part_probabilities = {0.5f};  // 1 prob, 2 parts
    d.kind = 'I';
    AlignmentRecord r;
    r.read_id = "x";
    r.status = AlignmentStatus::ChimericIntraRegion;
    r.chimeric = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// DarkNovel — cluster-internal anchor, no DB hit
// ===========================================================================

TEST(AlignmentRecordExt, DarkNovelValid) {
    DarkNovelDetail d;
    d.cluster_id = 42;
    d.cluster_anchor_id = "CLUSTER:42:representative";
    d.cluster_size = 12;
    d.cluster_anchor_ids = {"read_x1", "read_x2", "read_x3"};

    auto r = make_dark_novel("read_dark_1", 1'200, d);
    EXPECT_EQ(r.status, AlignmentStatus::DarkNovel);
    EXPECT_TRUE(r.is_lossless_consistent());
}

TEST(AlignmentRecordExt, DarkNovelEmptyClusterAnchorFails) {
    DarkNovelDetail d;
    d.cluster_id = 1;
    // d.cluster_anchor_id intentionally empty
    d.cluster_size = 5;
    AlignmentRecord r;
    r.read_id = "x";
    r.status = AlignmentStatus::DarkNovel;
    r.dark_novel = d;
    EXPECT_FALSE(r.is_lossless_consistent());
}

// ===========================================================================
// TranscriptKind round-trip + Name/Parse helpers
// ===========================================================================

TEST(TranscriptKind, NameRoundTripCoversAllEnumValues) {
    using K = core::TranscriptKind;
    constexpr K all_kinds[] = {
        K::Unknown, K::MatureMrna, K::PreMrna, K::SterileGermline,
        K::CircularRna, K::Lncrna,
        K::Snorna_CDbox, K::Snorna_HacaBox, K::Scarna,
        K::Mirna, K::Pirna, K::Sirna,
        K::Snrna_Major, K::Snrna_Minor,
        K::Trna, K::Rrna, K::Yrna, K::Vaultrna,
        K::Srprna_7sl, K::Rna_7sk, K::TercTerra, K::Rmrp_Rnasep,
        K::Erna, K::TirnaParna, K::Antisense, K::Intergenic,
        K::Mitochondrial, K::RepeatDerived, K::Viral,
        K::Fusion, K::NovelUnclassified,
    };
    for (auto k : all_kinds) {
        const char* name = core::TranscriptKindName(k);
        ASSERT_NE(name, nullptr);
        ASSERT_GT(std::string(name).size(), 0u);
        auto parsed = core::ParseTranscriptKind(name);
        ASSERT_TRUE(parsed.has_value()) << "name='" << name << "'";
        EXPECT_EQ(*parsed, k);
    }
}

TEST(TranscriptKind, ParseUnknownLabelReturnsNullopt) {
    EXPECT_FALSE(core::ParseTranscriptKind("not_a_real_kind").has_value());
    EXPECT_FALSE(core::ParseTranscriptKind("").has_value());
}

TEST(TranscriptKind, ClassPredicates) {
    using K = core::TranscriptKind;
    EXPECT_TRUE(core::IsCodingFamily(K::MatureMrna));
    EXPECT_TRUE(core::IsCodingFamily(K::PreMrna));
    EXPECT_TRUE(core::IsCodingFamily(K::SterileGermline));
    EXPECT_FALSE(core::IsCodingFamily(K::Lncrna));

    EXPECT_TRUE(core::IsSmallRna(K::Mirna));
    EXPECT_TRUE(core::IsSmallRna(K::Pirna));
    EXPECT_TRUE(core::IsSmallRna(K::Sirna));
    EXPECT_TRUE(core::IsSmallRna(K::Snorna_CDbox));
    EXPECT_FALSE(core::IsSmallRna(K::MatureMrna));

    EXPECT_TRUE(core::IsStructuralRna(K::Trna));
    EXPECT_TRUE(core::IsStructuralRna(K::Rrna));
    EXPECT_FALSE(core::IsStructuralRna(K::CircularRna));
}

// ===========================================================================
// AlignmentStatus name table
// ===========================================================================

TEST(AlignmentStatus, NameTableCoversAllValues) {
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::Mapped),               "MAPPED");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::MappedSterile),        "MAPPED_STERILE");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::MappedPreMrna),        "MAPPED_PREMRNA");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::MappedCircular),       "MAPPED_CIRCULAR");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::PartialMappedExact),   "PARTIAL");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::ChimericIntraRegion),  "CHIM_INTRA");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::ChimericInterChrom),   "CHIM_INTER");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::DarkNovel),            "DARK");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::Tentative),            "TENT");
    EXPECT_STREQ(AlignmentStatusName(AlignmentStatus::Unmapped),             "UNMAPPED");
}

TEST(AlignmentStatus, IsAnyMappedCoversCorrectStatuses) {
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::Mapped));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::MappedSterile));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::MappedPreMrna));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::MappedCircular));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::PartialMappedExact));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::ChimericIntraRegion));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::ChimericInterChrom));
    EXPECT_TRUE(IsAnyMapped(AlignmentStatus::DarkNovel));
    EXPECT_FALSE(IsAnyMapped(AlignmentStatus::Tentative));
    EXPECT_FALSE(IsAnyMapped(AlignmentStatus::Unmapped));
}
