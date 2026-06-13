// LLmap — chimera_provenance (class D) tests.
//
// Artifact chimera (cross-chrom / intra-region) → is_chimera (Layer-1 `chim`);
// VDJ class-switch → is_vdj_recombination (Layer-3 `vdj`, NOT a chimera artifact);
// a single part → neither. Reuses the Block-7 detector via the wrapper.

#include "provenance/chimera_provenance.h"
#include "chimera/chimera_detector.h"

#include <gtest/gtest.h>

#include <vector>

namespace {

using llmap::chimera::AlignedPart;
using llmap::chimera::VdjLocusMask;
using llmap::provenance::ClassifyChimera;

AlignedPart MakePart(std::string chrom, std::uint64_t rs, std::uint64_t re,
                     std::uint32_t qo, std::uint32_t ql,
                     bool switch_region = false) {
    AlignedPart p;
    p.ref_chrom = std::move(chrom);
    p.ref_start = rs;
    p.ref_end = re;
    p.read_offset = qo;
    p.read_length = ql;
    p.score = 100.0f;
    p.is_switch_region = switch_region;
    return p;
}

VdjLocusMask Mask() {
    VdjLocusMask m;
    m.LoadGrch38Defaults();
    return m;
}

TEST(ChimeraProvenance, CrossChromIsArtifactChimera) {
    const auto m = Mask();
    std::vector<AlignedPart> parts = {
        MakePart("chr1", 1000, 1500, 0, 500),
        MakePart("chr22", 5000, 5500, 500, 500),
    };
    const auto c = ClassifyChimera(parts, m);
    EXPECT_TRUE(c.is_chimera);
    EXPECT_FALSE(c.is_vdj_recombination);
    EXPECT_EQ(c.kind, 'X');
}

TEST(ChimeraProvenance, VdjClassSwitchIsNotArtifact) {
    const auto m = Mask();
    std::vector<AlignedPart> parts = {
        MakePart("chr14", 105'860'000, 105'860'400, 0, 400),       // IGHM-CH1
        MakePart("chr14", 105'625'700, 105'625'994, 400, 294, true),  // S-region
    };
    const auto c = ClassifyChimera(parts, m);
    // Real immune biology → Layer-3 `vdj`, never bucketed as a chimera artifact.
    EXPECT_TRUE(c.is_vdj_recombination);
    EXPECT_FALSE(c.is_chimera);
    EXPECT_EQ(c.kind, 'V');
}

TEST(ChimeraProvenance, SinglePartIsNeither) {
    const auto m = Mask();
    std::vector<AlignedPart> parts = {MakePart("chr14", 1000, 1500, 0, 500)};
    const auto c = ClassifyChimera(parts, m);
    EXPECT_FALSE(c.is_chimera);
    EXPECT_FALSE(c.is_vdj_recombination);
}

TEST(ChimeraProvenance, IntraRegionIsSvSpanningNotArtifact) {
    // Intra-chromosomal split = a read crossing a real structural breakpoint →
    // Layer-3 `bio:sv` (host biology), NOT a Layer-1 artifact chimera. (Matched
    // HG002: HiFi shows these, Illumina ~0 — the opposite of an artifact signal.)
    const auto m = Mask();
    std::vector<AlignedPart> parts = {
        MakePart("chr1", 100'000, 100'500, 0, 500),
        MakePart("chr1", 200'000, 200'500, 500, 500),  // 100 kb gap
    };
    const auto c = ClassifyChimera(parts, m);
    EXPECT_TRUE(c.is_sv_spanning);
    EXPECT_FALSE(c.is_chimera);
    EXPECT_FALSE(c.is_vdj_recombination);
    EXPECT_EQ(c.kind, 'I');
}

}  // namespace
