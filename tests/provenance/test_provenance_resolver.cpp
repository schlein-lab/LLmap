// Unit tests for the provenance resolver (evidence → ReadProvenance).

#include "provenance/provenance_resolver.h"

#include <gtest/gtest.h>

namespace llmap::provenance {
namespace {

ReadEvidence HostFit(float host_post) {
    ReadEvidence e;
    e.host_posterior = host_post;
    e.aligned_bases = 150;
    return e;
}

TEST(Resolver, NoEvidenceStaysHost) {
    auto e = HostFit(1.0f);
    const auto p = ResolveProvenance(e);
    EXPECT_EQ(p.origin, ProvenanceClass::Host);
    EXPECT_EQ(p.aligned_bases, 150u);
}

TEST(Resolver, ConfidentParalogBeatsHost) {
    auto e = HostFit(0.40f);                          // poor host fit
    e.mapping = {MappingConfusion::Paralog, 0.95f};   // strong paralog
    const auto p = ResolveProvenance(e);
    EXPECT_EQ(p.origin, ProvenanceClass::Paralog);
    EXPECT_EQ(p.detail, "para");
}

TEST(Resolver, WeakConfounderStaysHostConservative) {
    auto e = HostFit(0.80f);
    e.mapping = {MappingConfusion::Numt, 0.85f};      // only 0.05 over host < margin 0.15
    const auto p = ResolveProvenance(e);
    EXPECT_EQ(p.origin, ProvenanceClass::Host);        // not stripped
}

TEST(Resolver, ExogenousClaimsWhenStrong) {
    auto e = HostFit(0.30f);
    e.exogenous_posterior = 0.97f;
    e.exogenous_detail = "exo:phix";
    const auto p = ResolveProvenance(e);
    EXPECT_EQ(p.origin, ProvenanceClass::Exogenous);
    EXPECT_EQ(p.detail, "exo:phix");
}

TEST(Resolver, ChimeraAndDuplicateOnlyWhenHostLike) {
    // Chimera flag on an otherwise host-like read → Chimera.
    auto e = HostFit(0.95f);
    e.is_chimera = true;
    EXPECT_EQ(ResolveProvenance(e).origin, ProvenanceClass::Chimera);

    // Duplicate on host-like → Multiplicity.
    auto d = HostFit(0.95f);
    d.is_duplicate = true;
    EXPECT_EQ(ResolveProvenance(d).origin, ProvenanceClass::Multiplicity);

    // But a strong origin confounder takes precedence over the chimera flag.
    auto c = HostFit(0.20f);
    c.is_chimera = true;
    c.exogenous_posterior = 0.99f;
    c.exogenous_detail = "exo:ebv";
    EXPECT_EQ(ResolveProvenance(c).origin, ProvenanceClass::Exogenous);
}

TEST(Resolver, BioConfounderCarriedThroughRegardlessOfOrigin) {
    const std::uint16_t vdj = static_cast<std::uint16_t>(BioConfounder::Vdj);
    // Host read with VDJ overlay → stays Host, overlay preserved.
    auto h = HostFit(1.0f);
    h.bioconfounder = vdj;
    EXPECT_EQ(ResolveProvenance(h).origin, ProvenanceClass::Host);
    EXPECT_EQ(ResolveProvenance(h).bioconfounder, vdj);
    // Even when a confounder claims the read, the overlay rides along.
    auto p = HostFit(0.30f);
    p.bioconfounder = vdj;
    p.mapping = {MappingConfusion::Paralog, 0.95f};
    EXPECT_EQ(ResolveProvenance(p).origin, ProvenanceClass::Paralog);
    EXPECT_EQ(ResolveProvenance(p).bioconfounder, vdj);
}

}  // namespace
}  // namespace llmap::provenance
