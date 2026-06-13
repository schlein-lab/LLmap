// LLmap — Mapping-confusion provenance detector (Layer-1 origin classes).
//
// Provenance-layer detector for reads that are placed on the WRONG genomic copy
// — the operator's class C (mapping/reference artifacts that fake low-VAF
// variants). These are Layer-1 ORIGIN classes (partition members, Σ-invariant):
// the read's true origin is a different copy than where it landed.
//
//   Paralog    — high-identity paralog / segmental-duplication mis-placement
//                (PSV posterior is ambiguous; LLmap's IGHG4 / retina-segdup case)
//   Numt       — nuclear mtDNA segment read placed on (or stolen from) the mt
//                contig → the classic fake "heteroplasmy". The ARTIFACT case.
//   Pseudogene — processed (intronless) pseudogene vs its intron-bearing parent
//                (GBAP1/GBA1, PMS2/PMS2CL, SMN1/2); in transcript data the
//                pseudogene fakes intronless "expression".
//   Rdna       — rDNA / satellite / HOR-array read that cannot be uniquely placed
//
// Pure classifier over an explicit evidence struct (like rnamod::damage_profile)
// — no dependency on the partition enum; cmd/EM layer maps the result to the
// ProvenanceClass + the PV detail string. Evidence the pipeline cannot supply
// defaults to "no signal", so that sub-class simply returns None.
//
// NOTE the deliberate split from REAL mt heteroplasmy: a genuine mt-derived
// variant is host biology (Layer-3 bioconfounder `mthet`), NOT a Numt artifact.
// IsRealMtHeteroplasmy() exposes that decision over the same evidence so the
// bioconfounder layer and this Layer-1 classifier never double-count an mt read.

#pragma once

#include <cstdint>

namespace llmap::provenance {

enum class MappingConfusion : std::uint8_t {
    None = 0,    // uniquely / correctly placed — defer to other detectors (host)
    Paralog,
    Numt,
    Pseudogene,
    Rdna,
};

// PV detail string for the Layer-1 partition (e.g. "para", "numt", "pseudo",
// "rdna"); "none" for None. Matches the agreed provenance tag scheme.
[[nodiscard]] const char* MappingConfusionTag(MappingConfusion m) noexcept;

struct MappingConfusionEvidence {
    // --- Paralog (PSV) ---
    bool  has_psv{false};
    float psv_posterior{1.0f};        // posterior of the assigned copy; low ⇒ ambiguous

    // --- mt / NUMT discrimination (only when the read sits at an mt-homologous
    //     locus, i.e. mt sequence also present as a nuclear NUMT) ---
    bool  at_mt_homologous_locus{false};
    bool  placed_on_mt_contig{false};  // where the aligner actually put it
    float identity_to_mt{0.0f};        // read vs mitochondrial reference
    float identity_to_nuclear_numt{0.0f};  // read vs best nuclear NUMT copy

    // --- Pseudogene ---
    bool at_pseudogene_parent_locus{false};  // a known processed-pseudogene/parent pair
    bool read_is_spliced{false};             // has intron (N) ops → parent; intronless ⇒ pseudogene

    // --- repeat array ---
    bool         in_repeat_array{false};   // rDNA / α-satellite / HOR array
    std::uint32_t mapq{60};

    // Thresholds (overridable).
    float paralog_posterior_threshold{0.90f};
    float mt_identity_margin{0.01f};       // nuclear must beat mt by this to call Numt
};

struct MappingConfusionCall {
    MappingConfusion kind{MappingConfusion::None};
    float confidence{0.0f};
};

// Classify the read's Layer-1 origin confusion from the evidence. Priority:
// pseudogene (structural, unambiguous) → numt (mt-homologous + nuclear wins) →
// paralog (PSV-ambiguous) → rdna (repeat + low MAPQ) → None.
[[nodiscard]] MappingConfusionCall ClassifyMappingConfusion(
    const MappingConfusionEvidence& e);

// True iff the read is best explained as a GENUINE mitochondrial read carrying a
// real (heteroplasmic) variant — host biology, Layer-3 `mthet`, NOT a Numt
// artifact. Mutually exclusive with a Numt call over the same evidence.
[[nodiscard]] bool IsRealMtHeteroplasmy(const MappingConfusionEvidence& e);

}  // namespace llmap::provenance
