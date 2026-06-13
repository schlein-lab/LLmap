// LLmap — Read provenance: a THREE-LAYER model (Contamination / Provenance Mode).
//
// Layer 1 — read-ORIGIN partition (`ProvenanceClass`). Every read is in EXACTLY
//   one origin class, so the lossless Σ-invariant holds: Σ over classes ==
//   N_input. `Host` is the regular read. This is the partition the variant
//   caller filters on.
// Layer 2 — per-SITE base-artefact overlay (8-oxoG / FFPE / deamination /
//   bisulfite / A→I / C→U). A *position* property of an otherwise-Host read; it
//   does NOT move the read out of its origin class. Lives in a separate per-base
//   track (Agent 1's rnamod/damage_profile), NOT here.
// Layer 3 — per-READ biology overlay (`BioConfounder` bitflags). V(D)J / SHM /
//   class-switch / gene-conversion / true mtDNA-heteroplasmy / mosaicism are
//   REAL host biology that masquerades as variation — they co-occur with `Host`
//   and must NEVER be bucketed as non-host (that would strip real signal).
//
// Open-ended (Agent 1): the enum stays at mechanism-FAMILY granularity; specifics
// live in free-form `detail` strings (`exo:mouse`, `ref:ancestry`, `numt`, …) so
// the taxonomy grows without an ABI break.
// See docs/design/llmap_provenance_axes.md.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::provenance {

// ===========================================================================
// Layer 1 — read-origin partition (Σ == N_input).
// ===========================================================================
enum class ProvenanceClass : std::uint8_t {
    Host = 0,        // regular: clean host-germline observation

    // A. competing reference — origin is a different sequence (origin buckets)
    Exogenous,       // foreign organism: EBV, Mycoplasma, kitome, PhiX, vector
    CrossIndividual, // second mammalian genome: xenograft mouse, HeLa/cell-line
                     // cross-contam, sample-swap, transplant/materno-fetal
    Paralog,         // host, wrong paralog/SD copy
    Numt,            // nuclear-mtDNA segment (vs TRUE heteroplasmy → BioConfounder)
    Pseudogene,      // processed pseudogene (GBAP1/GBA1, PMS2CL, SMN2)
    Rdna,            // rDNA / satellite / decoy / unplaced contig
    MobileElement,   // somatic L1/Alu/SVA, HERV, ecDNA/eccDNA

    // E. the reference is wrong, not the read (still a partition member)
    RefArtefact,     // assembly collapse/gap, ancestry-mismatch vs GRCh38,
                     // liftover flip, reference-carries-minor-allele

    // C. structural split (read = two molecules)
    Chimera,         // ligation/PCR chimera, ONT fused/concatemer, RT switch

    // D. technical cross-sample / multiplicity
    CrossSample,     // index-hopping, within-run barcode collision
    Multiplicity,    // PCR/optical duplicate, WGA/MDA, PCR-jackpot, UMI collision

    Count
};

[[nodiscard]] const char* ProvenanceClassTag(ProvenanceClass c) noexcept;
[[nodiscard]] std::optional<ProvenanceClass> ParseProvenanceClass(
    std::string_view tag) noexcept;
// True for the competing-reference family (A): these become provenance buckets
// in the WaveCollapse EM. Host / RefArtefact / Chimera / CrossSample /
// Multiplicity are partition members but not "another sequence" buckets.
[[nodiscard]] bool IsOriginBucket(ProvenanceClass c) noexcept;

// ===========================================================================
// Layer 3 — per-read biology overlay (bitflags; co-occurs with Host, Σ-safe).
// ===========================================================================
enum class BioConfounder : std::uint16_t {
    None           = 0,
    Vdj            = 1u << 0,  // V(D)J recombination (IGH/TCR)
    Shm            = 1u << 1,  // somatic hypermutation (AID)
    ClassSwitch    = 1u << 2,  // IGH class-switch recombination
    GeneConversion = 1u << 3,  // paralog↔paralog real exchange
    MtHeteroplasmy = 1u << 4,  // TRUE mtDNA heteroplasmy (not a NUMT artefact)
    Mosaicism      = 1u << 5,  // tissue/developmental mosaic, LOY, mosaic CNV
    Imprinting     = 1u << 6,  // imprinting / X-inactivation allelic skew (RNA)
    SvSpanning     = 1u << 7,  // intra-chrom split = real SV breakpoint (not a
                               // ligation chimera) — Class-D `bio:sv`
};
[[nodiscard]] const char* BioConfounderTag(BioConfounder f) noexcept;

// ===========================================================================
// One read's resolved provenance.
// ===========================================================================
struct ReadProvenance {
    // Layer 1 — exactly one origin class (the partition; `Host` = regular).
    ProvenanceClass origin{ProvenanceClass::Host};
    float           posterior{1.0f};       // collapse confidence (retained)
    std::uint32_t   aligned_bases{0};      // for per-class base accounting
    std::string     detail;                // origin specifics: "exo:mouse", "numt"

    // Layer 3 — per-read biology overlay bitmask (OR of BioConfounder flags).
    // Co-occurs with origin (usually Host); does NOT affect the Σ-partition.
    std::uint16_t   bioconfounder{0};
};

}  // namespace llmap::provenance
