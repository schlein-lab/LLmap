// LLmap — Damage / chemistry / RNA-editing substitution classifier.
//
// Provenance-layer detector class B (contamination/artifact mode): given a
// single observed substitution (reference base → observed base) plus strand,
// flanking context and read-end distance, classify it as an artifact / editing
// signature vs a genuine variant candidate. These are the sources that fake
// low-VAF SNVs (operator taxonomy B): 8-oxoG (G>T, acoustic-shearing, FoxoG
// strand bias), FFPE deamination (C>T / G>A), ancient/degraded-DNA cytosine
// deamination (C>T concentrated at read ends), ADAR A-to-I editing (reads as
// A>G, dsRNA/Alu context), APOBEC C-to-U editing (C>T in TC context).
//
// Per-substitution evidence only — flags, never filters (lossless doctrine).
// The site-level strand-bias aggregate (the decisive 8-oxoG/FFPE signal) is
// refined downstream in the provenance layer; here we emit the candidate class,
// a base confidence, and whether the signature is strand-asymmetric. Reuses the
// sequence-context-prior idea from ModificationDb (same module family).

#pragma once

#include <cstdint>
#include <string_view>

namespace llmap::rnamod {

enum class SubstitutionProvenance : std::uint8_t {
    VariantCandidate = 0,  // no artifact/editing signature — treat as real
    Damage8oxoG,           // G>T (+) / C>A (−); oxidative, strand-biased
    DamageFfpeDeam,        // C>T / G>A; formalin deamination, strand-biased
    DamageAncientDeam,     // C>T / G>A concentrated at read ends (aDNA)
    EditAdarA2I,           // A>G (sense); ADAR A-to-I, dsRNA/Alu context
    EditApobecC2U,         // C>T in TC context; APOBEC C-to-U (RNA)
};

[[nodiscard]] const char* SubstitutionProvenanceName(
    SubstitutionProvenance p) noexcept;

struct SubstitutionContext {
    char ref{'N'};                 // reference base, uppercase, + strand
    char alt{'N'};                 // observed base, uppercase, + strand
    bool read_reverse{false};      // read aligned to the '-' strand
    std::string_view ctx5;         // 5 bp + strand reference context, site centred (index 2)
    std::uint32_t dist_from_read_end{1000};  // min(offset, len-offset) of the site in the read
    bool is_rna{false};            // RNA input → editing classes enabled
};

struct ProvenanceCall {
    SubstitutionProvenance kind{SubstitutionProvenance::VariantCandidate};
    float confidence{0.0f};        // [0,1] base confidence for the call
    bool strand_biased{false};     // signature is expected strand-asymmetric
};

// Classify a single substitution. Disambiguation when several signatures match
// the same base change: read-end proximity → ancient deamination; RNA + TC
// context → APOBEC; RNA + A>G → ADAR; otherwise the strand-biased DNA-damage
// classes (8-oxoG for G>T/C>A, FFPE for C>T/G>A). No signature → VariantCandidate.
[[nodiscard]] ProvenanceCall ClassifySubstitution(const SubstitutionContext& c);

}  // namespace llmap::rnamod
