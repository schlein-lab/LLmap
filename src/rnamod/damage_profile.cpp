// LLmap — Damage / chemistry / RNA-editing substitution classifier.

#include "rnamod/damage_profile.h"

#include <cctype>
#include <string>

namespace llmap::rnamod {

namespace {

char Up(char c) noexcept {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

char Comp(char c) noexcept {
    switch (Up(c)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

// Reverse-complement a short context (uppercased).
std::string RevComp(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) out.push_back(Comp(*it));
    return out;
}

}  // namespace

const char* SubstitutionProvenanceName(SubstitutionProvenance p) noexcept {
    switch (p) {
        case SubstitutionProvenance::VariantCandidate:  return "none";
        case SubstitutionProvenance::Damage8oxoG:       return "dmg:8oxoG";
        case SubstitutionProvenance::DamageFfpeDeam:    return "dmg:ffpe";
        case SubstitutionProvenance::DamageAncientDeam: return "dmg:deam";
        case SubstitutionProvenance::EditAdarA2I:       return "edit:adar";
        case SubstitutionProvenance::EditApobecC2U:     return "edit:apobec";
    }
    return "none";
}

ProvenanceCall ClassifySubstitution(const SubstitutionContext& c) {
    const char ref = Up(c.ref);
    const char alt = Up(c.alt);
    if (ref == 'N' || alt == 'N' || ref == alt) {
        return ProvenanceCall{};  // no substitution / unknown → VariantCandidate
    }

    // RNA editing is a biological process on the sense (transcript) strand —
    // evaluate the substitution as the molecule saw it. For a '-'-strand read
    // the sense substitution is the reverse-complement of the genomic one.
    if (c.is_rna) {
        const char s_ref = c.read_reverse ? Comp(ref) : ref;
        const char s_alt = c.read_reverse ? Comp(alt) : alt;
        // ADAR A-to-I: adenosine deaminated to inosine, read as guanosine.
        if (s_ref == 'A' && s_alt == 'G') {
            return ProvenanceCall{SubstitutionProvenance::EditAdarA2I, 0.6f, false};
        }
        // APOBEC C-to-U: prefers a 5' thymine (TC dinucleotide) on the sense
        // strand. ctx5 is the +-strand 5 bp centred on the site (index 2); take
        // the sense view so the 5' neighbour is well-defined on either strand.
        if (s_ref == 'C' && s_alt == 'T' && c.ctx5.size() == 5) {
            const std::string sense =
                c.read_reverse ? RevComp(c.ctx5) : std::string(c.ctx5);
            if (sense.size() == 5 && sense[1] == 'T' && sense[2] == 'C') {
                return ProvenanceCall{SubstitutionProvenance::EditApobecC2U, 0.6f,
                                      false};
            }
        }
    }

    // DNA damage signatures (genomic substitution type; strand-asymmetric — the
    // decisive bias is a site-level aggregate, refined in the provenance layer).
    const bool g_to_t = (ref == 'G' && alt == 'T') || (ref == 'C' && alt == 'A');
    if (g_to_t) {
        return ProvenanceCall{SubstitutionProvenance::Damage8oxoG, 0.5f, true};
    }
    const bool deam = (ref == 'C' && alt == 'T') || (ref == 'G' && alt == 'A');
    if (deam) {
        // Cytosine deamination concentrated at read ends → ancient/degraded DNA.
        if (c.dist_from_read_end <= 5) {
            return ProvenanceCall{SubstitutionProvenance::DamageAncientDeam, 0.6f,
                                  true};
        }
        return ProvenanceCall{SubstitutionProvenance::DamageFfpeDeam, 0.4f, true};
    }

    return ProvenanceCall{};  // transition/transversion with no artifact signature
}

}  // namespace llmap::rnamod
