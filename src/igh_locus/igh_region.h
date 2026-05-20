// LLmap — IGH locus genomic region definitions.
//
// The immunoglobulin heavy-chain constant region (IGHC) at 14q32.33 is a
// tandem array of near-identical paralogs (IGHG1/2/3/4, IGHM, IGHA1/2, IGHE,
// IGHD). Standard seed-chain-extend mappers collapse reads from one paralog
// onto whichever copy wins the chain score, so a post-hoc re-sort needs to know
// (a) whether an alignment falls inside this locus and (b) — for transcript
// references — whether a target name is an IGH copy at all.
//
// Header-only: just a small interval table plus membership tests.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::igh_locus {

// One genomic interval that belongs to the IGH constant region.
struct IghInterval {
    std::string contig;       // e.g. "chr14" or "14"
    std::uint64_t start{0};   // 0-based, inclusive
    std::uint64_t end{0};     // exclusive
};

// A small set of intervals covering the IGHC cluster on common references.
// Genomic re-sort is gated on these so non-IGH reads are never touched.
class IghRegion {
public:
    // Default table: GRCh38 + CHM13 IGHC clusters (constant-region span,
    // generous padding). Contig matching is suffix-based so "chr14", "14",
    // and "...#chr14" assembly-style names all resolve.
    static IghRegion Default() {
        IghRegion r;
        // GRCh38 primary assembly, IGHC cluster ~105.55-105.75 Mbp (minus strand).
        r.intervals_.push_back({"chr14", 105'550'000ULL, 105'760'000ULL});
        // T2T-CHM13 v2.0, IGHC cluster ~99.70-99.90 Mbp.
        r.intervals_.push_back({"chr14", 99'700'000ULL, 99'920'000ULL});
        return r;
    }

    void Add(IghInterval iv) { intervals_.push_back(std::move(iv)); }

    [[nodiscard]] bool Empty() const { return intervals_.empty(); }

    // True iff (target_id, pos) falls inside any IGH interval. Contig names are
    // matched by suffix to tolerate assembly-prefixed names like "HG002#1#chr14".
    [[nodiscard]] bool Contains(std::string_view target_id,
                                std::uint64_t pos) const {
        for (const auto& iv : intervals_) {
            if (ContigMatches(target_id, iv.contig) &&
                pos >= iv.start && pos < iv.end) {
                return true;
            }
        }
        return false;
    }

private:
    static bool ContigMatches(std::string_view target_id,
                              std::string_view contig) {
        if (target_id == contig) return true;
        // suffix match: target endswith ("#"|"_"|...) + contig, or endswith contig
        if (target_id.size() > contig.size() &&
            target_id.compare(target_id.size() - contig.size(),
                              contig.size(), contig) == 0) {
            const char sep = target_id[target_id.size() - contig.size() - 1];
            return sep == '#' || sep == '_' || sep == '.' || sep == '|' ||
                   sep == '/' || sep == ':';
        }
        return false;
    }

    std::vector<IghInterval> intervals_;
};

}  // namespace llmap::igh_locus
