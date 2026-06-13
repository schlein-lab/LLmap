// LLmap — Transposable-element (TE) family catalog (Block 2, Provenance/TE).
//
// Genome-intrinsic TE annotation: which reference positions belong to which TE
// family/subfamily (Alu → AluY/AluYa5/AluSx…, LINE-1 → L1HS/L1PA…, SVA, HERV).
// ~50% of the genome is TE (~1.1M Alu, ~500k L1), so a read landing in a TE is
// the rule, not a corner case. The matched HG002 test showed ~2–6.5% of reads
// are repeat-ambiguous — but raw MAPQ measures that aligner-dependently. This
// catalog replaces the raw-MAPQ proxy with a genome-intrinsic ground truth:
//   * Lookup(chrom,pos) → the TE family/subfamily covering that position,
//   * feeding (a) the `mei` Layer-1 classification (a read confidently inside a
//     TE family is ambiguous across its copies, not uniquely placed), and
//   * (b) the WaveCollapse spread-mass: keep probability mass over the family's
//     copies instead of faking a unique hit; collapse via a unique flank.
//
// Dependency-light (stdlib only). Loadable from a RepeatMasker .out / Dfam BED
// (the production source); a tiny built-in starter lets unit tests + smoke runs
// work without the multi-GB annotation staged.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::provenance {

enum class TeClass : std::uint8_t { Alu, Line1, Sva, Herv, OtherRepeat };

[[nodiscard]] const char* TeClassTag(TeClass c) noexcept;  // "alu","l1","sva","herv","rep"

struct TeLocus {
    std::string   chrom;
    std::uint64_t start{0};       // half-open [start, end)
    std::uint64_t end{0};
    TeClass       cls{TeClass::OtherRepeat};
    std::string   subfamily;      // e.g. "AluYa5", "L1HS", "SVA_E"
    float         divergence{0.0f};  // % from consensus (age proxy; young = ambiguous)
    char          strand{'+'};
};

class TeFamilyCatalog {
public:
    // Tiny in-genome starter (a few representative Alu/L1/SVA loci, GRCh38
    // approximate, divergence-flagged); replace with LoadRepeatMasker/LoadBed.
    void LoadBuiltinStarter();

    // RepeatMasker .out or a BED-like TSV: chrom,start,end,class,subfamily,
    // divergence,strand. Returns false on read error.
    [[nodiscard]] bool LoadBed(const std::string& path);

    // The TE locus covering (chrom,pos), or nullptr if the position is unique.
    [[nodiscard]] const TeLocus* Lookup(std::string_view chrom,
                                        std::uint64_t pos) const;

    [[nodiscard]] std::size_t Size() const noexcept { return loci_.size(); }

private:
    std::vector<TeLocus> loci_;   // sorted by (chrom,start) for binary search
};

}  // namespace llmap::provenance
