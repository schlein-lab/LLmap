// LLmap — JunctionDb: multi-source splice-junction evidence aggregator.
//
// One JunctionEvidence record per (chrom, donor_pos, acceptor_pos)
// triple. Sources are pulled in from disk separately and the
// JunctionDb merges them. Each source contributes a boolean flag
// (present/absent) and optionally a tissue-frequency map (GTEx only).
//
// Consumers:
//   - The Multi-Signal Fusion engine (Block 4.5) reads JunctionEvidence
//     when computing L_junction for a candidate junction. A junction
//     present in GENCODE + GTEx + ChessDB gets a higher likelihood
//     boost than one only present in GENCODE — but the floor never
//     goes to zero.
//   - The k-mer-index builder (Block 3) can opt to skip junction-
//     spanning k-mers across junctions that NO source supports
//     (cuts down index size). Off by default to preserve lossless
//     semantics; the option is `--strict-junctions`.
//
// Storage: nested unordered_map<chrom, unordered_map<(donor,acceptor),
// JunctionEvidence>>. For human autosomes + GENCODE v46 + GTEx v8 the
// expected memory footprint is ~80 MB; sufficient. If the catalog
// grows past ~200 M junctions a fast-flat-hashmap migration is easy.

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::annot {

struct JunctionKey {
    std::int64_t donor;      ///< genomic 0-based
    std::int64_t acceptor;   ///< genomic 0-based; donor>acceptor for back-splice
    bool operator==(const JunctionKey& o) const noexcept {
        return donor == o.donor && acceptor == o.acceptor;
    }
};

struct JunctionKeyHash {
    std::size_t operator()(const JunctionKey& k) const noexcept {
        // Mix donor + acceptor with a 64-bit avalanche so adjacent pairs
        // don't cluster in the bucket array. Cheap; called per lookup.
        const std::uint64_t a = static_cast<std::uint64_t>(k.donor);
        const std::uint64_t b = static_cast<std::uint64_t>(k.acceptor);
        std::uint64_t x = a * 0x9E3779B97F4A7C15ULL + b;
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;
        return static_cast<std::size_t>(x);
    }
};

struct JunctionEvidence {
    bool in_gencode{false};
    bool in_mane{false};
    bool in_refseq{false};
    bool in_chessdb{false};
    bool in_gtex{false};
    bool in_circ_db{false};       ///< CIRCpedia / circBase
    std::uint32_t gtex_sample_count{0};
    /// Per-tissue freq: tissue label → freq ∈ [0,1].
    std::unordered_map<std::string, float> gtex_tissue_freq;
};

class JunctionDb {
public:
    // ----- Loaders (each best-effort, returns true on parse success) -----
    bool LoadGencode(const std::filesystem::path& gff);
    bool LoadGtexJunctions(const std::filesystem::path& bed);
    bool LoadCircRnaDb(const std::filesystem::path& bed);
    bool LoadChessDbJunctions(const std::filesystem::path& gtf);

    // ----- Inspection -----------------------------------------------------

    /// Direct lookup. Returns evidence struct with all flags = false
    /// when the junction isn't recorded.
    JunctionEvidence Lookup(std::string_view chrom,
                            std::int64_t donor,
                            std::int64_t acceptor) const;

    /// True iff at least one source has this junction.
    bool HasAnyEvidence(std::string_view chrom,
                         std::int64_t donor,
                         std::int64_t acceptor) const;

    /// Aggregate count across all chroms; used by the lossless summary.
    std::size_t TotalJunctions() const;

    /// Per-chrom count (telemetry).
    std::size_t JunctionsOnChrom(std::string_view chrom) const;

private:
    // chrom → ( (donor,acceptor) → evidence ).
    std::unordered_map<std::string,
        std::unordered_map<JunctionKey, JunctionEvidence,
                           JunctionKeyHash>> by_chrom_;

    /// Internal: get-or-create an evidence record for mutation.
    JunctionEvidence& GetOrInsert(std::string_view chrom,
                                    std::int64_t donor,
                                    std::int64_t acceptor);
};

}  // namespace llmap::annot
