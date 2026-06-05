// LLmap — junction_hunter: reference-anchored routing gate for
// large genome-scale junction-hunt inputs.
//
// Loads a minimizer index over GRCh38 (or any reference) and an
// AnnotationIndex of NAHR-pair LCR arms. Per-read, it does a cheap
// seed-cluster lookup against the reference, then decides whether the
// read needs the expensive multi-tier NAHR cache scan or can be
// classified as boring (mapped well outside any NAHR pair). About 99 %
// of reads from an HPRC-grade assembly are SKIP, the remaining 1 % go
// through ROUTE_TO_NAHR.
//
// Output verdict is one of:
//   Skip          — confidently mapped to a non-NAHR ref region
//   RouteToNahr   — mapped or partially mapped to an NAHR pair arm
//                   (pair_id field is populated)
//   Unmapped      — too few seeds to call (rare, random or novel)

#pragma once

#include "annot/interval_tree.h"
#include "classical/minimizer_index.h"

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>

namespace llmap::junction_hunter {

struct GateResult {
    enum Verdict : std::uint8_t { Skip = 0, RouteToNahr = 1, Unmapped = 2 };
    Verdict       verdict{Unmapped};
    std::uint32_t ref_id{0};
    std::uint32_t ref_start{0};       ///< best seed bucket start (bp)
    std::uint32_t ref_end{0};         ///< best seed bucket end (bp)
    std::uint16_t n_seeds_in_bucket{0};
    std::string   pair_id;            ///< populated only on RouteToNahr
    std::string   arm;                ///< "up" or "dn" — only on RouteToNahr
};

struct ReferenceGateConfig {
    std::uint16_t min_seeds{8};       ///< minimum seeds in best bucket
    std::uint32_t cluster_window{50000};  ///< bp window for seed bucketing
};

class ReferenceGate {
public:
    ReferenceGate() = default;

    /// Load the minimizer index (.llmi) and the NAHR-arm BED.
    /// BED format: chrom\\tstart\\tend\\tpair_id\\tarm
    bool Load(const std::string& llmi_path,
              const std::string& nahr_bed_path);

    /// Last error from Load() if it returned false.
    const std::string& LastError() const noexcept { return last_error_; }

    /// Classify one read. Pass the raw nucleotide sequence (no '>' header).
    GateResult Classify(std::string_view read_seq) const;

    /// Number of NAHR-arm intervals indexed; 0 = Load() not run yet.
    std::size_t NumNahrIntervals() const noexcept;

    /// Map a contig name (as in the reference FASTA) to its ref_id, or
    /// std::numeric_limits<uint32_t>::max() if unknown.
    std::uint32_t LookupRefId(std::string_view chrom) const noexcept;

    void SetConfig(const ReferenceGateConfig& cfg) noexcept { cfg_ = cfg; }
    const ReferenceGateConfig& GetConfig() const noexcept { return cfg_; }

private:
    std::unique_ptr<classical::MinimizerIndex>    idx_;
    annot::AnnotationIndex                        nahr_iv_;
    std::unordered_map<std::string, std::uint32_t> chrom_to_ref_id_;
    ReferenceGateConfig                           cfg_;
    std::string                                   last_error_;
};

}  // namespace llmap::junction_hunter
