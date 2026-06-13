// LLmap — Contamination / provenance spectrum aggregator.
//
// Accumulates per-read provenance into the sample's spectrum: per-class read
// count, read fraction, base count, and mean collapse posterior. Enforces the
// lossless Σ-invariant — every input read is counted in exactly one class, so
// Σ over classes == N_input. A violation means a read was dropped = a bug.
//
// Output: a flat spectrum table (`class  n_reads  fraction  bases  mean_post`),
// directly plottable as the sample's contamination profile ("0.8% PhiX, 2.1%
// NUMT-collapse, …"). TSV here; the BAM/Parquet writers consume the same struct.

#pragma once

#include "provenance/provenance_class.h"

#include <array>
#include <cstdint>
#include <string>

namespace llmap::provenance {

struct ClassStat {
    std::uint64_t n_reads{0};
    double        fraction{0.0};       // n_reads / total
    std::uint64_t bases{0};
    double        mean_posterior{0.0};
};

class ContaminationSpectrum {
public:
    // Accumulate one resolved read: its Layer-1 origin class (partition) plus,
    // separately, its Layer-3 bioconfounder overlay flags (which co-occur with
    // the origin — usually Host — and do NOT enter the partition / Σ).
    void Add(const ReadProvenance& rp) noexcept;

    [[nodiscard]] std::uint64_t TotalReads() const noexcept { return total_; }
    [[nodiscard]] ClassStat Stat(ProvenanceClass c) const noexcept;
    // Overlay: reads carrying a given bioconfounder flag (can overlap Host).
    [[nodiscard]] std::uint64_t BioConfounderReads(BioConfounder f) const noexcept;

    // Σ-invariant: every input read landed in exactly one ORIGIN class. The
    // overlay flags are orthogonal and explicitly NOT summed here.
    [[nodiscard]] bool CheckLossless(std::uint64_t n_input) const noexcept {
        return total_ == n_input;
    }

    // `class<TAB>n_reads<TAB>fraction<TAB>bases<TAB>mean_posterior`, header + one
    // row per non-empty class (Host first). Returns false on write error.
    [[nodiscard]] bool WriteTsv(const std::string& path) const;
    [[nodiscard]] std::string ToString() const;

private:
    static constexpr std::size_t kN = static_cast<std::size_t>(ProvenanceClass::Count);
    std::array<std::uint64_t, kN> n_reads_{};
    std::array<std::uint64_t, kN> bases_{};
    std::array<double, kN>        post_sum_{};
    std::uint64_t total_{0};

    // Layer-3 overlay: count of reads carrying each bioconfounder bit (1<<0..6).
    static constexpr std::size_t kBio = 8;
    std::array<std::uint64_t, kBio> bio_{};
};

}  // namespace llmap::provenance
