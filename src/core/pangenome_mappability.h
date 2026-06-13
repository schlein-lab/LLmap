// LLmap — Pangenome mapping-determinism prior M(pos).
//
// Operator idea (2026-06-13): a placement is more probable when the SAME mapping
// is reproducibly sensible across many pangenome samples, and honestly less
// certain when it is similarly blurry everywhere. This is the MAPPING twin of
// the splice-determinism D(pos): a per-locus mapping reproducibility across
// hundreds of genomes from different populations.
//
//   M(pos) = cross-sample mean mapping confidence at a locus.
//   * high  → reproducibly uniquely-mappable (conserved unique locus) → boost
//   * low   → reproducibly blurry EVERYWHERE (real segdup/TE/paralog property) →
//             honest low confidence, NOT a penalty — the ambiguity is reliable
//   * a sample whose LOCAL confidence is low where M(pos) is high = an ANOMALY
//     (artifact/contamination/sample-specific) — that comparison is done by the
//     consumer (te_spread_mass), not here.
//
// Pure aggregator (stdlib only, lives in core so both the mapping EM and the
// provenance layer can read it without a cross-module dependency). The pangenome
// run feeds ONE per-sample mapping-confidence per locus-window via AddSample;
// PopulationSupport() returns M for a locus. Stored as a binned track (a single
// per-base map would be genome×samples — too large), window default 100 bp since
// mappability is a locus property, not a per-base one.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <string_view>

namespace llmap::core {

class PangenomeMappability {
public:
    explicit PangenomeMappability(std::uint32_t window_bp = 100)
        : window_bp_(window_bp == 0 ? 1 : window_bp) {}

    // Add ONE sample's mapping confidence (in [0,1]) at a locus. Call once per
    // (sample, locus) — the run computes a per-sample-per-window confidence
    // (e.g. mean MAPQ→[0,1] or collapse confidence) and adds it here; M is the
    // mean over those per-sample values.
    void AddSample(std::string_view ref_id, std::uint64_t pos, float confidence);

    // population_support ∈ [0,1] for a locus = M(window). Returns `neutral` for
    // a never-seen window (no pangenome info ⇒ don't bias the EM either way).
    [[nodiscard]] float PopulationSupport(std::string_view ref_id,
                                          std::uint64_t pos,
                                          float neutral = 0.5f) const;

    // Number of samples that contributed to the window covering pos.
    [[nodiscard]] std::uint32_t SampleCount(std::string_view ref_id,
                                            std::uint64_t pos) const;

    // BedGraph-like track I/O: "ref<TAB>start<TAB>end<TAB>M<TAB>n_samples".
    [[nodiscard]] bool Save(const std::string& path) const;
    [[nodiscard]] bool Load(const std::string& path);

    [[nodiscard]] std::uint32_t window_bp() const noexcept { return window_bp_; }
    [[nodiscard]] std::size_t n_windows() const noexcept;

private:
    struct Cell {
        double        sum{0.0};   // Σ per-sample confidences
        std::uint32_t n{0};       // contributing samples
    };
    std::uint32_t window_bp_;
    // ref_id → (window index → cell)
    std::map<std::string, std::map<std::uint64_t, Cell>> track_;
};

}  // namespace llmap::core
