// LLmap — Splice-site PWMs + branch-point + polypyrimidine scoring.
//
// Textbook (Burge & Sharp 1997, Sheth et al 2006, Roca et al 2013, Pangolin
// 2023) splice-site frequencies, hand-encoded as position-weight matrices.
// Used at two stages:
//   (a) When the GENCODE loader synthesises ExonBoundary records, the
//       PWM scores annotate each junction so the k-mer index (Block 3)
//       knows which junction-spanning k-mers can be trusted.
//   (b) The Multi-Signal Fusion engine (Block 4.5) feeds L_junction
//       into the Wave-Collapse likelihood; canonical junctions score
//       higher, non-canonical scores lower but are NEVER zero (floor
//       preserved by the floor-clamp in fusion::LikelihoodFactors).
//
// All scores ∈ [0,1]. PWMs are normalised log-likelihood ratios scaled
// to that range; the underlying frequencies trace back to the references
// cited above.
//
// Design choices:
//
//   * U2 (major) vs U12 (minor) discrimination is explicit. The
//     proportions in mammals are ~99.6 % U2 / ~0.4 % U12, but the U12
//     class matters disproportionately at certain IGH loci.
//
//   * Back-splice (circRNA) detection is a structural check, not a PWM:
//     the donor must come AFTER the acceptor in linear genome coords.
//     We still verify the 2-bp motifs (the back-splice donor/acceptor
//     are still GT/AG canonical in most cases).
//
//   * Branch-point detection runs a sliding 7-mer scan over the last
//     50 bp of the intron; the strongest YNYURAC match wins, returned
//     as the offset from the acceptor AG (negative number).
//
// Reference frequencies (donor U2, position relative to donor cut):
//
//     pos:    -3  -2  -1  | +1 +2 +3 +4 +5 +6
//     base:    A   G   G  |  G  T  R  A  G  T   (R = A/G)
//
//     A:     0.34 0.04 0.04 0.00 0.00 0.45 0.74 0.06 0.00
//     C:     0.36 0.09 0.07 0.00 0.00 0.04 0.05 0.06 0.00
//     G:     0.18 0.79 0.81 1.00 0.00 0.50 0.07 0.81 0.00
//     T:     0.12 0.08 0.08 0.00 1.00 0.00 0.14 0.07 1.00

#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string_view>

namespace llmap::annot {

/// Result of scoring a candidate junction.
struct SpliceScoreResult {
    /// PWM score at the donor in [0,1]. 1.0 = perfect canonical.
    float donor_score{0.0f};
    /// PWM score at the acceptor in [0,1].
    float acceptor_score{0.0f};
    /// Spliceosome class:
    ///   0 = U2 major   (GT-AG / GC-AG)
    ///   1 = U12 minor  (AT-AC)
    ///   2 = non-canonical
    ///   3 = back-splice (still GT-AG in motif, but acceptor < donor)
    std::uint8_t spliceosome_class{2};
    /// Branch-point adenosine offset from the acceptor AG, or -1.
    std::int32_t branch_point_offset{-1};
    /// Length of the polypyrimidine tract before the acceptor.
    std::uint8_t polypyrimidine_len{0};
};

class SpliceSiteDb {
public:
    /// Load the built-in PWMs. Always succeeds.
    void LoadDefaults();

    /// Optional refinement from a hexamer override file
    /// (one line: "<6-mer>\t<score>"). Returns false on parse failure.
    bool LoadHexamerOverrides(const std::filesystem::path& path);

    /// Score a candidate junction.
    ///
    ///   donor_2bp   — exactly 2 bp on the INTRON side of the donor cut
    ///                  (e.g. "GT" for U2 canonical, "GC" for U2 minor
    ///                  canonical, "AT" for U12)
    ///   acceptor_2bp — exactly 2 bp on the INTRON side of the acceptor
    ///                  cut (e.g. "AG" for U2, "AC" for U12)
    ///   intron_3p_region — last up-to-50 bp of the intron before the
    ///                       acceptor — used for branch-point + Py-tract
    ///                       detection. May be shorter; the scorer
    ///                       degrades gracefully.
    ///   intron_5p_region — first up-to-8 bp of the intron after the
    ///                       donor — used for refining donor score.
    SpliceScoreResult ScoreJunction(std::string_view donor_2bp,
                                     std::string_view acceptor_2bp,
                                     std::string_view intron_3p_region,
                                     std::string_view intron_5p_region) const;

    /// Convenience predicate — true iff donor/acceptor positions are
    /// consistent with a back-splice (acceptor < donor in linear coords)
    /// AND the 2-bp motifs are still GT/AG-class.
    bool IsBackSpliceConsistent(std::uint64_t donor_pos,
                                 std::uint64_t acceptor_pos,
                                 std::string_view donor_motif,
                                 std::string_view acceptor_motif) const;

    /// Mean-content polypyrimidine score of a region. Public so callers
    /// can re-use it for QC outside the splice-scoring path.
    /// Returns Y-content fraction ∈ [0,1].
    static float PolypyrimidineScore(std::string_view region) noexcept;

private:
    /// 4×9 donor PWM (U2 major), index order [A,C,G,T] × position.
    std::array<std::array<float, 9>, 4> donor_u2_pwm_{};
    /// 4×6 acceptor PWM (U2 major). Window covers the last 5 bp of
    /// the intron through the +1 first exon base.
    std::array<std::array<float, 6>, 4> acceptor_u2_pwm_{};
    /// 4×8 donor PWM (U12 minor).
    std::array<std::array<float, 8>, 4> donor_u12_pwm_{};
    /// 4×8 acceptor PWM (U12 minor).
    std::array<std::array<float, 8>, 4> acceptor_u12_pwm_{};
    /// 4×7 branch-point PWM (mammalian YNYURAC consensus).
    std::array<std::array<float, 7>, 4> branch_point_pwm_{};

    bool defaults_loaded_{false};
};

}  // namespace llmap::annot
