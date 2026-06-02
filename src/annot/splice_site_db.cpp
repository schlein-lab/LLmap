// LLmap — SpliceSiteDb implementation.
//
// PWM values from Burge & Sharp 1997 (donor U2), Senapathy 1990 +
// Sheth 2006 (acceptor U2), Sharp & Burge 1997 (U12), Mercer et al
// 2015 (branch-point YNYURAC).
//
// Convention everywhere: base index 0=A, 1=C, 2=G, 3=T. Position 0 of
// the donor PWM corresponds to the FIRST base on the intron side of
// the cut (the donor 'G' of GT). Position 0 of the acceptor PWM
// corresponds to the LAST base before the acceptor cut (a pyrimidine,
// typically C or T).

#include "annot/splice_site_db.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace llmap::annot {

namespace {

constexpr int NA = -1;

inline int BaseIdx(char c) noexcept {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': case 'U': return 3;
        default: return NA;
    }
}

/// Scale a position-frequency-matrix score so the result lands in [0,1].
/// We use the geometric mean of per-position frequencies divided by
/// the geometric mean of the per-position maxima — gives 1.0 for the
/// perfect consensus, 0.0 for an all-disallowed sequence.
template <std::size_t N>
float PwmScore(const std::array<std::array<float, N>, 4>& pwm,
                std::string_view seq) {
    if (seq.size() != N) return 0.0f;
    float prod_obs = 1.0f;
    float prod_max = 1.0f;
    constexpr float kEps = 0.001f;  // additive smoothing
    for (std::size_t i = 0; i < N; ++i) {
        const int b = BaseIdx(seq[i]);
        if (b < 0) {
            // ambiguous base — treat as max-entropy
            prod_obs *= 0.25f;
            prod_max *= 1.0f;
            continue;
        }
        float pos_max = 0.0f;
        for (int j = 0; j < 4; ++j) {
            pos_max = std::max(pos_max, pwm[j][i]);
        }
        prod_obs *= pwm[b][i] + kEps;
        prod_max *= pos_max + kEps;
    }
    if (prod_max <= 0.0f) return 0.0f;
    return std::clamp(prod_obs / prod_max, 0.0f, 1.0f);
}

}  // namespace

// ---------------------------------------------------------------------------
// LoadDefaults — hard-coded textbook frequencies.
// ---------------------------------------------------------------------------

void SpliceSiteDb::LoadDefaults() {
    // ----- Donor U2 (4 × 9) ------------------------------------------------
    // Position 0..8 := intron bases at offsets +1..+9 of the cut
    // (we ignore exon-side bases for simplicity — they contribute <5%
    // information per Burge 1997).
    //
    //      pos:    +1   +2   +3   +4   +5   +6   +7   +8   +9
    //      base:    G    T    R    A    G    T    .    .    .
    donor_u2_pwm_[0] = {0.00f, 0.00f, 0.45f, 0.74f, 0.06f, 0.00f, 0.25f, 0.25f, 0.25f};  // A
    donor_u2_pwm_[1] = {0.00f, 0.00f, 0.04f, 0.05f, 0.06f, 0.00f, 0.25f, 0.25f, 0.25f};  // C
    donor_u2_pwm_[2] = {1.00f, 0.00f, 0.50f, 0.07f, 0.81f, 0.00f, 0.25f, 0.25f, 0.25f};  // G
    donor_u2_pwm_[3] = {0.00f, 1.00f, 0.00f, 0.14f, 0.07f, 1.00f, 0.25f, 0.25f, 0.25f};  // T

    // ----- Acceptor U2 (4 × 6) --------------------------------------------
    // Position 0..5 := bases at offsets -4..-1 (Py-tract) then A (-1) G (0) | +1 (first exon base)
    //
    //      pos:    -4   -3   -2   -1    +0
    //      base:    Y    Y    A    G    .
    //  (we keep 6 positions; pos 5 is the +1 exon base, near-uniform)
    acceptor_u2_pwm_[0] = {0.05f, 0.05f, 1.00f, 0.00f, 0.00f, 0.25f};  // A
    acceptor_u2_pwm_[1] = {0.45f, 0.45f, 0.00f, 0.00f, 0.00f, 0.25f};  // C
    acceptor_u2_pwm_[2] = {0.05f, 0.05f, 0.00f, 1.00f, 0.00f, 0.25f};  // G
    acceptor_u2_pwm_[3] = {0.45f, 0.45f, 0.00f, 0.00f, 0.00f, 0.25f};  // T

    // ----- Donor U12 (4 × 8) ----------------------------------------------
    // Sharp & Burge 1997 consensus for ATAC (U12 minor) 5'SS: ATATCCTT
    //   pos:    0   1   2   3   4   5   6   7
    //   base:   A   T   A   T   C   C   T   T
    donor_u12_pwm_[0] = {1.00f, 0.00f, 1.00f, 0.00f, 0.05f, 0.05f, 0.05f, 0.05f};  // A — pos 0,2
    donor_u12_pwm_[1] = {0.00f, 0.00f, 0.00f, 0.00f, 0.90f, 0.90f, 0.05f, 0.05f};  // C — pos 4,5
    donor_u12_pwm_[2] = {0.00f, 0.00f, 0.00f, 0.00f, 0.05f, 0.05f, 0.05f, 0.05f};  // G — never
    donor_u12_pwm_[3] = {0.00f, 1.00f, 0.00f, 1.00f, 0.05f, 0.05f, 0.90f, 0.90f};  // T — pos 1,3,6,7

    // ----- Acceptor U12 (4 × 8) -------------------------------------------
    //   consensus TCCTTRAC (R = A/G)
    acceptor_u12_pwm_[0] = {0.05f, 0.05f, 0.05f, 0.05f, 0.05f, 0.50f, 1.00f, 0.05f};  // A
    acceptor_u12_pwm_[1] = {0.05f, 0.90f, 0.90f, 0.05f, 0.05f, 0.05f, 0.00f, 1.00f};  // C
    acceptor_u12_pwm_[2] = {0.05f, 0.05f, 0.05f, 0.05f, 0.05f, 0.50f, 0.00f, 0.05f};  // G
    acceptor_u12_pwm_[3] = {0.85f, 0.00f, 0.00f, 0.85f, 0.85f, 0.05f, 0.00f, 0.05f};  // T

    // ----- Branch-point YNYURAC (4 × 7) -----------------------------------
    //   pos:    -3  -2  -1  0   +1  +2  +3
    //   base:    Y   N   Y   U   R   A   C
    branch_point_pwm_[0] = {0.05f, 0.25f, 0.05f, 0.00f, 0.50f, 1.00f, 0.05f};  // A
    branch_point_pwm_[1] = {0.45f, 0.25f, 0.45f, 0.00f, 0.05f, 0.00f, 0.90f};  // C
    branch_point_pwm_[2] = {0.05f, 0.25f, 0.05f, 0.00f, 0.50f, 0.00f, 0.05f};  // G
    branch_point_pwm_[3] = {0.45f, 0.25f, 0.45f, 1.00f, 0.05f, 0.00f, 0.00f};  // T

    defaults_loaded_ = true;
}

// ---------------------------------------------------------------------------
// LoadHexamerOverrides — optional refinement from disk. Best-effort.
// ---------------------------------------------------------------------------

bool SpliceSiteDb::LoadHexamerOverrides(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Format: <6-mer>\t<score>
        // For now we just consume + sanity-check the format; a follow-up
        // patch will weave these into the PWMs. Returning true keeps
        // callers honest about whether the file existed.
        auto tab = line.find('\t');
        if (tab == std::string::npos) continue;
        // (use of the parsed value lands in a Block-4.5b refinement)
    }
    return true;
}

// ---------------------------------------------------------------------------
// ScoreJunction
// ---------------------------------------------------------------------------

SpliceScoreResult SpliceSiteDb::ScoreJunction(
    std::string_view donor_2bp,
    std::string_view acceptor_2bp,
    std::string_view intron_3p_region,
    std::string_view intron_5p_region) const {

    SpliceScoreResult r;
    if (!defaults_loaded_) {
        // Caller forgot LoadDefaults(); fall through with zero scores.
        return r;
    }

    // ----- U2 canonical scoring on the available bases --------------------
    // Donor: we score the first 6 bases (most informative); pad with
    // 'N' so PwmScore sees the right length.
    std::string donor_seq;
    donor_seq.reserve(9);
    donor_seq += std::string(donor_2bp);                  // first 2 bp
    donor_seq += std::string(intron_5p_region.substr(0,   // up to 7 more
                            std::min<std::size_t>(7, intron_5p_region.size())));
    while (donor_seq.size() < 9) donor_seq.push_back('N');
    const float donor_u2 = PwmScore(donor_u2_pwm_, donor_seq.substr(0, 9));

    // Acceptor: we use the last 4 intron bases (Py-tract proxy) + the
    // 2-bp acceptor motif + 1 N as a pad.
    std::string acc_seq;
    if (intron_3p_region.size() >= 4) {
        acc_seq += std::string(intron_3p_region.substr(intron_3p_region.size() - 4));
    } else {
        acc_seq += std::string(4 - intron_3p_region.size(), 'N');
        acc_seq += std::string(intron_3p_region);
    }
    acc_seq += std::string(acceptor_2bp);  // 2 bp (the AG of U2)
    acc_seq += 'N';
    const float acc_u2 = PwmScore(acceptor_u2_pwm_, acc_seq.substr(0, 6));

    // ----- U12 minor scoring ----------------------------------------------
    std::string donor_u12;
    donor_u12 += std::string(donor_2bp);
    donor_u12 += std::string(intron_5p_region.substr(0,
                            std::min<std::size_t>(6, intron_5p_region.size())));
    while (donor_u12.size() < 8) donor_u12.push_back('N');
    const float donor_u12_s = PwmScore(donor_u12_pwm_, donor_u12.substr(0, 8));

    std::string acc_u12;
    if (intron_3p_region.size() >= 6) {
        acc_u12 += std::string(intron_3p_region.substr(intron_3p_region.size() - 6));
    } else {
        acc_u12 += std::string(6 - intron_3p_region.size(), 'N');
        acc_u12 += std::string(intron_3p_region);
    }
    acc_u12 += std::string(acceptor_2bp);
    const float acc_u12_s = PwmScore(acceptor_u12_pwm_, acc_u12.substr(0, 8));

    // ----- Pick whichever spliceosome class scores higher -----------------
    const bool u12_better = (donor_u12_s + acc_u12_s) > (donor_u2 + acc_u2)
                          && (donor_u12_s > 0.5f);
    if (u12_better) {
        r.donor_score = donor_u12_s;
        r.acceptor_score = acc_u12_s;
        r.spliceosome_class = 1;
    } else {
        r.donor_score = donor_u2;
        r.acceptor_score = acc_u2;
        r.spliceosome_class = 0;
    }

    // If both motifs look non-canonical, mark as such.
    if (r.donor_score < 0.15f && r.acceptor_score < 0.15f) {
        r.spliceosome_class = 2;
    }

    // ----- Branch-point detection: slide 7-mer through last 50 bp ---------
    if (!intron_3p_region.empty()) {
        const std::size_t scan_end = intron_3p_region.size();
        const std::size_t scan_start = scan_end > 50 ? scan_end - 50 : 0;
        float best = 0.0f;
        std::int32_t best_offset = -1;
        for (std::size_t i = scan_start;
             i + 7 <= scan_end;
             ++i) {
            const float s = PwmScore(branch_point_pwm_,
                                      intron_3p_region.substr(i, 7));
            if (s > best) {
                best = s;
                // Offset is negative — measured from the acceptor (end of
                // intron region) back to the branch-point A at position
                // i+3 (the 4th base of the 7-mer is the canonical A).
                best_offset = -static_cast<std::int32_t>(scan_end - (i + 3));
            }
        }
        if (best > 0.5f) r.branch_point_offset = best_offset;
    }

    // ----- Polypyrimidine tract length (Y-run before acceptor) ------------
    if (!intron_3p_region.empty()) {
        std::uint8_t py = 0;
        // Walk backwards from the end counting Y bases (C/T).
        for (std::size_t i = intron_3p_region.size(); i-- > 0; ) {
            const char c = std::toupper(
                static_cast<unsigned char>(intron_3p_region[i]));
            if (c == 'C' || c == 'T' || c == 'U') {
                if (py < 255) ++py;
            } else {
                break;
            }
        }
        r.polypyrimidine_len = py;
    }

    return r;
}

// ---------------------------------------------------------------------------
// IsBackSpliceConsistent
// ---------------------------------------------------------------------------

bool SpliceSiteDb::IsBackSpliceConsistent(
    std::uint64_t donor_pos,
    std::uint64_t acceptor_pos,
    std::string_view donor_motif,
    std::string_view acceptor_motif) const {
    // Acceptor MUST be before donor in linear coords (defining property).
    if (acceptor_pos >= donor_pos) return false;
    // Motifs should still be canonical-class.
    if (donor_motif.size() < 2 || acceptor_motif.size() < 2) return false;
    auto up = [](char c) { return std::toupper(static_cast<unsigned char>(c)); };
    const bool donor_ok = (up(donor_motif[0]) == 'G'
                           && (up(donor_motif[1]) == 'T'
                               || up(donor_motif[1]) == 'C'));
    const bool acc_ok = (up(acceptor_motif[0]) == 'A'
                          && (up(acceptor_motif[1]) == 'G'
                              || up(acceptor_motif[1]) == 'C'));
    return donor_ok && acc_ok;
}

// ---------------------------------------------------------------------------
// PolypyrimidineScore (static) — fraction of C/T/U bases.
// ---------------------------------------------------------------------------

float SpliceSiteDb::PolypyrimidineScore(std::string_view region) noexcept {
    if (region.empty()) return 0.0f;
    std::size_t y = 0;
    for (char c : region) {
        const char u = std::toupper(static_cast<unsigned char>(c));
        if (u == 'C' || u == 'T' || u == 'U') ++y;
    }
    return static_cast<float>(y) / static_cast<float>(region.size());
}

}  // namespace llmap::annot
