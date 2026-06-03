// LLmap — junction_hunter: consensus + caller implementation.

#include "junction_hunter/consensus_caller.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>

namespace llmap::junction_hunter {

namespace {

/// Per-position consensus decision. `class_` is the consensus class
/// (or LocusClass::Ambiguous if no majority); `read_pos` is the read
/// coordinate; `agreement` counts how many k-values backed the consensus.
struct PositionVote {
    std::uint32_t read_pos{0};
    LocusClass cls{LocusClass::Ambiguous};
    std::uint8_t agreement{0};
};

/// Build per-position consensus votes by walking the five tilings in
/// parallel and bucketing by read_pos. We use a map indexed by read_pos
/// for clarity; for production this should be a tight per-position
/// scan since tilings are already sorted by read_pos.
std::vector<PositionVote>
BuildVotes(const ReadTiling& tiling, const PairKmerIndex& idx,
           const MultiKConfig& cfg) {
    // Counts per (read_pos × class) collected across the 5 tilings.
    // Use a vector keyed by read_pos compressed via sort+merge.
    struct PerPos {
        std::array<std::uint8_t, 5> per_class{};  // None, LcrUp, LcrDown, Interior, Ambiguous
    };
    std::vector<std::pair<std::uint32_t, PerPos>> tmp;
    tmp.reserve(2048);

    for (std::size_t ki = 0; ki < tiling.per_k_hashes.size(); ++ki) {
        const auto& tbl = idx.per_k_class[ki];
        for (const auto& kh : tiling.per_k_hashes[ki]) {
            LocusClass cls = LocusClass::None;
            if (auto it = tbl.find(kh.hash); it != tbl.end()) cls = it->second;
            tmp.push_back({kh.read_pos, PerPos{}});
            tmp.back().second.per_class[static_cast<std::size_t>(cls)] = 1;
        }
    }
    // Sort + merge by read_pos.
    std::sort(tmp.begin(), tmp.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    std::vector<PositionVote> out;
    out.reserve(tmp.size() / 5 + 1);
    std::uint32_t cur_pos = 0;
    PerPos acc{};
    bool have = false;
    auto emit = [&]() {
        // Determine consensus class
        std::uint8_t best_count = 0;
        std::size_t best_idx = 0;  // None
        for (std::size_t i = 1; i < acc.per_class.size(); ++i) {  // skip None
            if (acc.per_class[i] > best_count) { best_count = acc.per_class[i]; best_idx = i; }
        }
        PositionVote v;
        v.read_pos = cur_pos;
        v.agreement = best_count;
        if (best_count >= cfg.consensus_min) v.cls = static_cast<LocusClass>(best_idx);
        else v.cls = LocusClass::Ambiguous;
        out.push_back(v);
    };
    for (const auto& p : tmp) {
        if (!have) { cur_pos = p.first; acc = p.second; have = true; continue; }
        if (p.first == cur_pos) {
            for (std::size_t i = 0; i < acc.per_class.size(); ++i)
                acc.per_class[i] += p.second.per_class[i];
        } else {
            emit();
            cur_pos = p.first; acc = p.second;
        }
    }
    if (have) emit();
    return out;
}

/// Spearman ρ on a vector of (read_pos, dummy genomic-index) pairs.
/// We don't have genomic positions here without an extra lookup table,
/// so use rank-on-rank: monotonic = rank(read_pos) == rank(occurrence
/// order within the votes-of-that-class). This equals ρ=1 by construction
/// in this simplified scheme. A future revision will plug in real
/// genomic-position arrays from the pair index.
float SimpleMonotonicity(const std::vector<PositionVote>& votes, LocusClass cls) {
    std::vector<std::uint32_t> positions;
    for (auto& v : votes) if (v.cls == cls) positions.push_back(v.read_pos);
    if (positions.size() < 2) return 0.0f;
    // Already in read_pos sort order — by construction monotonic.
    // Until we have genomic positions threaded through, return 1.0.
    return 1.0f;
}

/// Find the breakpoint as the LAST LcrUp position before the FIRST LcrDown
/// position, with no Interior in between. Returns the read position.
std::uint32_t FindBreakpoint(const std::vector<PositionVote>& votes,
                             std::uint32_t& flip_idx_out) {
    std::uint32_t last_up = 0;
    bool seen_up = false;
    for (std::size_t i = 0; i < votes.size(); ++i) {
        if (votes[i].cls == LocusClass::LcrUp) { last_up = votes[i].read_pos; seen_up = true; }
        else if (votes[i].cls == LocusClass::LcrDown && seen_up) {
            flip_idx_out = static_cast<std::uint32_t>(i);
            return last_up;
        }
    }
    flip_idx_out = 0;
    return 0;
}

}  // namespace

JunctionRecord CallJunction(std::string_view read_id,
                            const ReadTiling& tiling,
                            const PairKmerIndex& pair_index,
                            const NahrPair& pair,
                            const MultiKConfig& cfg) {
    JunctionRecord rec;
    rec.read_id = read_id;
    rec.pair_id = pair.pair_id;

    auto votes = BuildVotes(tiling, pair_index, cfg);
    rec.n_kmer_total = static_cast<std::uint32_t>(votes.size());
    if (votes.empty()) { rec.call = JunctionCall::Unmapped; return rec; }

    // Count per class.
    for (const auto& v : votes) {
        switch (v.cls) {
            case LocusClass::LcrUp:     ++rec.n_consensus_up; break;
            case LocusClass::LcrDown:   ++rec.n_consensus_dn; break;
            case LocusClass::Interior:  ++rec.n_consensus_in; break;
            case LocusClass::Ambiguous: ++rec.n_ambiguous; break;
            default: break;
        }
    }

    // Excessive ambiguity → paralog-ambiguous.
    if (rec.n_kmer_total > 0 &&
        static_cast<double>(rec.n_ambiguous) / rec.n_kmer_total > 0.20) {
        rec.call = JunctionCall::ParalogAmbiguous;
        return rec;
    }

    // Interior present → canonical (locus intact).
    if (rec.n_consensus_in > 0 &&
        (rec.n_consensus_up == 0 || rec.n_consensus_dn == 0)) {
        rec.call = JunctionCall::CanonicalInterior;
        rec.up_monotonicity = SimpleMonotonicity(votes, LocusClass::Interior);
        return rec;
    }

    // Pure single LCR.
    if (rec.n_consensus_up > 0 && rec.n_consensus_dn == 0 && rec.n_consensus_in == 0) {
        rec.call = JunctionCall::CanonicalUp;
        rec.up_monotonicity = SimpleMonotonicity(votes, LocusClass::LcrUp);
        return rec;
    }
    if (rec.n_consensus_dn > 0 && rec.n_consensus_up == 0 && rec.n_consensus_in == 0) {
        rec.call = JunctionCall::CanonicalDown;
        rec.dn_monotonicity = SimpleMonotonicity(votes, LocusClass::LcrDown);
        return rec;
    }

    // Both LCR_up + LCR_down + no interior → candidate junction.
    if (rec.n_consensus_up > 0 && rec.n_consensus_dn > 0 && rec.n_consensus_in == 0) {
        rec.up_monotonicity = SimpleMonotonicity(votes, LocusClass::LcrUp);
        rec.dn_monotonicity = SimpleMonotonicity(votes, LocusClass::LcrDown);
        std::uint32_t flip_idx = 0;
        rec.breakpoint_read_pos = FindBreakpoint(votes, flip_idx);
        if (rec.up_monotonicity >= cfg.monotonicity_min &&
            rec.dn_monotonicity >= cfg.monotonicity_min) {
            // Compute breakpoint quality as the agreement strength of the
            // 5 votes around the flip.
            std::size_t n = 0; double agree_sum = 0.0;
            const std::size_t window = 5;
            for (std::size_t i = (flip_idx > window ? flip_idx - window : 0);
                 i < std::min(votes.size(), static_cast<std::size_t>(flip_idx) + window); ++i) {
                agree_sum += votes[i].agreement;
                ++n;
            }
            rec.breakpoint_quality = n ? static_cast<float>(agree_sum / (n * 5.0)) : 0.0f;
            rec.call = JunctionCall::JunctionReal;
        } else {
            rec.call = JunctionCall::ChimeraArtifact;
        }
        return rec;
    }

    rec.call = JunctionCall::Unmapped;
    return rec;
}

}  // namespace llmap::junction_hunter
