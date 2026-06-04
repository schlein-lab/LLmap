// LLmap — junction_hunter: consensus + caller implementation.
//
// Per read position we collapse the multi-k votes into one consensus
// class, carrying TWO extra signals the earlier prototype discarded:
//   - the in-region genomic offset of the winning k-mer (for the real
//     outside-in monotonicity test), and
//   - the longest k that backed the winning class (for PSV-grade gating).
//
// A both-ends read is only JunctionReal when BOTH copies are traversed
// monotonically (|Spearman ρ| ≥ monotonicity_min on the in-region
// offsets) AND both copies carry ≥ min_psv_switches PSV-grade (long-k)
// consensus positions. Monotonic-but-thin → ParalogAmbiguous; geometry
// broken → ChimeraArtifact. This is the discrimination minimap2 cannot
// do: a chimeric alignment artefact scrambles the offset order and is
// rejected here even though its k-mer membership looks junction-like.

#include "junction_hunter/consensus_caller.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace llmap::junction_hunter {

namespace {

/// Per-position consensus decision.
struct PositionVote {
    std::uint32_t read_pos{0};
    LocusClass cls{LocusClass::Ambiguous};
    std::uint8_t agreement{0};   ///< how many k backed the winning class
    std::uint8_t best_k{0};      ///< longest k backing the winning class
    std::uint32_t offset{0};     ///< in-region offset of the longest-k winner
};

/// One raw (read_pos, class, k, offset) emission before grouping.
struct RawVote {
    std::uint32_t read_pos{0};
    LocusClass cls{LocusClass::None};
    std::uint8_t k{0};
    std::uint32_t offset{0};
};

/// Walk all tilings, resolve every k-mer hit to its (class, offset) via
/// the pair index, then group by read_pos into per-position consensus.
std::vector<PositionVote>
BuildVotes(const ReadTiling& tiling, const PairKmerIndex& idx,
           const MultiKConfig& cfg) {
    std::vector<RawVote> raw;
    raw.reserve(2048);
    for (std::size_t ki = 0; ki < tiling.per_k_hashes.size(); ++ki) {
        const auto& tbl = idx.per_k_class[ki];
        const std::uint8_t k = ki < tiling.k_values.size() ? tiling.k_values[ki] : 0;
        for (const auto& kh : tiling.per_k_hashes[ki]) {
            RawVote rv;
            rv.read_pos = kh.read_pos;
            rv.k = k;
            if (auto it = tbl.find(kh.hash); it != tbl.end()) {
                rv.cls = it->second.cls;
                rv.offset = it->second.offset;
            } else {
                rv.cls = LocusClass::None;
            }
            raw.push_back(rv);
        }
    }
    std::stable_sort(raw.begin(), raw.end(),
        [](const RawVote& a, const RawVote& b){ return a.read_pos < b.read_pos; });

    std::vector<PositionVote> out;
    out.reserve(raw.size() / 5 + 1);

    std::size_t i = 0;
    while (i < raw.size()) {
        std::size_t j = i;
        const std::uint32_t pos = raw[i].read_pos;
        std::array<std::uint8_t, 5> per_class{};   // None,LcrUp,LcrDown,Interior,Ambiguous
        while (j < raw.size() && raw[j].read_pos == pos) {
            per_class[static_cast<std::size_t>(raw[j].cls)] += 1;
            ++j;
        }
        // Winning class: max count over non-None classes.
        std::uint8_t best_count = 0;
        std::size_t best_idx = 0;
        for (std::size_t c = 1; c < per_class.size(); ++c)
            if (per_class[c] > best_count) { best_count = per_class[c]; best_idx = c; }

        // Index hits = votes that landed in any real region (not a miss).
        const std::uint32_t index_hits =
            per_class[1] + per_class[2] + per_class[3] + per_class[4];

        PositionVote v;
        v.read_pos = pos;
        v.agreement = best_count;
        if (best_count >= cfg.consensus_min) {
            v.cls = static_cast<LocusClass>(best_idx);
            // Among raw votes of the winning class at this position, take
            // the longest k and its offset — the long k carries the PSV.
            for (std::size_t r = i; r < j; ++r) {
                if (raw[r].cls == v.cls && raw[r].k >= v.best_k) {
                    v.best_k = raw[r].k;
                    v.offset = raw[r].offset;
                }
            }
        } else if (index_hits > 0) {
            // Real index hits but no majority class → genuinely ambiguous
            // (k-mer shared between copies at this position).
            v.cls = LocusClass::Ambiguous;
        } else {
            // No region hit at all — a true miss, not paralog-ambiguity.
            // Kept in the stream (counts toward n_kmer_total) but does not
            // contribute signal, so an unrelated read stays Unmapped.
            v.cls = LocusClass::None;
        }
        out.push_back(v);
        i = j;
    }
    return out;
}

/// |Spearman ρ| between read positions (already sorted, distinct → ranks
/// 0..n-1) and the in-region offsets of one class. Strand-agnostic: an
/// antisense walk gives ρ≈-1, still |ρ|≈1. Ties in offset get arbitrary
/// adjacent ranks (rare for unique long k); n<2 → 0 (cannot assess).
float SpearmanAbs(const std::vector<std::pair<std::uint32_t, std::uint32_t>>& pts) {
    const std::size_t n = pts.size();
    if (n < 2) return 0.0f;
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [&](std::size_t a, std::size_t b){ return pts[a].second < pts[b].second; });
    std::vector<double> off_rank(n);
    for (std::size_t r = 0; r < n; ++r) off_rank[order[r]] = static_cast<double>(r);
    double d2 = 0.0;
    for (std::size_t r = 0; r < n; ++r) {
        const double d = static_cast<double>(r) - off_rank[r];   // read-rank = r
        d2 += d * d;
    }
    const double denom = static_cast<double>(n) * (static_cast<double>(n) * n - 1.0);
    const double rho = 1.0 - 6.0 * d2 / denom;
    return static_cast<float>(std::fabs(rho));
}

/// Collect (read_pos, offset) pairs of one class, in read order.
std::vector<std::pair<std::uint32_t, std::uint32_t>>
PointsOf(const std::vector<PositionVote>& votes, LocusClass cls) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pts;
    for (const auto& v : votes)
        if (v.cls == cls) pts.emplace_back(v.read_pos, v.offset);
    return pts;
}

/// Last LcrUp position before the first following LcrDown (no interior
/// between). Emits the flanking offsets and the flip index.
std::uint32_t FindBreakpoint(const std::vector<PositionVote>& votes,
                             std::uint32_t& flip_idx_out,
                             std::uint32_t& up_off_out,
                             std::uint32_t& dn_off_out) {
    std::uint32_t last_up = 0, last_up_off = 0;
    bool seen_up = false;
    for (std::size_t i = 0; i < votes.size(); ++i) {
        if (votes[i].cls == LocusClass::LcrUp) {
            last_up = votes[i].read_pos; last_up_off = votes[i].offset; seen_up = true;
        } else if (votes[i].cls == LocusClass::LcrDown && seen_up) {
            flip_idx_out = static_cast<std::uint32_t>(i);
            up_off_out = last_up_off;
            dn_off_out = votes[i].offset;
            return last_up;
        }
    }
    flip_idx_out = 0; up_off_out = 0; dn_off_out = 0;
    return 0;
}

/// PSV-grade consensus positions of one class (winning class backed by a
/// k ≥ psv_k_min). These carry paralog-discriminating variants.
std::uint32_t CountPsv(const std::vector<PositionVote>& votes,
                       LocusClass cls, std::uint8_t psv_k_min) {
    std::uint32_t n = 0;
    for (const auto& v : votes)
        if (v.cls == cls && v.best_k >= psv_k_min) ++n;
    return n;
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

    for (const auto& v : votes) {
        switch (v.cls) {
            case LocusClass::LcrUp:     ++rec.n_consensus_up; break;
            case LocusClass::LcrDown:   ++rec.n_consensus_dn; break;
            case LocusClass::Interior:  ++rec.n_consensus_in; break;
            case LocusClass::Ambiguous: ++rec.n_ambiguous;    break;
            default: break;
        }
    }
    rec.n_psv_up = CountPsv(votes, LocusClass::LcrUp,   cfg.psv_k_min);
    rec.n_psv_dn = CountPsv(votes, LocusClass::LcrDown, cfg.psv_k_min);

    // Ambiguous positions are NOISE, not a verdict: at LCR identity <99 %
    // a large fraction of short-k matches are identical between the copies.
    // Disambiguation comes from the long k that carry the PSVs. We only
    // fall to ParalogAmbiguous when there is NO unambiguous signal at all.
    const std::uint32_t n_signal =
        rec.n_consensus_up + rec.n_consensus_dn + rec.n_consensus_in;
    if (n_signal == 0 && rec.n_ambiguous > 0) {
        rec.call = JunctionCall::ParalogAmbiguous;
        return rec;
    }

    // Interior present without a clean both-ends split → locus intact.
    if (rec.n_consensus_in > 0 &&
        (rec.n_consensus_up == 0 || rec.n_consensus_dn == 0)) {
        rec.call = JunctionCall::CanonicalInterior;
        rec.up_monotonicity = SpearmanAbs(PointsOf(votes, LocusClass::Interior));
        return rec;
    }

    // Single LCR copy only.
    if (rec.n_consensus_up > 0 && rec.n_consensus_dn == 0 && rec.n_consensus_in == 0) {
        rec.call = JunctionCall::CanonicalUp;
        rec.up_monotonicity = SpearmanAbs(PointsOf(votes, LocusClass::LcrUp));
        return rec;
    }
    if (rec.n_consensus_dn > 0 && rec.n_consensus_up == 0 && rec.n_consensus_in == 0) {
        rec.call = JunctionCall::CanonicalDown;
        rec.dn_monotonicity = SpearmanAbs(PointsOf(votes, LocusClass::LcrDown));
        return rec;
    }

    // Both copies, no interior → junction candidate. Now the real test.
    if (rec.n_consensus_up > 0 && rec.n_consensus_dn > 0 && rec.n_consensus_in == 0) {
        rec.up_monotonicity = SpearmanAbs(PointsOf(votes, LocusClass::LcrUp));
        rec.dn_monotonicity = SpearmanAbs(PointsOf(votes, LocusClass::LcrDown));

        std::uint32_t flip_idx = 0, up_off = 0, dn_off = 0;
        rec.breakpoint_read_pos = FindBreakpoint(votes, flip_idx, up_off, dn_off);
        rec.breakpoint_genomic_up = pair.lcr_up_start   + up_off;
        rec.breakpoint_genomic_dn = pair.lcr_down_start + dn_off;

        const bool monotonic =
            rec.up_monotonicity >= cfg.monotonicity_min &&
            rec.dn_monotonicity >= cfg.monotonicity_min;
        const bool psv_both =
            rec.n_psv_up >= cfg.min_psv_switches &&
            rec.n_psv_dn >= cfg.min_psv_switches;

        if (!monotonic) {
            // Geometry broken — the membership looks junction-like but the
            // offsets do not walk each copy in order. Classic chimeric
            // alignment artefact (the minimap2 failure mode at LCR<99 %).
            rec.call = JunctionCall::ChimeraArtifact;
        } else if (!psv_both) {
            // Monotonic but one copy lacks PSV-grade (long-k) evidence —
            // cannot confirm both copies were genuinely traversed.
            rec.call = JunctionCall::ParalogAmbiguous;
        } else {
            // Breakpoint quality = mean multi-k agreement in the flip window.
            std::size_t n = 0; double agree_sum = 0.0;
            const std::size_t window = 5;
            const std::size_t lo = flip_idx > window ? flip_idx - window : 0;
            const std::size_t hi =
                std::min(votes.size(), static_cast<std::size_t>(flip_idx) + window);
            for (std::size_t i = lo; i < hi; ++i) { agree_sum += votes[i].agreement; ++n; }
            rec.breakpoint_quality =
                n ? static_cast<float>(agree_sum / (n * 5.0)) : 0.0f;
            rec.call = JunctionCall::JunctionReal;
        }
        return rec;
    }

    rec.call = JunctionCall::Unmapped;
    return rec;
}

}  // namespace llmap::junction_hunter
