// LLmap — short-read pipeline implementation.
//
// Three pure helpers:
//   BuildEquivalenceClass — Salmon/Kallisto-style soft assignment
//   ImpliedInsertSize     — paired-end fragment-length estimate
//   BuildOutcome          — compiles the above into a ShortReadOutcome
//
// No I/O, no global state, no threads. The actual FASTQ reading +
// alignment happens upstream; we only finalise the per-record outcome
// for the lossless pipeline.

#include "io/short_read_handler.h"

#include <algorithm>
#include <numeric>

namespace llmap::io {

// ---------------------------------------------------------------------------
// BuildEquivalenceClass
// ---------------------------------------------------------------------------

std::vector<EquivalenceClassMember>
BuildEquivalenceClass(
    const std::vector<std::pair<std::string, float>>& hits,
    const ShortReadConfig& cfg) {

    if (hits.empty()) return {};

    // Sum scores; normalise to weights.
    float total = 0.0f;
    for (const auto& [_, s] : hits) total += std::max(0.0f, s);
    if (total <= 0.0f) return {};

    std::vector<EquivalenceClassMember> members;
    members.reserve(hits.size());
    for (const auto& [tx, s] : hits) {
        if (s <= 0.0f) continue;
        const float w = s / total;
        if (w < cfg.min_member_weight) continue;
        members.push_back({tx, w});
    }

    // Re-normalise after dropping low-weight members so weights still
    // sum to 1.0.
    if (!members.empty()) {
        float kept = 0.0f;
        for (const auto& m : members) kept += m.weight;
        if (kept > 0.0f) {
            for (auto& m : members) m.weight /= kept;
        }
    }
    // Sort descending by weight so consumers can take top-K easily.
    std::sort(members.begin(), members.end(),
              [](const auto& a, const auto& b) {
                  return a.weight > b.weight;
              });
    return members;
}

// ---------------------------------------------------------------------------
// ImpliedInsertSize
// ---------------------------------------------------------------------------

std::uint32_t ImpliedInsertSize(std::uint32_t pos1, std::uint32_t pos2) {
    if (pos1 == 0 || pos2 == 0) return 0;
    return (pos1 > pos2) ? (pos1 - pos2) : (pos2 - pos1);
}

// ---------------------------------------------------------------------------
// BuildOutcome
// ---------------------------------------------------------------------------

ShortReadOutcome BuildOutcome(
    const FastqRecordPair& rec,
    const std::vector<std::pair<std::string, float>>& hits_r1,
    const std::vector<std::pair<std::string, float>>& hits_r2,
    std::optional<std::uint32_t> implied_insert,
    bool junction_seen,
    const ShortReadConfig& cfg) {

    ShortReadOutcome out;
    out.read_id = rec.read_id;

    // Merge hit lists by transcript_id, summing scores. The merged
    // map gives each transcript its combined evidence across mates.
    std::vector<std::pair<std::string, float>> merged;
    auto add_or_merge = [&](const std::string& tx, float s) {
        for (auto& [k, v] : merged) {
            if (k == tx) { v += s; return; }
        }
        merged.emplace_back(tx, s);
    };
    for (const auto& [tx, s] : hits_r1) add_or_merge(tx, s);
    for (const auto& [tx, s] : hits_r2) add_or_merge(tx, s);

    if (cfg.emit_equivalence_classes) {
        out.members = BuildEquivalenceClass(merged, cfg);
    } else {
        // Pick the single best.
        if (!merged.empty()) {
            auto best = std::max_element(
                merged.begin(), merged.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            out.members.push_back({best->first, 1.0f});
        }
    }

    // Insert size + junction inference.
    out.implied_insert_size = implied_insert.value_or(0);
    const bool insert_too_long =
        cfg.expected_insert_size > 0
        && out.implied_insert_size
            > static_cast<std::uint32_t>(
                  cfg.expected_insert_size * 3 / 2);  // 1.5× rule
    out.junction_inferred = junction_seen || insert_too_long;

    return out;
}

}  // namespace llmap::io
