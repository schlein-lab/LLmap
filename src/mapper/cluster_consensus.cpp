// LLmap — Cluster → consensus mini-contig implementation (V1, forward strand).

#include "mapper/cluster_consensus.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace llmap::mapper {

namespace {

// Modal offset (backbone_pos - read_pos) over shared k-mers, and how many
// k-mers voted for it. offset places the read in the backbone coordinate frame.
struct Anchor {
    std::int64_t offset{0};
    std::uint32_t votes{0};
};

Anchor AnchorRead(std::string_view read,
                  const std::unordered_map<std::string, std::uint32_t>& bkmer,
                  std::uint32_t k) {
    if (read.size() < k) return {};
    std::map<std::int64_t, std::uint32_t> tally;
    for (std::size_t i = 0; i + k <= read.size(); ++i) {
        const auto it = bkmer.find(std::string(read.substr(i, k)));
        if (it == bkmer.end()) continue;
        const std::int64_t off =
            static_cast<std::int64_t>(it->second) - static_cast<std::int64_t>(i);
        ++tally[off];
    }
    Anchor best;
    for (const auto& [off, v] : tally) {
        if (v > best.votes) {
            best.votes = v;
            best.offset = off;
        }
    }
    return best;
}

char MajorityBase(std::uint32_t a, std::uint32_t c, std::uint32_t g,
                  std::uint32_t t) {
    std::uint32_t m = a;
    char b = 'A';
    if (c > m) { m = c; b = 'C'; }
    if (g > m) { m = g; b = 'G'; }
    if (t > m) { m = t; b = 'T'; }
    return m == 0 ? 'N' : b;
}

}  // namespace

ConsensusContig AssembleConsensus(std::span<const std::string> reads,
                                  const ConsensusConfig& cfg) {
    ConsensusContig out;
    if (reads.empty()) return out;
    if (reads.size() == 1) {
        out.sequence = reads[0];
        out.members.push_back({0, 0, true});
        out.n_anchored = 1;
        return out;
    }

    // Backbone = longest read.
    std::size_t bb = 0;
    for (std::size_t i = 1; i < reads.size(); ++i) {
        if (reads[i].size() > reads[bb].size()) bb = i;
    }
    const std::string& backbone = reads[bb];
    const std::uint32_t k = cfg.k == 0 ? 15 : cfg.k;

    std::unordered_map<std::string, std::uint32_t> bkmer;
    if (backbone.size() >= k) {
        for (std::size_t i = 0; i + k <= backbone.size(); ++i) {
            bkmer.emplace(backbone.substr(i, k), static_cast<std::uint32_t>(i));
        }
    }

    // Place each read in the backbone frame (backbone itself at offset 0).
    std::vector<std::int64_t> offset(reads.size(), 0);
    std::vector<bool> anchored(reads.size(), false);
    std::int64_t lo = 0, hi = static_cast<std::int64_t>(backbone.size());
    for (std::size_t i = 0; i < reads.size(); ++i) {
        if (i == bb) {
            anchored[i] = true;
            continue;
        }
        const Anchor a = AnchorRead(reads[i], bkmer, k);
        if (a.votes >= cfg.min_shared_kmers) {
            offset[i] = a.offset;
            anchored[i] = true;
            lo = std::min(lo, a.offset);
            hi = std::max(hi, a.offset + static_cast<std::int64_t>(reads[i].size()));
        }
    }

    // Shift so the contig starts at 0, then column-majority over anchored reads.
    const std::size_t len = static_cast<std::size_t>(hi - lo);
    std::vector<std::array<std::uint32_t, 4>> col(len, {0, 0, 0, 0});
    auto base_idx = [](char ch) -> int {
        switch (ch) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };
    for (std::size_t i = 0; i < reads.size(); ++i) {
        if (!anchored[i]) continue;
        const std::int64_t start = offset[i] - lo;
        for (std::size_t j = 0; j < reads[i].size(); ++j) {
            const int bi = base_idx(reads[i][j]);
            if (bi >= 0) ++col[static_cast<std::size_t>(start) + j][static_cast<std::size_t>(bi)];
        }
    }

    out.sequence.resize(len);
    for (std::size_t p = 0; p < len; ++p) {
        out.sequence[p] = MajorityBase(col[p][0], col[p][1], col[p][2], col[p][3]);
    }
    for (std::size_t i = 0; i < reads.size(); ++i) {
        out.members.push_back(
            {i, anchored[i] ? offset[i] - lo : 0, anchored[i]});
        if (anchored[i]) ++out.n_anchored;
    }
    return out;
}

}  // namespace llmap::mapper
