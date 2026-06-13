// LLmap — Cluster → mini-contig implementation (Line A: probe + lossless layout).

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

std::string RevComp(std::string_view s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'a': c = 't'; break;
            case 'C': c = 'G'; break;
            case 'c': c = 'g'; break;
            case 'G': c = 'C'; break;
            case 'g': c = 'c'; break;
            case 'T': c = 'A'; break;
            case 't': c = 'a'; break;
            default: break;
        }
    }
    return r;
}

// Modal offset (backbone_pos - read_pos) over shared k-mers, and the vote count.
struct Anchor {
    std::int64_t  offset{0};
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

// Representative base for the TRANSIENT probe only (region-finding). This is the
// one place smoothing is allowed — it never touches the output layout.
char ProbeBase(std::uint32_t a, std::uint32_t c, std::uint32_t g,
               std::uint32_t t) {
    std::uint32_t m = a;
    char b = 'A';
    if (c > m) { m = c; b = 'C'; }
    if (g > m) { m = g; b = 'G'; }
    if (t > m) { m = t; b = 'T'; }
    return m == 0 ? 'N' : b;
}

int BaseIdx(char ch) {
    switch (ch) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

}  // namespace

ConsensusContig AssembleConsensus(std::span<const std::string> reads,
                                  const ConsensusConfig& cfg,
                                  std::span<const std::uint8_t> reverse) {
    ConsensusContig out;
    if (reads.empty()) return out;

    auto is_rev = [&](std::size_t i) -> bool {
        return i < reverse.size() && reverse[i] != 0;
    };
    // Orient every read into the contig frame (revcomp the reverse-strand ones).
    // The original strand is recorded per member; nothing is lost.
    std::vector<std::string> oriented(reads.size());
    for (std::size_t i = 0; i < reads.size(); ++i) {
        oriented[i] = is_rev(i) ? RevComp(reads[i]) : reads[i];
    }

    if (reads.size() == 1) {
        out.probe = oriented[0];
        out.length = out.probe.size();
        out.members.push_back({0, 0, true, is_rev(0)});
        out.depth.assign(out.length, 1);
        out.n_anchored = 1;
        return out;
    }

    // Backbone = longest oriented read.
    std::size_t bb = 0;
    for (std::size_t i = 1; i < oriented.size(); ++i) {
        if (oriented[i].size() > oriented[bb].size()) bb = i;
    }
    const std::string& backbone = oriented[bb];
    const std::uint32_t k = cfg.k == 0 ? 15 : cfg.k;

    std::unordered_map<std::string, std::uint32_t> bkmer;
    if (backbone.size() >= k) {
        for (std::size_t i = 0; i + k <= backbone.size(); ++i) {
            bkmer.emplace(backbone.substr(i, k), static_cast<std::uint32_t>(i));
        }
    }

    // Place each read in the backbone frame (backbone itself at offset 0).
    std::vector<std::int64_t> offset(oriented.size(), 0);
    std::vector<bool> anchored(oriented.size(), false);
    std::int64_t lo = 0, hi = static_cast<std::int64_t>(backbone.size());
    for (std::size_t i = 0; i < oriented.size(); ++i) {
        if (i == bb) { anchored[i] = true; continue; }
        const Anchor a = AnchorRead(oriented[i], bkmer, k);
        if (a.votes >= cfg.min_shared_kmers) {
            offset[i] = a.offset;
            anchored[i] = true;
            lo = std::min(lo, a.offset);
            hi = std::max(hi, a.offset + static_cast<std::int64_t>(oriented[i].size()));
        }
    }

    // Column tallies over anchored reads → depth (lossless) + probe (transient).
    const std::size_t len = static_cast<std::size_t>(hi - lo);
    std::vector<std::array<std::uint32_t, 4>> col(len, {0, 0, 0, 0});
    out.depth.assign(len, 0);
    for (std::size_t i = 0; i < oriented.size(); ++i) {
        if (!anchored[i]) continue;
        const std::int64_t start = offset[i] - lo;
        for (std::size_t j = 0; j < oriented[i].size(); ++j) {
            const std::size_t p = static_cast<std::size_t>(start) + j;
            ++out.depth[p];                       // coverage, always counted
            const int bi = BaseIdx(oriented[i][j]);
            if (bi >= 0) ++col[p][static_cast<std::size_t>(bi)];
        }
    }

    // Probe = representative base per column. TRANSIENT — region-finding only.
    out.probe.resize(len);
    for (std::size_t p = 0; p < len; ++p) {
        out.probe[p] = ProbeBase(col[p][0], col[p][1], col[p][2], col[p][3]);
    }
    out.length = len;

    // Lossless layout: every read, its contig offset, its original strand.
    for (std::size_t i = 0; i < oriented.size(); ++i) {
        out.members.push_back(
            {i, anchored[i] ? offset[i] - lo : 0, anchored[i], is_rev(i)});
        if (anchored[i]) ++out.n_anchored;
    }
    return out;
}

}  // namespace llmap::mapper
