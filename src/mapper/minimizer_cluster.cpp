// LLmap — Minimizer-overlap CPU read clustering implementation.

#include "mapper/minimizer_cluster.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace llmap::mapper {

namespace {

std::uint8_t Code(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0xFF;
    }
}

// One canonical minimizer occurrence in a read.
struct Mini {
    std::uint64_t canon{0};  // canonical hash (min of fwd / revcomp)
    std::uint32_t pos{0};    // read position
    bool          rc{false}; // canonical came from the reverse-complement strand
};

// Canonical minimizers of a read (window w over k-mers). Dedups consecutive.
std::vector<Mini> Minimizers(std::string_view read, std::uint32_t k,
                             std::uint32_t w) {
    std::vector<Mini> out;
    if (read.size() < k) return out;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    std::uint64_t fwd = 0, rev = 0;
    std::vector<Mini> kmers;
    std::uint32_t valid = 0;
    for (std::uint32_t i = 0; i < read.size(); ++i) {
        const std::uint8_t c = Code(read[i]);
        if (c == 0xFF) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | c) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - c) << (2 * (k - 1)));
        if (++valid >= k) {
            const std::uint32_t pos = i - k + 1;
            const bool rc = rev < fwd;
            kmers.push_back({rc ? rev : fwd, pos, rc});
        }
    }
    // Minimizer per window of w consecutive k-mers.
    std::uint64_t last = ~0ULL;
    for (std::size_t i = 0; i + w <= kmers.size(); ++i) {
        std::size_t best = i;
        for (std::size_t j = i + 1; j < i + w; ++j) {
            if (kmers[j].canon < kmers[best].canon) best = j;
        }
        if (kmers[best].canon != last) {
            out.push_back(kmers[best]);
            last = kmers[best].canon;
        }
    }
    return out;
}

// Signed (parity) union-find: tracks each node's orientation relative to its root.
struct SignedDsu {
    std::vector<std::uint32_t> parent;
    std::vector<std::uint8_t> rel;  // orientation relative to parent (0 same, 1 flip)
    explicit SignedDsu(std::size_t n) : parent(n), rel(n, 0) {
        for (std::uint32_t i = 0; i < n; ++i) parent[i] = i;
    }
    // Returns {root, parity-to-root}.
    std::pair<std::uint32_t, std::uint8_t> find(std::uint32_t x) {
        std::uint8_t p = 0;
        while (parent[x] != x) {
            p ^= rel[x];
            x = parent[x];
        }
        return {x, p};
    }
    void unite(std::uint32_t a, std::uint32_t b, std::uint8_t flip) {
        auto [ra, pa] = find(a);
        auto [rb, pb] = find(b);
        if (ra == rb) return;
        parent[ra] = rb;
        rel[ra] = pa ^ pb ^ flip;  // orientation of ra relative to rb
    }
};

}  // namespace

std::vector<ReadCluster> ClusterByMinimizers(std::span<const std::string> reads,
                                             const MinimizerClusterConfig& cfg) {
    const std::size_t n = reads.size();
    std::vector<ReadCluster> result(n);
    if (n == 0) return result;
    const std::uint32_t k = cfg.k == 0 ? 15 : cfg.k;
    const std::uint32_t w = cfg.w == 0 ? 10 : cfg.w;

    // Minimizers per read + inverted index canon → [(read, rc)].
    std::vector<std::vector<Mini>> mins(n);
    std::unordered_map<std::uint64_t, std::vector<std::pair<std::uint32_t, bool>>> idx;
    for (std::size_t i = 0; i < n; ++i) {
        mins[i] = Minimizers(reads[i], k, w);
        for (const auto& m : mins[i]) {
            idx[m.canon].push_back({static_cast<std::uint32_t>(i), m.rc});
        }
    }

    // Frequency cap (self-frequency = read-set repeat proxy): drop minimizers in
    // too many reads — they over-merge paralog/repeat families.
    const std::size_t cap = std::max<std::size_t>(
        2, static_cast<std::size_t>(cfg.max_freq_frac * static_cast<double>(n)));

    // Pairwise shared count + orientation-flip tally.
    struct PairAcc { std::uint32_t shared{0}; std::uint32_t flips{0}; };
    std::map<std::pair<std::uint32_t, std::uint32_t>, PairAcc> pairs;
    for (const auto& [canon, posting] : idx) {
        if (posting.size() < 2 || posting.size() > cap) continue;  // skip repetitive
        for (std::size_t a = 0; a < posting.size(); ++a) {
            for (std::size_t b = a + 1; b < posting.size(); ++b) {
                std::uint32_t ra = posting[a].first, rb = posting[b].first;
                if (ra == rb) continue;
                if (ra > rb) std::swap(ra, rb);
                auto& acc = pairs[{ra, rb}];
                ++acc.shared;
                if (posting[a].second != posting[b].second) ++acc.flips;
            }
        }
    }

    SignedDsu dsu(n);
    // Overlap-span proxy: each shared minimizer covers ~w bp; require both the
    // shared-count floor and the bp-span floor (configurable, not hardcoded).
    for (const auto& [pr, acc] : pairs) {
        if (acc.shared < cfg.min_shared) continue;
        if (acc.shared * w < cfg.min_overlap_bp) continue;
        const std::uint8_t flip = (acc.flips * 2 > acc.shared) ? 1 : 0;  // majority
        dsu.unite(pr.first, pr.second, flip);
    }

    // Dense cluster ids + per-read orientation relative to the cluster root.
    std::unordered_map<std::uint32_t, std::uint32_t> dense;
    for (std::size_t i = 0; i < n; ++i) {
        auto [root, parity] = dsu.find(static_cast<std::uint32_t>(i));
        auto it = dense.find(root);
        std::uint32_t cid;
        if (it == dense.end()) {
            cid = static_cast<std::uint32_t>(dense.size());
            dense.emplace(root, cid);
        } else {
            cid = it->second;
        }
        result[i] = ReadCluster{cid, parity != 0};
    }
    return result;
}

std::vector<std::vector<std::size_t>> GroupClusters(
    const std::vector<ReadCluster>& assignment) {
    std::uint32_t maxid = 0;
    for (const auto& a : assignment) maxid = std::max(maxid, a.cluster_id);
    std::vector<std::vector<std::size_t>> groups(assignment.empty() ? 0
                                                                    : maxid + 1);
    for (std::size_t i = 0; i < assignment.size(); ++i) {
        groups[assignment[i].cluster_id].push_back(i);
    }
    return groups;
}

}  // namespace llmap::mapper
