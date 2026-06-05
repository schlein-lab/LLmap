// LLmap — Mode-6 Taxbin engine implementation. See taxbin.h.
#include "classify/taxbin.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace llmap::classify {

std::vector<float> Taxbin::RawLikelihood(const std::vector<float>& conf,
                                         float bind_threshold) {
    const std::size_t S = conf.size();
    std::vector<float> raw(S + 1, 0.0f);
    float cmax = 0.0f;
    for (std::size_t s = 0; s < S; ++s) {
        float c = conf[s] < 0.0f ? 0.0f : conf[s];
        raw[s] = c;
        cmax = std::max(cmax, c);
    }
    // Residual disbelief → NOVEL. When the best species is below the binding
    // threshold the read is, to that extent, unexplained by the panel.
    raw[S] = std::max(0.0f, bind_threshold - cmax);
    float sum = std::accumulate(raw.begin(), raw.end(), 0.0f);
    if (sum <= 0.0f) {
        // No evidence at all → entirely NOVEL.
        std::fill(raw.begin(), raw.end(), 0.0f);
        raw[S] = 1.0f;
        return raw;
    }
    for (float& v : raw) v /= sum;
    return raw;
}

namespace {

// argmax with deterministic tie-break (lowest index wins) and the runner-up,
// returning {top_idx, top_val, margin}.
struct TopTwo { int idx; float val; float margin; };
TopTwo ArgTop2(const std::vector<float>& v) {
    int best = 0; float bestv = v.empty() ? 0.0f : v[0];
    float second = -1.0f;
    for (std::size_t i = 1; i < v.size(); ++i) {
        if (v[i] > bestv) { second = bestv; bestv = v[i]; best = static_cast<int>(i); }
        else if (v[i] > second) { second = v[i]; }
    }
    if (second < 0.0f) second = 0.0f;
    return {best, bestv, bestv - second};
}

}  // namespace

TaxbinResult Taxbin::Run(const TaxbinInput& in) const {
    TaxbinResult out;
    const std::size_t N = in.read_ids.size();
    const std::size_t S = in.species_labels.size();
    const int NOVEL = static_cast<int>(S);
    out.num_species = S;

    const bool have_clusters =
        cfg_.enable_cluster_collapse && in.cluster_ids.size() == N && N > 0;

    // 1. Per-read raw likelihoods over (species..., NOVEL).
    std::vector<std::vector<float>> own(N);
    for (std::size_t r = 0; r < N; ++r) {
        const std::vector<float>& conf =
            (r < in.per_read_species_conf.size()) ? in.per_read_species_conf[r]
                                                   : std::vector<float>(S, 0.0f);
        own[r] = RawLikelihood(conf, cfg_.bind_threshold);
    }

    // 2. Cluster aggregate likelihoods (sum of member raw vectors → normalize).
    //    Summing pre-normalized-per-read vectors weights every read equally;
    //    a cluster's collective mass decides its consensus even when single
    //    reads are individually ambiguous.
    std::unordered_map<std::uint32_t, std::vector<float>> cluster_agg;
    std::unordered_map<std::uint32_t, std::size_t> cluster_size;
    if (have_clusters) {
        for (std::size_t r = 0; r < N; ++r) {
            auto cid = in.cluster_ids[r];
            auto& agg = cluster_agg[cid];
            if (agg.empty()) agg.assign(S + 1, 0.0f);
            for (std::size_t k = 0; k <= S; ++k) agg[k] += own[r][k];
            cluster_size[cid]++;
        }
        for (auto& [cid, agg] : cluster_agg) {
            float sum = std::accumulate(agg.begin(), agg.end(), 0.0f);
            if (sum > 0.0f) for (float& v : agg) v /= sum;
        }
    }

    // 3. Per-read posterior = (1-w)*own + w*cluster_aggregate, normalized.
    //    Then the hard call is argmax of the posterior; NOVEL is just index S.
    out.reads.reserve(N);
    for (std::size_t r = 0; r < N; ++r) {
        TaxbinReadResult rr;
        rr.read_id = in.read_ids[r];
        rr.cluster_id = have_clusters ? in.cluster_ids[r] : 0u;

        std::vector<float> post = own[r];
        if (have_clusters) {
            const auto& agg = cluster_agg[in.cluster_ids[r]];
            const float w = cfg_.cluster_weight;
            for (std::size_t k = 0; k <= S; ++k)
                post[k] = (1.0f - w) * own[r][k] + w * agg[k];
            float sum = std::accumulate(post.begin(), post.end(), 0.0f);
            if (sum > 0.0f) for (float& v : post) v /= sum;
        }
        rr.likelihood = post;

        TopTwo own_top = ArgTop2(own[r]);
        TopTwo post_top = ArgTop2(post);
        rr.top_prob = post_top.val;
        rr.margin = post_top.margin;

        // A confident species call needs (a) argmax != NOVEL, (b) enough margin.
        if (post_top.idx == NOVEL || post_top.margin < cfg_.min_margin) {
            // Ambiguous or unexplained.
            if (post_top.idx == NOVEL) { rr.novel = true; rr.top_species = -1; }
            else { rr.top_species = post_top.idx; }  // low-margin species call retained but flagged via margin
        } else {
            rr.top_species = post_top.idx;
        }
        if (rr.top_species == -1 && post_top.idx != NOVEL) rr.top_species = post_top.idx;

        rr.novel = (post_top.idx == NOVEL);
        rr.by_cluster = (!rr.novel) && (own_top.idx != post_top.idx);
        if (rr.novel) out.num_novel++;
        if (rr.by_cluster) out.num_collapsed_by_cluster++;
        out.reads.push_back(std::move(rr));
    }

    // 4. Cluster summary (consensus species + purity).
    if (have_clusters) {
        for (auto& [cid, agg] : cluster_agg) {
            TaxbinClusterResult cr;
            cr.cluster_id = cid;
            cr.size = cluster_size[cid];
            cr.aggregate_likelihood = agg;
            TopTwo t = ArgTop2(agg);
            cr.consensus_species = (t.idx == NOVEL) ? -1 : t.idx;
            std::size_t agree = 0;
            for (std::size_t r = 0; r < N; ++r) {
                if (in.cluster_ids[r] != cid) continue;
                TopTwo ot = ArgTop2(own[r]);
                int own_call = (ot.idx == NOVEL) ? -1 : ot.idx;
                if (own_call == cr.consensus_species) agree++;
            }
            cr.purity = cr.size ? static_cast<float>(agree) / static_cast<float>(cr.size) : 0.0f;
            out.clusters.push_back(std::move(cr));
        }
        std::sort(out.clusters.begin(), out.clusters.end(),
                  [](const TaxbinClusterResult& a, const TaxbinClusterResult& b) {
                      return a.cluster_id < b.cluster_id;
                  });
    }
    return out;
}

}  // namespace llmap::classify
