// LLmap — EM Iterator CPU implementation.

#include "reference_collapse/em_iterator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

namespace llmap {

std::span<const std::uint32_t>
BucketNeighborhood::GetNeighbors(std::uint32_t bucket_id) const {
    if (bucket_id >= NumBuckets()) {
        return {};
    }
    const auto start = offsets[bucket_id];
    const auto end = offsets[bucket_id + 1];
    return {neighbors.data() + start, end - start};
}

std::span<const float>
BucketNeighborhood::GetWeights(std::uint32_t bucket_id) const {
    if (bucket_id >= NumBuckets()) {
        return {};
    }
    const auto start = offsets[bucket_id];
    const auto end = offsets[bucket_id + 1];
    return {weights.data() + start, end - start};
}

BucketNeighborhood BuildNeighborhood(
    std::span<const std::uint64_t> positions,
    const SmoothingKernelConfig& config)
{
    BucketNeighborhood hood;
    const auto n = positions.size();
    if (n == 0) {
        hood.offsets.push_back(0);
        return hood;
    }

    hood.offsets.reserve(n + 1);
    hood.offsets.push_back(0);

    const float inv_2sigma_sq = 1.0f / (2.0f * config.sigma_genome_bp * config.sigma_genome_bp);
    const auto max_dist = static_cast<std::uint64_t>(
        config.sigma_genome_bp * std::sqrt(-2.0f * std::log(config.min_weight)));

    for (std::size_t b = 0; b < n; ++b) {
        const auto pos_b = positions[b];
        std::vector<std::pair<std::uint32_t, float>> candidates;

        for (std::size_t b2 = 0; b2 < n; ++b2) {
            if (b == b2) continue;

            const auto pos_b2 = positions[b2];
            const auto dist = (pos_b > pos_b2) ? (pos_b - pos_b2) : (pos_b2 - pos_b);

            if (dist > max_dist) continue;

            const float dist_f = static_cast<float>(dist);
            const float weight = std::exp(-dist_f * dist_f * inv_2sigma_sq);

            if (weight >= config.min_weight) {
                candidates.emplace_back(static_cast<std::uint32_t>(b2), weight);
            }
        }

        std::sort(candidates.begin(), candidates.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        const auto keep = std::min(candidates.size(),
            static_cast<std::size_t>(config.max_neighbors));

        for (std::size_t i = 0; i < keep; ++i) {
            hood.neighbors.push_back(candidates[i].first);
            hood.weights.push_back(candidates[i].second);
        }

        hood.offsets.push_back(static_cast<std::uint32_t>(hood.neighbors.size()));
    }

    return hood;
}

EmIterator::EmIterator(const EmIteratorConfig& config)
    : config_(config)
{
}

EmIterator::~EmIterator() = default;

EmIterator::EmIterator(EmIterator&&) noexcept = default;
EmIterator& EmIterator::operator=(EmIterator&&) noexcept = default;

void EmIterator::SetLikelihoodFn(LikelihoodFn fn) {
    likelihood_fn_ = std::move(fn);
}

void EmIterator::SetAiPriorFn(AiPriorFn fn) {
    ai_prior_fn_ = std::move(fn);
}

void EmIterator::SetBiologyPrior(std::span<const float> bio_prior) {
    bio_prior_.assign(bio_prior.begin(), bio_prior.end());
}

void EmIterator::SetNeighborhood(const BucketNeighborhood& neighborhood) {
    neighborhood_ = &neighborhood;
}

std::vector<float> EmIterator::ComputeCoveragePrior(
    const WaveState& state,
    std::uint32_t num_buckets) const
{
    std::vector<float> lambda(num_buckets, 0.0f);

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) continue;

        const auto buckets = state.bucket_indices_for_read(r);
        const auto probs = state.probabilities_for_read(r);

        for (std::size_t i = 0; i < buckets.size(); ++i) {
            if (buckets[i] < num_buckets) {
                lambda[buckets[i]] += probs[i];
            }
        }
    }

    return lambda;
}

float EmIterator::ComputeNeighborContribution(
    const WaveState& state,
    std::uint32_t read_idx,
    std::uint32_t bucket_id) const
{
    if (!neighborhood_ || !config_.apply_smoothing) {
        return 1.0f;
    }

    const auto neighbors = neighborhood_->GetNeighbors(bucket_id);
    const auto weights = neighborhood_->GetWeights(bucket_id);

    if (neighbors.empty()) {
        return 1.0f;
    }

    float sum = 0.0f;
    for (std::size_t i = 0; i < neighbors.size(); ++i) {
        const float p = state.get_probability(read_idx, neighbors[i]);
        sum += weights[i] * p;
    }

    return 1.0f + sum;
}

void EmIterator::ComputeReadUpdate(
    const WaveState& state,
    std::uint32_t read_idx,
    const std::vector<float>& coverage_prior,
    std::vector<float>& out_probs) const
{
    const auto buckets = state.bucket_indices_for_read(read_idx);
    const auto old_probs = state.probabilities_for_read(read_idx);

    out_probs.resize(buckets.size());

    for (std::size_t i = 0; i < buckets.size(); ++i) {
        const auto b = buckets[i];
        const float p_old = old_probs[i];

        float L = 1.0f;
        if (likelihood_fn_) {
            L = likelihood_fn_(read_idx, b);
        }
        L = std::pow(L, config_.weight_seq_likelihood);

        float lambda = 1.0f;
        if (b < coverage_prior.size()) {
            lambda = std::max(coverage_prior[b], 0.001f);
        }
        lambda = std::pow(lambda, config_.weight_coverage);

        float pi_ai = 1.0f;
        if (ai_prior_fn_) {
            pi_ai = ai_prior_fn_(read_idx, b);
        }
        pi_ai = std::pow(pi_ai, config_.weight_ai_prior);

        float pi_bio = 1.0f;
        if (b < bio_prior_.size()) {
            pi_bio = bio_prior_[b];
        }
        pi_bio = std::pow(pi_bio, config_.weight_bio_prior);

        const float K_contrib = ComputeNeighborContribution(state, read_idx, b);

        const float update_term = L * lambda * pi_ai * pi_bio * K_contrib;

        out_probs[i] = (1.0f - config_.gamma) * p_old + config_.gamma * update_term;
    }

    NormalizeProbabilities(out_probs);
}

void EmIterator::ApplyKSmoothing(
    WaveState& state,
    const std::vector<float>& /*coverage_prior*/) const
{
    if (!neighborhood_ || !config_.apply_smoothing) {
        return;
    }

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) continue;

        const auto buckets = state.bucket_indices_for_read(r);
        const auto probs = state.probabilities_for_read(r);

        std::vector<float> smoothed(buckets.size());

        for (std::size_t i = 0; i < buckets.size(); ++i) {
            float contrib = ComputeNeighborContribution(state, r, buckets[i]);
            smoothed[i] = probs[i] * contrib;
        }

        NormalizeProbabilities(smoothed);

        for (std::size_t i = 0; i < buckets.size(); ++i) {
            state.update_probability(r, buckets[i], smoothed[i]);
        }
    }
}

EmIterationStats EmIterator::Step(
    WaveState& state,
    std::uint32_t num_buckets)
{
    auto start = std::chrono::high_resolution_clock::now();

    EmIterationStats stats;
    stats.reads_processed = 0;
    stats.reads_collapsed = 0;
    stats.max_prob_delta = 0.0f;

    auto coverage = ComputeCoveragePrior(state, num_buckets);

    std::vector<float> new_probs;
    float entropy_sum = 0.0f;
    std::uint32_t entropy_count = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) {
            ++stats.reads_collapsed;
            continue;
        }

        ++stats.reads_processed;

        const auto old_probs = state.probabilities_for_read(r);
        ComputeReadUpdate(state, r, coverage, new_probs);

        float max_delta = 0.0f;
        const auto buckets = state.bucket_indices_for_read(r);
        for (std::size_t i = 0; i < buckets.size(); ++i) {
            const float delta = std::abs(new_probs[i] - old_probs[i]);
            max_delta = std::max(max_delta, delta);
            state.update_probability(r, buckets[i], new_probs[i]);
        }

        stats.max_prob_delta = std::max(stats.max_prob_delta, max_delta);

        entropy_sum += ComputeEntropy(new_probs);
        ++entropy_count;
    }

    stats.reads_active = stats.reads_processed;
    stats.mean_entropy = (entropy_count > 0) ? (entropy_sum / entropy_count) : 0.0f;

    auto end = std::chrono::high_resolution_clock::now();
    stats.iteration_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return stats;
}

std::uint32_t EmIterator::CheckAndCollapse(WaveState& state) {
    std::uint32_t newly_collapsed = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) continue;

        if (state.check_collapse(r, config_.tau_collapse)) {
            state.collapse_read(r);
            ++newly_collapsed;
        }
    }

    return newly_collapsed;
}

float ComputeEntropy(std::span<const float> probs) {
    float entropy = 0.0f;
    for (const float p : probs) {
        if (p > 1e-10f) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

void NormalizeProbabilities(std::span<float> probs) {
    const float sum = std::accumulate(probs.begin(), probs.end(), 0.0f);
    if (sum > 1e-10f) {
        const float inv_sum = 1.0f / sum;
        for (float& p : probs) {
            p *= inv_sum;
        }
    } else if (!probs.empty()) {
        const float uniform = 1.0f / static_cast<float>(probs.size());
        for (float& p : probs) {
            p = uniform;
        }
    }
}

}  // namespace llmap
