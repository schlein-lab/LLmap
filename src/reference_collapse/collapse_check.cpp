// LLmap — Collapse Check implementation.

#include "reference_collapse/collapse_check.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace llmap {

CollapseChecker::CollapseChecker(const CollapseCheckConfig& config)
    : config_(config)
{
}

ReadCollapseInfo CollapseChecker::AnalyzeRead(
    const WaveState& state,
    std::uint32_t read_idx) const
{
    ReadCollapseInfo info;
    info.read_idx = read_idx;

    const auto probs = state.probabilities_for_read(read_idx);
    const auto buckets = state.bucket_indices_for_read(read_idx);

    // Handle empty distribution
    if (probs.empty()) {
        info.dropout_type = DropoutType::Orphaned;
        return info;
    }

    // Find top-2 probabilities and their buckets
    float max_p = 0.0f;
    float second_p = 0.0f;
    std::uint32_t max_bucket = 0;
    std::uint32_t second_bucket = 0;

    for (std::size_t i = 0; i < probs.size(); ++i) {
        if (probs[i] > max_p) {
            second_p = max_p;
            second_bucket = max_bucket;
            max_p = probs[i];
            max_bucket = buckets[i];
        } else if (probs[i] > second_p) {
            second_p = probs[i];
            second_bucket = buckets[i];
        }
    }

    info.max_probability = max_p;
    info.second_probability = second_p;
    info.top_bucket = max_bucket;
    info.second_bucket = second_bucket;

    // Compute entropy
    info.entropy = ComputeReadEntropy(probs);

    // Classify dropout type
    info.dropout_type = ClassifyDropout(info);

    return info;
}

DropoutType CollapseChecker::ClassifyDropout(const ReadCollapseInfo& info) const {
    // Orphaned: no candidates or max_p is essentially zero
    if (info.max_probability < 1e-9f) {
        return DropoutType::Orphaned;
    }

    // Converged: max probability above collapse threshold
    if (info.max_probability >= config_.tau_collapse) {
        return DropoutType::Converged;
    }

    // Ambiguous: gap between top-2 is too small
    const float gap = info.max_probability - info.second_probability;
    if (gap < config_.tau_ambiguous && info.second_probability > 0.0f) {
        return DropoutType::Ambiguous;
    }

    // Low confidence: max probability too low
    if (info.max_probability < config_.tau_low_confidence) {
        return DropoutType::LowConfidence;
    }

    // Still active
    return DropoutType::None;
}

CollapseCheckStats CollapseChecker::Check(WaveState& state) {
    CollapseCheckStats stats;
    stats.total_reads = state.n_reads();

    float sum_max_prob = 0.0f;
    float sum_entropy = 0.0f;
    std::uint32_t analyzed_count = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) {
            ++stats.already_collapsed;
            continue;
        }

        const auto info = AnalyzeRead(state, r);
        ++analyzed_count;

        // Update probability stats
        sum_max_prob += info.max_probability;
        sum_entropy += info.entropy;
        stats.min_max_prob = std::min(stats.min_max_prob, info.max_probability);
        stats.max_max_prob = std::max(stats.max_max_prob, info.max_probability);

        switch (info.dropout_type) {
            case DropoutType::Converged:
                if (config_.auto_collapse) {
                    state.collapse_read(r);
                }
                ++stats.newly_collapsed;
                break;

            case DropoutType::Ambiguous:
                ++stats.ambiguous;
                break;

            case DropoutType::LowConfidence:
                ++stats.low_confidence;
                break;

            case DropoutType::Orphaned:
                ++stats.orphaned;
                break;

            case DropoutType::None:
            default:
                ++stats.active_after;
                break;
        }
    }

    // Compute averages
    if (analyzed_count > 0) {
        stats.mean_max_prob = sum_max_prob / analyzed_count;
        stats.mean_entropy = sum_entropy / analyzed_count;
    }

    // Fix min for edge case
    if (analyzed_count == 0) {
        stats.min_max_prob = 0.0f;
    }

    return stats;
}

CollapseCheckStats CollapseChecker::CheckDetailed(
    WaveState& state,
    std::vector<ReadCollapseInfo>& read_info)
{
    CollapseCheckStats stats;
    stats.total_reads = state.n_reads();

    read_info.clear();
    read_info.reserve(state.n_reads());

    float sum_max_prob = 0.0f;
    float sum_entropy = 0.0f;
    std::uint32_t analyzed_count = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) {
            ++stats.already_collapsed;
            // Add placeholder for collapsed reads
            ReadCollapseInfo info;
            info.read_idx = r;
            info.dropout_type = DropoutType::Converged;
            info.max_probability = 1.0f;
            info.top_bucket = state.collapsed_bucket(r);
            read_info.push_back(info);
            continue;
        }

        auto info = AnalyzeRead(state, r);
        read_info.push_back(info);
        ++analyzed_count;

        // Update probability stats
        sum_max_prob += info.max_probability;
        sum_entropy += info.entropy;
        stats.min_max_prob = std::min(stats.min_max_prob, info.max_probability);
        stats.max_max_prob = std::max(stats.max_max_prob, info.max_probability);

        switch (info.dropout_type) {
            case DropoutType::Converged:
                if (config_.auto_collapse) {
                    state.collapse_read(r);
                }
                ++stats.newly_collapsed;
                break;

            case DropoutType::Ambiguous:
                ++stats.ambiguous;
                break;

            case DropoutType::LowConfidence:
                ++stats.low_confidence;
                break;

            case DropoutType::Orphaned:
                ++stats.orphaned;
                break;

            case DropoutType::None:
            default:
                ++stats.active_after;
                break;
        }
    }

    // Compute averages
    if (analyzed_count > 0) {
        stats.mean_max_prob = sum_max_prob / analyzed_count;
        stats.mean_entropy = sum_entropy / analyzed_count;
    }

    // Fix min for edge case
    if (analyzed_count == 0) {
        stats.min_max_prob = 0.0f;
    }

    return stats;
}

// Convenience functions

CollapseCheckStats CheckAndCollapse(
    WaveState& state,
    float tau_collapse)
{
    CollapseCheckConfig config;
    config.tau_collapse = tau_collapse;
    config.auto_collapse = true;

    CollapseChecker checker(config);
    return checker.Check(state);
}

std::vector<DropoutType> GetDropoutTypes(
    const WaveState& state,
    float tau_collapse,
    float tau_ambiguous)
{
    CollapseCheckConfig config;
    config.tau_collapse = tau_collapse;
    config.tau_ambiguous = tau_ambiguous;
    config.auto_collapse = false;

    CollapseChecker checker(config);

    std::vector<DropoutType> types;
    types.reserve(state.n_reads());

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) {
            types.push_back(DropoutType::Converged);
        } else {
            const auto info = checker.AnalyzeRead(state, r);
            types.push_back(info.dropout_type);
        }
    }

    return types;
}

float ComputeReadEntropy(std::span<const float> probs) {
    float entropy = 0.0f;
    for (const float p : probs) {
        if (p > 1e-10f) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

std::vector<std::uint32_t> GetReadsByDropoutType(
    const std::vector<DropoutType>& types,
    DropoutType target)
{
    std::vector<std::uint32_t> indices;
    for (std::uint32_t i = 0; i < types.size(); ++i) {
        if (types[i] == target) {
            indices.push_back(i);
        }
    }
    return indices;
}

}  // namespace llmap
