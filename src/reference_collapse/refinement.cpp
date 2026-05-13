// LLmap — Refinement implementation.

#include "reference_collapse/refinement.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>

namespace llmap {

// ========== ChildIndex ==========

ChildIndex ChildIndex::BuildL0ToL1(const BucketPyramid& pyramid) {
    ChildIndex index;
    const std::size_t num_l0 = pyramid.l0_count();
    const std::size_t num_l1 = pyramid.l1_count();

    if (num_l0 == 0) {
        return index;
    }

    // Count children per L0 bucket
    std::vector<std::uint32_t> child_counts(num_l0, 0);
    for (std::size_t l1_id = 0; l1_id < num_l1; ++l1_id) {
        std::uint32_t parent = pyramid.l1_parent(static_cast<std::uint32_t>(l1_id));
        if (parent < num_l0) {
            ++child_counts[parent];
        }
    }

    // Build CSR offsets
    index.offsets_.resize(num_l0 + 1);
    index.offsets_[0] = 0;
    for (std::size_t i = 0; i < num_l0; ++i) {
        index.offsets_[i + 1] = index.offsets_[i] + child_counts[i];
    }

    // Fill children
    index.children_.resize(index.offsets_[num_l0]);
    std::vector<std::uint32_t> fill_pos = index.offsets_;
    fill_pos.resize(num_l0);

    for (std::size_t l1_id = 0; l1_id < num_l1; ++l1_id) {
        std::uint32_t parent = pyramid.l1_parent(static_cast<std::uint32_t>(l1_id));
        if (parent < num_l0) {
            index.children_[fill_pos[parent]++] = static_cast<std::uint32_t>(l1_id);
        }
    }

    return index;
}

ChildIndex ChildIndex::BuildL1ToL2(const BucketPyramid& pyramid) {
    ChildIndex index;
    const std::size_t num_l1 = pyramid.l1_count();
    const std::size_t num_l2 = pyramid.l2_count();

    if (num_l1 == 0) {
        return index;
    }

    // Count children per L1 bucket
    std::vector<std::uint32_t> child_counts(num_l1, 0);
    for (std::size_t l2_id = 0; l2_id < num_l2; ++l2_id) {
        std::uint32_t parent = pyramid.l2_parent(static_cast<std::uint32_t>(l2_id));
        if (parent < num_l1) {
            ++child_counts[parent];
        }
    }

    // Build CSR offsets
    index.offsets_.resize(num_l1 + 1);
    index.offsets_[0] = 0;
    for (std::size_t i = 0; i < num_l1; ++i) {
        index.offsets_[i + 1] = index.offsets_[i] + child_counts[i];
    }

    // Fill children
    index.children_.resize(index.offsets_[num_l1]);
    std::vector<std::uint32_t> fill_pos = index.offsets_;
    fill_pos.resize(num_l1);

    for (std::size_t l2_id = 0; l2_id < num_l2; ++l2_id) {
        std::uint32_t parent = pyramid.l2_parent(static_cast<std::uint32_t>(l2_id));
        if (parent < num_l1) {
            index.children_[fill_pos[parent]++] = static_cast<std::uint32_t>(l2_id);
        }
    }

    return index;
}

std::span<const std::uint32_t> ChildIndex::GetChildren(std::uint32_t parent_id) const {
    if (parent_id >= NumParents()) {
        return {};
    }
    const std::uint32_t start = offsets_[parent_id];
    const std::uint32_t end = offsets_[parent_id + 1];
    return std::span<const std::uint32_t>(children_.data() + start, end - start);
}

// ========== Refinement ==========

Refinement::Refinement(const RefinementConfig& config)
    : config_(config) {}

void Refinement::SetAiPrior(std::span<const float> ai_prior) {
    ai_prior_.assign(ai_prior.begin(), ai_prior.end());
}

std::uint32_t Refinement::CountRefinementCandidates(
    const WaveState& state,
    WaveLevel from_level) const {

    std::uint32_t count = 0;
    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r)) continue;
        if (state.get_level(r) != from_level) continue;

        // Check if any bucket is above expansion threshold
        auto probs = state.probabilities_for_read(r);
        for (float p : probs) {
            if (p >= config_.expansion_threshold) {
                ++count;
                break;
            }
        }
    }
    return count;
}

std::vector<BucketProb> Refinement::ExpandBucket(
    std::uint32_t parent_id,
    float parent_prob,
    const ChildIndex& child_index) const {

    auto children = child_index.GetChildren(parent_id);
    if (children.empty()) {
        return {};
    }

    std::vector<BucketProb> result;
    result.reserve(children.size());

    if (config_.weight_by_ai_prior && !ai_prior_.empty()) {
        // Weight by AI prior
        float total_prior = 0.0f;
        for (std::uint32_t child : children) {
            if (child < ai_prior_.size()) {
                total_prior += ai_prior_[child];
            } else {
                total_prior += 1.0f;
            }
        }

        if (total_prior > 0.0f) {
            for (std::uint32_t child : children) {
                float prior = (child < ai_prior_.size()) ? ai_prior_[child] : 1.0f;
                float child_prob = config_.preserve_relative_probs
                    ? parent_prob * (prior / total_prior)
                    : prior / total_prior;
                result.push_back({child, child_prob});
            }
        }
    } else {
        // Uniform distribution among children
        float child_prob = config_.preserve_relative_probs
            ? parent_prob / static_cast<float>(children.size())
            : 1.0f / static_cast<float>(children.size());

        for (std::uint32_t child : children) {
            result.push_back({child, child_prob});
        }
    }

    return result;
}

std::uint32_t Refinement::RefineRead(
    WaveState& state,
    std::uint32_t read_idx,
    WaveLevel from_level,
    const ChildIndex& child_index) {

    if (state.is_collapsed(read_idx)) return 0;
    if (state.get_level(read_idx) != from_level) return 0;

    auto buckets = state.bucket_indices_for_read(read_idx);
    auto probs = state.probabilities_for_read(read_idx);

    if (buckets.empty()) return 0;

    // Collect expanded candidates
    std::vector<BucketProb> new_candidates;
    new_candidates.reserve(config_.max_candidates);

    for (std::size_t i = 0; i < buckets.size(); ++i) {
        float prob = probs[i];
        if (prob < config_.expansion_threshold) continue;

        auto expanded = ExpandBucket(buckets[i], prob, child_index);
        for (const auto& bp : expanded) {
            new_candidates.push_back(bp);
        }
    }

    if (new_candidates.empty()) return 0;

    // Sort by bucket ID for CSR format
    std::sort(new_candidates.begin(), new_candidates.end(),
              [](const BucketProb& a, const BucketProb& b) {
                  return a.bucket_id < b.bucket_id;
              });

    // Merge duplicates (same child from multiple parents)
    std::vector<BucketProb> merged;
    merged.reserve(new_candidates.size());

    for (const auto& bp : new_candidates) {
        if (!merged.empty() && merged.back().bucket_id == bp.bucket_id) {
            merged.back().probability += bp.probability;
        } else {
            merged.push_back(bp);
        }
    }

    // Limit to max_candidates (keep highest probability)
    if (merged.size() > config_.max_candidates) {
        std::partial_sort(merged.begin(),
                          merged.begin() + config_.max_candidates,
                          merged.end(),
                          [](const BucketProb& a, const BucketProb& b) {
                              return a.probability > b.probability;
                          });
        merged.resize(config_.max_candidates);

        // Re-sort by bucket ID
        std::sort(merged.begin(), merged.end(),
                  [](const BucketProb& a, const BucketProb& b) {
                      return a.bucket_id < b.bucket_id;
                  });
    }

    // Normalize
    float total = std::accumulate(merged.begin(), merged.end(), 0.0f,
                                  [](float sum, const BucketProb& bp) {
                                      return sum + bp.probability;
                                  });
    if (total > 0.0f) {
        for (auto& bp : merged) {
            bp.probability /= total;
        }
    }

    // Update state
    state.set_read_candidates(read_idx, merged);
    state.set_level(read_idx, NextLevel(from_level));

    return static_cast<std::uint32_t>(merged.size());
}

RefinementStats Refinement::Refine(
    WaveState& state,
    WaveLevel from_level,
    const ChildIndex& child_index) {

    auto start_time = std::chrono::steady_clock::now();

    RefinementStats stats;
    std::uint64_t total_children = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.is_collapsed(r) || state.get_level(r) != from_level) {
            ++stats.reads_skipped;
            continue;
        }

        std::uint32_t children_added = RefineRead(state, r, from_level, child_index);
        if (children_added > 0) {
            ++stats.reads_refined;
            total_children += children_added;

            // Count buckets that were expanded
            auto old_buckets = state.bucket_indices_for_read(r);
            stats.buckets_expanded += static_cast<std::uint32_t>(old_buckets.size());
        } else {
            ++stats.reads_skipped;
        }
    }

    stats.total_children_added = static_cast<std::uint32_t>(total_children);
    if (stats.buckets_expanded > 0) {
        stats.avg_expansion_factor =
            static_cast<float>(total_children) / stats.buckets_expanded;
    }

    auto end_time = std::chrono::steady_clock::now();
    stats.refinement_time_ms =
        std::chrono::duration<float, std::milli>(end_time - start_time).count();

    return stats;
}

RefinementStats Refinement::RefineL0ToL1(
    WaveState& state,
    const ChildIndex& child_index) {
    return Refine(state, WaveLevel::L0, child_index);
}

RefinementStats Refinement::RefineL1ToL2(
    WaveState& state,
    const ChildIndex& child_index) {
    return Refine(state, WaveLevel::L1, child_index);
}

// ========== Free functions ==========

bool ShouldRefine(
    const WaveState& state,
    WaveLevel current_level,
    float collapsed_fraction_threshold,
    float high_prob_threshold) {

    if (state.n_reads() == 0) return false;

    std::uint32_t at_level = 0;
    std::uint32_t collapsed_or_high_prob = 0;

    for (std::uint32_t r = 0; r < state.n_reads(); ++r) {
        if (state.get_level(r) != current_level) continue;
        ++at_level;

        if (state.is_collapsed(r)) {
            ++collapsed_or_high_prob;
            continue;
        }

        // Check max probability
        auto probs = state.probabilities_for_read(r);
        if (!probs.empty()) {
            float max_prob = *std::max_element(probs.begin(), probs.end());
            if (max_prob >= high_prob_threshold) {
                ++collapsed_or_high_prob;
            }
        }
    }

    if (at_level == 0) return false;

    float fraction = static_cast<float>(collapsed_or_high_prob) / at_level;
    return fraction >= collapsed_fraction_threshold;
}

WaveLevel NextLevel(WaveLevel level) {
    switch (level) {
        case WaveLevel::L0: return WaveLevel::L1;
        case WaveLevel::L1: return WaveLevel::L2;
        case WaveLevel::L2: return WaveLevel::L3;
        case WaveLevel::L3: return WaveLevel::L3;
    }
    return WaveLevel::L3;
}

WaveLevel PrevLevel(WaveLevel level) {
    switch (level) {
        case WaveLevel::L0: return WaveLevel::L0;
        case WaveLevel::L1: return WaveLevel::L0;
        case WaveLevel::L2: return WaveLevel::L1;
        case WaveLevel::L3: return WaveLevel::L2;
    }
    return WaveLevel::L0;
}

}  // namespace llmap
