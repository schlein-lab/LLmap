// LLmap — Reference Index core implementation.
// Accessors and query methods.

#include "reference_collapse/reference_index.h"
#include "reference_collapse/reference_index_impl.h"

#include <algorithm>

namespace llmap {

// --- ReferenceIndex implementation ---

ReferenceIndex::ReferenceIndex(std::unique_ptr<ReferenceIndexImpl> impl)
    : impl_(std::move(impl)) {}

ReferenceIndex::~ReferenceIndex() = default;

ReferenceIndex::ReferenceIndex(ReferenceIndex&&) noexcept = default;
ReferenceIndex& ReferenceIndex::operator=(ReferenceIndex&&) noexcept = default;

const BucketPyramid& ReferenceIndex::Pyramid() const {
    return impl_->pyramid;
}

std::span<const ReferenceTarget> ReferenceIndex::Targets() const {
    return impl_->targets;
}

std::optional<ReferenceTarget> ReferenceIndex::FindTarget(
    std::string_view name) const {
    for (const auto& target : impl_->targets) {
        if (target.name == name) {
            return target;
        }
    }
    return std::nullopt;
}

std::size_t ReferenceIndex::NumTargets() const {
    return impl_->targets.size();
}

std::span<const float> ReferenceIndex::L0Embeddings() const {
    return impl_->l0_embeddings;
}

std::span<const float> ReferenceIndex::L1Embeddings() const {
    return impl_->l1_embeddings;
}

std::span<const float> ReferenceIndex::L2Embeddings() const {
    return impl_->l2_embeddings;
}

std::size_t ReferenceIndex::EmbeddingDim() const {
    return impl_->config.embedding_dim;
}

bool ReferenceIndex::HasEmbeddings() const {
    return !impl_->l0_embeddings.empty() &&
           !impl_->l1_embeddings.empty() &&
           !impl_->l2_embeddings.empty();
}

std::span<const float> ReferenceIndex::GetL0Embedding(std::uint32_t bucket_id) const {
    if (impl_->l0_embeddings.empty()) {
        return {};
    }
    auto dim = impl_->config.embedding_dim;
    auto offset = static_cast<std::size_t>(bucket_id) * dim;
    if (offset + dim > impl_->l0_embeddings.size()) {
        return {};
    }
    return std::span<const float>(
        impl_->l0_embeddings.data() + offset, dim);
}

std::span<const float> ReferenceIndex::GetL1Embedding(std::uint32_t bucket_id) const {
    if (impl_->l1_embeddings.empty()) {
        return {};
    }
    auto dim = impl_->config.embedding_dim;
    auto offset = static_cast<std::size_t>(bucket_id) * dim;
    if (offset + dim > impl_->l1_embeddings.size()) {
        return {};
    }
    return std::span<const float>(
        impl_->l1_embeddings.data() + offset, dim);
}

std::span<const float> ReferenceIndex::GetL2Embedding(std::uint32_t bucket_id) const {
    if (impl_->l2_embeddings.empty()) {
        return {};
    }
    auto dim = impl_->config.embedding_dim;
    auto offset = static_cast<std::size_t>(bucket_id) * dim;
    if (offset + dim > impl_->l2_embeddings.size()) {
        return {};
    }
    return std::span<const float>(
        impl_->l2_embeddings.data() + offset, dim);
}

std::optional<std::uint32_t> ReferenceIndex::FindL2Bucket(
    std::string_view target_name,
    std::uint64_t position) const {

    auto it = impl_->target_to_l2_buckets.find(std::string(target_name));
    if (it == impl_->target_to_l2_buckets.end()) {
        return std::nullopt;
    }

    const auto& bucket_ids = it->second;
    auto buckets = impl_->pyramid.l2_buckets();

    for (auto id : bucket_ids) {
        const auto& bucket = buckets[id];
        if (position >= bucket.start && position < bucket.end) {
            return id;
        }
    }
    return std::nullopt;
}

std::vector<std::uint32_t> ReferenceIndex::FindL2BucketsInRange(
    std::string_view target_name,
    std::uint64_t start,
    std::uint64_t end) const {

    std::vector<std::uint32_t> result;

    auto it = impl_->target_to_l2_buckets.find(std::string(target_name));
    if (it == impl_->target_to_l2_buckets.end()) {
        return result;
    }

    const auto& bucket_ids = it->second;
    auto buckets = impl_->pyramid.l2_buckets();

    for (auto id : bucket_ids) {
        const auto& bucket = buckets[id];
        // Check for overlap: [bucket.start, bucket.end) ∩ [start, end)
        if (bucket.start < end && start < bucket.end) {
            result.push_back(id);
        }
    }

    // Sort by start position
    std::sort(result.begin(), result.end(), [&](auto a, auto b) {
        return buckets[a].start < buckets[b].start;
    });

    return result;
}

const ReferenceIndexConfig& ReferenceIndex::Config() const {
    return impl_->config;
}

const ReferenceIndexStats& ReferenceIndex::Stats() const {
    return impl_->stats;
}

std::string ReferenceIndex::ReferenceVersion() const {
    return impl_->config.reference_version;
}

// --- Utility functions ---

std::uint64_t ComputeTotalLength(std::span<const ReferenceTarget> targets) {
    std::uint64_t total = 0;
    for (const auto& target : targets) {
        total += target.length;
    }
    return total;
}

std::pair<std::size_t, std::size_t> EstimateBucketCounts(
    std::uint64_t total_length,
    const ReferenceIndexConfig& config) {

    auto l1_count = (total_length + config.l1_granularity - 1) / config.l1_granularity;
    auto l2_count = (total_length + config.l2_granularity - 1) / config.l2_granularity;

    return {static_cast<std::size_t>(l1_count),
            static_cast<std::size_t>(l2_count)};
}

}  // namespace llmap
