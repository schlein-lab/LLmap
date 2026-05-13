// LLmap — Reference Index Builder implementation.
// Builder pattern for constructing ReferenceIndex from targets and embeddings.

#include "reference_collapse/reference_index.h"
#include "reference_collapse/reference_index_impl.h"

#include <algorithm>
#include <chrono>

namespace llmap {

// --- Builder implementation ---

ReferenceIndex::Builder::Builder(const ReferenceIndexConfig& config)
    : config_(config) {}

ReferenceIndex::Builder::~Builder() = default;

ReferenceIndex::Builder::Builder(Builder&&) noexcept = default;
ReferenceIndex::Builder& ReferenceIndex::Builder::operator=(Builder&&) noexcept = default;

ReferenceIndex::Builder& ReferenceIndex::Builder::AddTarget(
    const ReferenceTarget& target) {
    targets_.push_back(target);
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::AddTargets(
    std::span<const ReferenceTarget> targets) {
    targets_.insert(targets_.end(), targets.begin(), targets.end());
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::SetL0Embeddings(
    std::span<const float> embeddings) {
    l0_embeddings_.assign(embeddings.begin(), embeddings.end());
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::SetL1Embeddings(
    std::span<const float> embeddings) {
    l1_embeddings_.assign(embeddings.begin(), embeddings.end());
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::SetL2Embeddings(
    std::span<const float> embeddings) {
    l2_embeddings_.assign(embeddings.begin(), embeddings.end());
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::SetBiologyPrior(
    const std::filesystem::path& /*json_path*/) {
    // TODO: Parse Claude Session A JSON in Phase 7
    return *this;
}

ReferenceIndex::Builder& ReferenceIndex::Builder::SetBiologyHint(
    std::uint32_t l2_bucket_id,
    const BiologyHint& hint) {
    biology_hints_[l2_bucket_id] = hint;
    return *this;
}

std::unique_ptr<ReferenceIndex> ReferenceIndex::Builder::Build() {
    auto start = std::chrono::steady_clock::now();

    if (targets_.empty()) {
        last_error_ = "No reference targets provided";
        return nullptr;
    }

    auto impl = std::make_unique<ReferenceIndexImpl>();
    impl->config = config_;
    impl->targets = targets_;

    // Build bucket pyramid from targets
    std::uint32_t l0_id = 0;
    std::uint32_t l1_id = 0;
    std::uint32_t l2_id = 0;

    for (const auto& target : targets_) {
        // L0: one bucket per chromosome/target
        L0Bucket l0_bucket;
        l0_bucket.id = l0_id;
        l0_bucket.name = target.name;
        l0_bucket.total_span = target.length;
        impl->pyramid.add_l0_bucket(l0_bucket);

        // L1: 5 MB windows
        std::uint64_t pos = 0;
        while (pos < target.length) {
            std::uint64_t end = std::min(pos + config_.l1_granularity, target.length);
            L1Bucket l1_bucket;
            l1_bucket.id = l1_id;
            l1_bucket.target_id = target.name;
            l1_bucket.start = pos;
            l1_bucket.end = end;
            impl->pyramid.add_l1_bucket(l1_bucket, l0_id);

            // L2: 50 kb windows within this L1
            std::uint64_t l2_pos = pos;
            while (l2_pos < end) {
                std::uint64_t l2_end = std::min(l2_pos + config_.l2_granularity, end);
                L2Bucket l2_bucket;
                l2_bucket.id = l2_id;
                l2_bucket.target_id = target.name;
                l2_bucket.start = l2_pos;
                l2_bucket.end = l2_end;
                impl->pyramid.add_l2_bucket(l2_bucket, l1_id);

                l2_pos = l2_end;
                ++l2_id;
            }

            pos = end;
            ++l1_id;
        }
        ++l0_id;
    }

    // Set biology hints
    for (const auto& [bucket_id, hint] : biology_hints_) {
        impl->pyramid.set_biology_hint(bucket_id, hint);
    }

    // Copy embeddings
    impl->l0_embeddings = std::move(l0_embeddings_);
    impl->l1_embeddings = std::move(l1_embeddings_);
    impl->l2_embeddings = std::move(l2_embeddings_);

    // Build spatial index
    impl->BuildSpatialIndex();

    // Validate pyramid
    if (!impl->pyramid.validate()) {
        last_error_ = "Bucket pyramid validation failed";
        return nullptr;
    }

    // Compute stats
    auto end = std::chrono::steady_clock::now();
    stats_.num_targets = targets_.size();
    stats_.total_length = ComputeTotalLength(targets_);
    stats_.l0_buckets = impl->pyramid.l0_count();
    stats_.l1_buckets = impl->pyramid.l1_count();
    stats_.l2_buckets = impl->pyramid.l2_count();
    stats_.build_time_seconds = std::chrono::duration<float>(end - start).count();

    impl->stats = stats_;

    return std::unique_ptr<ReferenceIndex>(
        new ReferenceIndex(std::move(impl)));
}

}  // namespace llmap
