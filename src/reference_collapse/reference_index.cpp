#include "reference_collapse/reference_index.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace llmap {

namespace {

constexpr std::uint32_t INDEX_MAGIC = 0x4C4C4958;  // "LLIX" (LLmap IndeX)
constexpr std::uint32_t INDEX_VERSION = 1;

template <typename T>
void WriteValue(std::ostream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
void ReadValue(std::istream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
}

void WriteString(std::ostream& out, const std::string& str) {
    auto size = static_cast<std::uint32_t>(str.size());
    WriteValue(out, size);
    out.write(str.data(), size);
}

std::string ReadString(std::istream& in) {
    std::uint32_t size = 0;
    ReadValue(in, size);
    std::string str(size, '\0');
    in.read(str.data(), size);
    return str;
}

void WriteVector(std::ostream& out, const std::vector<float>& vec) {
    auto size = static_cast<std::uint64_t>(vec.size());
    WriteValue(out, size);
    out.write(reinterpret_cast<const char*>(vec.data()),
              static_cast<std::streamsize>(size * sizeof(float)));
}

std::vector<float> ReadVector(std::istream& in) {
    std::uint64_t size = 0;
    ReadValue(in, size);
    std::vector<float> vec(size);
    in.read(reinterpret_cast<char*>(vec.data()),
            static_cast<std::streamsize>(size * sizeof(float)));
    return vec;
}

}  // namespace

// PIMPL implementation
class ReferenceIndexImpl {
public:
    ReferenceIndexConfig config;
    ReferenceIndexStats stats;
    BucketPyramid pyramid;
    std::vector<ReferenceTarget> targets;
    std::vector<float> l0_embeddings;
    std::vector<float> l1_embeddings;
    std::vector<float> l2_embeddings;

    // Spatial index for fast bucket lookup
    std::unordered_map<std::string, std::vector<std::uint32_t>> target_to_l2_buckets;

    void BuildSpatialIndex() {
        target_to_l2_buckets.clear();
        auto buckets = pyramid.l2_buckets();
        for (std::size_t i = 0; i < buckets.size(); ++i) {
            target_to_l2_buckets[buckets[i].target_id].push_back(
                static_cast<std::uint32_t>(i));
        }
    }
};

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

// --- ReferenceIndex implementation ---

ReferenceIndex::ReferenceIndex(std::unique_ptr<ReferenceIndexImpl> impl)
    : impl_(std::move(impl)) {}

ReferenceIndex::~ReferenceIndex() = default;

ReferenceIndex::ReferenceIndex(ReferenceIndex&&) noexcept = default;
ReferenceIndex& ReferenceIndex::operator=(ReferenceIndex&&) noexcept = default;

std::unique_ptr<ReferenceIndex> ReferenceIndex::Load(
    const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return nullptr;
    }

    std::uint32_t magic = 0;
    std::uint32_t version = 0;
    ReadValue(in, magic);
    ReadValue(in, version);

    if (magic != INDEX_MAGIC) {
        return nullptr;
    }
    if (version != INDEX_VERSION) {
        return nullptr;
    }

    auto impl = std::make_unique<ReferenceIndexImpl>();

    // Read config
    ReadValue(in, impl->config.l1_granularity);
    ReadValue(in, impl->config.l2_granularity);
    ReadValue(in, impl->config.embedding_dim);
    ReadValue(in, impl->config.include_embeddings);
    impl->config.reference_version = ReadString(in);

    // Read targets
    std::uint32_t num_targets = 0;
    ReadValue(in, num_targets);
    impl->targets.resize(num_targets);
    for (auto& target : impl->targets) {
        target.name = ReadString(in);
        ReadValue(in, target.length);
        target.md5 = ReadString(in);
    }

    // Read pyramid (temporary file path for deserialization)
    auto pyramid_path = path;
    pyramid_path.replace_extension(".pyramid.tmp");

    // Read pyramid size and data
    std::uint64_t pyramid_size = 0;
    ReadValue(in, pyramid_size);

    // Write pyramid data to temp file
    {
        std::ofstream pyramid_out(pyramid_path, std::ios::binary);
        std::vector<char> buffer(pyramid_size);
        in.read(buffer.data(), static_cast<std::streamsize>(pyramid_size));
        pyramid_out.write(buffer.data(), static_cast<std::streamsize>(pyramid_size));
    }

    impl->pyramid = BucketPyramid::deserialize(pyramid_path);
    std::filesystem::remove(pyramid_path);

    // Read embeddings
    impl->l0_embeddings = ReadVector(in);
    impl->l1_embeddings = ReadVector(in);
    impl->l2_embeddings = ReadVector(in);

    // Read stats
    ReadValue(in, impl->stats.num_targets);
    ReadValue(in, impl->stats.total_length);
    ReadValue(in, impl->stats.l0_buckets);
    ReadValue(in, impl->stats.l1_buckets);
    ReadValue(in, impl->stats.l2_buckets);
    ReadValue(in, impl->stats.build_time_seconds);
    ReadValue(in, impl->stats.embedding_time_seconds);
    ReadValue(in, impl->stats.index_size_bytes);

    // Rebuild spatial index
    impl->BuildSpatialIndex();

    return std::unique_ptr<ReferenceIndex>(
        new ReferenceIndex(std::move(impl)));
}

bool ReferenceIndex::Save(const std::filesystem::path& path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        return false;
    }

    WriteValue(out, INDEX_MAGIC);
    WriteValue(out, INDEX_VERSION);

    // Write config
    WriteValue(out, impl_->config.l1_granularity);
    WriteValue(out, impl_->config.l2_granularity);
    WriteValue(out, impl_->config.embedding_dim);
    WriteValue(out, impl_->config.include_embeddings);
    WriteString(out, impl_->config.reference_version);

    // Write targets
    auto num_targets = static_cast<std::uint32_t>(impl_->targets.size());
    WriteValue(out, num_targets);
    for (const auto& target : impl_->targets) {
        WriteString(out, target.name);
        WriteValue(out, target.length);
        WriteString(out, target.md5);
    }

    // Write pyramid (serialize to temp file, then embed)
    auto pyramid_path = path;
    pyramid_path.replace_extension(".pyramid.tmp");
    impl_->pyramid.serialize(pyramid_path);

    auto pyramid_size = static_cast<std::uint64_t>(
        std::filesystem::file_size(pyramid_path));
    WriteValue(out, pyramid_size);

    {
        std::ifstream pyramid_in(pyramid_path, std::ios::binary);
        std::vector<char> buffer(pyramid_size);
        pyramid_in.read(buffer.data(), static_cast<std::streamsize>(pyramid_size));
        out.write(buffer.data(), static_cast<std::streamsize>(pyramid_size));
    }
    std::filesystem::remove(pyramid_path);

    // Write embeddings
    WriteVector(out, impl_->l0_embeddings);
    WriteVector(out, impl_->l1_embeddings);
    WriteVector(out, impl_->l2_embeddings);

    // Write stats
    WriteValue(out, impl_->stats.num_targets);
    WriteValue(out, impl_->stats.total_length);
    WriteValue(out, impl_->stats.l0_buckets);
    WriteValue(out, impl_->stats.l1_buckets);
    WriteValue(out, impl_->stats.l2_buckets);
    WriteValue(out, impl_->stats.build_time_seconds);
    WriteValue(out, impl_->stats.embedding_time_seconds);
    WriteValue(out, impl_->stats.index_size_bytes);

    return out.good();
}

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
