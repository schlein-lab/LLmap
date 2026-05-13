// LLmap — Reference Index I/O operations.
// Serialization helpers and Load/Save implementations.

#include "reference_collapse/reference_index.h"
#include "reference_collapse/reference_index_impl.h"

#include <cstring>
#include <fstream>

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

}  // namespace llmap
