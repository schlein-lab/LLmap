#include "core/bucket_pyramid.h"

#include <fstream>
#include <stdexcept>

namespace llmap {

void BucketPyramid::add_l0_bucket(L0Bucket bucket) {
    bucket.id = static_cast<std::uint32_t>(l0_.size());
    l0_.push_back(std::move(bucket));
}

void BucketPyramid::add_l1_bucket(L1Bucket bucket, std::uint32_t parent_l0_id) {
    bucket.id = static_cast<std::uint32_t>(l1_.size());
    l1_.push_back(std::move(bucket));
    l1_to_l0_.push_back(parent_l0_id);
}

void BucketPyramid::add_l2_bucket(L2Bucket bucket, std::uint32_t parent_l1_id) {
    bucket.id = static_cast<std::uint32_t>(l2_.size());
    l2_.push_back(std::move(bucket));
    l2_to_l1_.push_back(parent_l1_id);
}

void BucketPyramid::set_biology_hint(std::uint32_t l2_bucket_id, BiologyHint hint) {
    biology_hints_[l2_bucket_id] = std::move(hint);
}

std::span<const L0Bucket> BucketPyramid::l0_buckets() const noexcept {
    return l0_;
}

std::span<const L1Bucket> BucketPyramid::l1_buckets() const noexcept {
    return l1_;
}

std::span<const L2Bucket> BucketPyramid::l2_buckets() const noexcept {
    return l2_;
}

std::uint32_t BucketPyramid::l1_parent(std::uint32_t l1_id) const {
    if (l1_id >= l1_to_l0_.size()) {
        throw std::out_of_range("l1_id out of range");
    }
    return l1_to_l0_[l1_id];
}

std::uint32_t BucketPyramid::l2_parent(std::uint32_t l2_id) const {
    if (l2_id >= l2_to_l1_.size()) {
        throw std::out_of_range("l2_id out of range");
    }
    return l2_to_l1_[l2_id];
}

std::optional<BiologyHint> BucketPyramid::get_biology_hint(std::uint32_t l2_bucket_id) const {
    auto it = biology_hints_.find(l2_bucket_id);
    if (it != biology_hints_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool BucketPyramid::validate() const noexcept {
    // Check l1→l0 indices valid
    for (std::size_t i = 0; i < l1_to_l0_.size(); ++i) {
        if (l1_to_l0_[i] >= l0_.size()) {
            return false;
        }
    }
    // Check l2→l1 indices valid
    for (std::size_t i = 0; i < l2_to_l1_.size(); ++i) {
        if (l2_to_l1_[i] >= l1_.size()) {
            return false;
        }
    }
    // Check size consistency
    if (l1_.size() != l1_to_l0_.size()) return false;
    if (l2_.size() != l2_to_l1_.size()) return false;

    return true;
}

namespace {

template<typename T>
void write_pod(std::ostream& os, const T& val) {
    os.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template<typename T>
T read_pod(std::istream& is) {
    T val;
    is.read(reinterpret_cast<char*>(&val), sizeof(T));
    return val;
}

void write_string(std::ostream& os, const std::string& s) {
    auto len = static_cast<std::uint32_t>(s.size());
    write_pod(os, len);
    os.write(s.data(), static_cast<std::streamsize>(len));
}

std::string read_string(std::istream& is) {
    auto len = read_pod<std::uint32_t>(is);
    std::string s(len, '\0');
    is.read(s.data(), static_cast<std::streamsize>(len));
    return s;
}

}  // namespace

void BucketPyramid::serialize(const std::filesystem::path& path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open file for writing: " + path.string());
    }

    // Header
    write_pod(ofs, MAGIC);
    write_pod(ofs, VERSION);

    // Counts
    write_pod(ofs, static_cast<std::uint32_t>(l0_.size()));
    write_pod(ofs, static_cast<std::uint32_t>(l1_.size()));
    write_pod(ofs, static_cast<std::uint32_t>(l2_.size()));
    write_pod(ofs, static_cast<std::uint32_t>(biology_hints_.size()));

    // L0 buckets
    for (const auto& b : l0_) {
        write_pod(ofs, b.id);
        write_string(ofs, b.name);
        write_pod(ofs, b.total_span);
    }

    // L1 buckets + parent indices
    for (const auto& b : l1_) {
        write_pod(ofs, b.id);
        write_string(ofs, b.target_id);
        write_pod(ofs, b.start);
        write_pod(ofs, b.end);
    }
    for (auto idx : l1_to_l0_) {
        write_pod(ofs, idx);
    }

    // L2 buckets + parent indices
    for (const auto& b : l2_) {
        write_pod(ofs, b.id);
        write_string(ofs, b.target_id);
        write_pod(ofs, b.start);
        write_pod(ofs, b.end);
    }
    for (auto idx : l2_to_l1_) {
        write_pod(ofs, idx);
    }

    // Biology hints
    for (const auto& [bucket_id, hint] : biology_hints_) {
        write_pod(ofs, bucket_id);
        write_pod(ofs, hint.prior_weight);
        bool has_annotation = hint.annotation.has_value();
        write_pod(ofs, has_annotation);
        if (has_annotation) {
            write_string(ofs, *hint.annotation);
        }
        bool has_partner = hint.paralog_partner_bucket.has_value();
        write_pod(ofs, has_partner);
        if (has_partner) {
            write_pod(ofs, *hint.paralog_partner_bucket);
        }
        write_pod(ofs, hint.expected_coverage_multiplier);
    }
}

BucketPyramid BucketPyramid::deserialize(const std::filesystem::path& path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading: " + path.string());
    }

    // Header
    auto magic = read_pod<std::uint32_t>(ifs);
    if (magic != MAGIC) {
        throw std::runtime_error("Invalid magic number in bucket pyramid file");
    }
    auto version = read_pod<std::uint32_t>(ifs);
    if (version != VERSION) {
        throw std::runtime_error("Unsupported bucket pyramid version: " + std::to_string(version));
    }

    // Counts
    auto l0_count = read_pod<std::uint32_t>(ifs);
    auto l1_count = read_pod<std::uint32_t>(ifs);
    auto l2_count = read_pod<std::uint32_t>(ifs);
    auto hint_count = read_pod<std::uint32_t>(ifs);

    BucketPyramid pyramid;

    // L0 buckets
    pyramid.l0_.reserve(l0_count);
    for (std::uint32_t i = 0; i < l0_count; ++i) {
        L0Bucket b;
        b.id = read_pod<std::uint32_t>(ifs);
        b.name = read_string(ifs);
        b.total_span = read_pod<std::uint64_t>(ifs);
        pyramid.l0_.push_back(std::move(b));
    }

    // L1 buckets
    pyramid.l1_.reserve(l1_count);
    for (std::uint32_t i = 0; i < l1_count; ++i) {
        L1Bucket b;
        b.id = read_pod<std::uint32_t>(ifs);
        b.target_id = read_string(ifs);
        b.start = read_pod<std::uint64_t>(ifs);
        b.end = read_pod<std::uint64_t>(ifs);
        pyramid.l1_.push_back(std::move(b));
    }
    pyramid.l1_to_l0_.reserve(l1_count);
    for (std::uint32_t i = 0; i < l1_count; ++i) {
        pyramid.l1_to_l0_.push_back(read_pod<std::uint32_t>(ifs));
    }

    // L2 buckets
    pyramid.l2_.reserve(l2_count);
    for (std::uint32_t i = 0; i < l2_count; ++i) {
        L2Bucket b;
        b.id = read_pod<std::uint32_t>(ifs);
        b.target_id = read_string(ifs);
        b.start = read_pod<std::uint64_t>(ifs);
        b.end = read_pod<std::uint64_t>(ifs);
        pyramid.l2_.push_back(std::move(b));
    }
    pyramid.l2_to_l1_.reserve(l2_count);
    for (std::uint32_t i = 0; i < l2_count; ++i) {
        pyramid.l2_to_l1_.push_back(read_pod<std::uint32_t>(ifs));
    }

    // Biology hints
    for (std::uint32_t i = 0; i < hint_count; ++i) {
        auto bucket_id = read_pod<std::uint32_t>(ifs);
        BiologyHint hint;
        hint.prior_weight = read_pod<float>(ifs);
        bool has_annotation = read_pod<bool>(ifs);
        if (has_annotation) {
            hint.annotation = read_string(ifs);
        }
        bool has_partner = read_pod<bool>(ifs);
        if (has_partner) {
            hint.paralog_partner_bucket = read_pod<std::uint32_t>(ifs);
        }
        hint.expected_coverage_multiplier = read_pod<float>(ifs);
        pyramid.biology_hints_[bucket_id] = std::move(hint);
    }

    return pyramid;
}

bool BucketPyramid::operator==(const BucketPyramid& other) const noexcept {
    if (l0_.size() != other.l0_.size()) return false;
    if (l1_.size() != other.l1_.size()) return false;
    if (l2_.size() != other.l2_.size()) return false;
    if (biology_hints_.size() != other.biology_hints_.size()) return false;

    // L0
    for (std::size_t i = 0; i < l0_.size(); ++i) {
        if (l0_[i].id != other.l0_[i].id) return false;
        if (l0_[i].name != other.l0_[i].name) return false;
        if (l0_[i].total_span != other.l0_[i].total_span) return false;
    }

    // L1
    for (std::size_t i = 0; i < l1_.size(); ++i) {
        if (l1_[i].id != other.l1_[i].id) return false;
        if (l1_[i].target_id != other.l1_[i].target_id) return false;
        if (l1_[i].start != other.l1_[i].start) return false;
        if (l1_[i].end != other.l1_[i].end) return false;
    }
    if (l1_to_l0_ != other.l1_to_l0_) return false;

    // L2
    for (std::size_t i = 0; i < l2_.size(); ++i) {
        if (l2_[i].id != other.l2_[i].id) return false;
        if (l2_[i].target_id != other.l2_[i].target_id) return false;
        if (l2_[i].start != other.l2_[i].start) return false;
        if (l2_[i].end != other.l2_[i].end) return false;
    }
    if (l2_to_l1_ != other.l2_to_l1_) return false;

    // Biology hints (unordered_map comparison)
    for (const auto& [k, v] : biology_hints_) {
        auto it = other.biology_hints_.find(k);
        if (it == other.biology_hints_.end()) return false;
        const auto& ov = it->second;
        if (v.prior_weight != ov.prior_weight) return false;
        if (v.annotation != ov.annotation) return false;
        if (v.paralog_partner_bucket != ov.paralog_partner_bucket) return false;
        if (v.expected_coverage_multiplier != ov.expected_coverage_multiplier) return false;
    }

    return true;
}

}  // namespace llmap
