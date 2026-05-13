#include "classical/minimizer_index.h"
#include "classical/minimizer_index_impl.h"

#include <algorithm>
#include <limits>

namespace llmap::classical {

MinimizerIndex::MinimizerIndex() : impl_(std::make_unique<Impl>()) {}

MinimizerIndex::MinimizerIndex(const MinimizerConfig& config)
    : config_(config), impl_(std::make_unique<Impl>()) {}

MinimizerIndex::~MinimizerIndex() = default;

MinimizerIndex::MinimizerIndex(MinimizerIndex&&) noexcept = default;
MinimizerIndex& MinimizerIndex::operator=(MinimizerIndex&&) noexcept = default;

// Builder implementation
MinimizerIndex::Builder::Builder(const MinimizerConfig& config)
    : impl_(std::make_unique<Impl>()) {
    impl_->config = config;
}

MinimizerIndex::Builder::~Builder() = default;

MinimizerIndex::Builder::Builder(Builder&&) noexcept = default;
MinimizerIndex::Builder& MinimizerIndex::Builder::operator=(Builder&&) noexcept = default;

MinimizerIndex::Builder& MinimizerIndex::Builder::AddSequence(
    std::string_view name,
    std::string_view sequence)
{
    impl_->names.emplace_back(name);
    impl_->sequences.emplace_back(sequence);
    return *this;
}

MinimizerIndex::Builder& MinimizerIndex::Builder::AddSequences(
    std::span<const std::string> names,
    std::span<const std::string> sequences)
{
    for (size_t i = 0; i < names.size() && i < sequences.size(); ++i) {
        impl_->names.push_back(names[i]);
        impl_->sequences.push_back(sequences[i]);
    }
    return *this;
}

std::unique_ptr<MinimizerIndex> MinimizerIndex::Builder::Build() {
    auto index = std::make_unique<MinimizerIndex>(impl_->config);

    // Extract and index minimizers from all sequences
    for (size_t ref_id = 0; ref_id < impl_->sequences.size(); ++ref_id) {
        const auto& seq = impl_->sequences[ref_id];
        const auto& name = impl_->names[ref_id];

        // Store sequence info
        IndexedSequence info;
        info.name = name;
        info.length = static_cast<uint32_t>(seq.size());
        index->impl_->sequences.push_back(info);

        // Extract minimizers
        auto minimizers = ExtractMinimizers(seq, impl_->config);
        index->impl_->stats.total_kmers += seq.size() - impl_->config.k + 1;
        index->impl_->stats.unique_minimizers += minimizers.size();

        // Add to index
        for (const auto& m : minimizers) {
            index->impl_->index[m.hash].emplace_back(
                static_cast<uint32_t>(ref_id),
                m.pos);
        }
    }

    // Compute average spacing
    if (!impl_->sequences.empty()) {
        size_t total_len = 0;
        for (const auto& seq : impl_->sequences) {
            total_len += seq.size();
        }
        if (index->impl_->stats.unique_minimizers > 0) {
            index->impl_->stats.avg_minimizer_spacing =
                static_cast<float>(total_len) /
                static_cast<float>(index->impl_->stats.unique_minimizers);
        }
    }

    return index;
}

std::vector<MinimizerHit> MinimizerIndex::Query(std::string_view sequence) const {
    return Query(sequence, std::numeric_limits<size_t>::max());
}

std::vector<MinimizerHit> MinimizerIndex::Query(
    std::string_view sequence,
    size_t max_hits) const
{
    std::vector<MinimizerHit> hits;

    auto minimizers = ExtractMinimizers(sequence, config_);

    for (const auto& m : minimizers) {
        auto it = impl_->index.find(m.hash);
        if (it == impl_->index.end()) continue;

        const auto& positions = it->second;

        // Skip high-occurrence minimizers
        if (positions.size() > config_.max_occ) {
            continue;
        }

        for (const auto& [ref_id, ref_pos] : positions) {
            MinimizerHit hit;
            hit.ref_id = ref_id;
            hit.ref_pos = ref_pos;
            hit.query_pos = m.pos;
            hit.same_strand = !m.is_reverse;
            hits.push_back(hit);

            if (hits.size() >= max_hits) {
                goto done;
            }
        }
    }
done:

    // Sort by (ref_id, ref_pos) for chaining
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
        return a.ref_pos < b.ref_pos;
    });

    return hits;
}

size_t MinimizerIndex::GetOccurrenceCount(uint64_t hash) const {
    auto it = impl_->index.find(hash);
    if (it == impl_->index.end()) return 0;
    return it->second.size();
}

const std::vector<IndexedSequence>& MinimizerIndex::GetSequences() const {
    return impl_->sequences;
}

const MinimizerStats& MinimizerIndex::GetStats() const {
    return impl_->stats;
}

bool MinimizerIndex::Empty() const {
    return impl_->index.empty();
}

size_t MinimizerIndex::Size() const {
    size_t count = 0;
    for (const auto& [hash, positions] : impl_->index) {
        count += positions.size();
    }
    return count;
}

}  // namespace llmap::classical
