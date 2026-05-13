// LLmap — Memory-mapped FASTA reader sequence access methods.

#include "io/mmap_fasta.h"
#include "io/mmap_fasta_internal.h"

#include <cctype>
#include <sys/mman.h>

namespace llmap::io {

std::string MmapFastaImpl::ExtractSequence(size_t data_offset, size_t data_end) const {
    if (mapped == nullptr) return {};

    const char* data = Data();
    std::string result;
    result.reserve(data_end - data_offset);

    for (size_t i = data_offset; i < data_end; ++i) {
        char c = data[i];
        if (!std::isspace(static_cast<unsigned char>(c))) {
            result.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
        }
    }

    return result;
}

std::string MmapFastaImpl::ExtractSubsequence(size_t data_offset, size_t data_end,
                                               size_t start, size_t len) const {
    if (mapped == nullptr) return {};

    const char* data = Data();
    std::string result;
    result.reserve(len);

    size_t base_idx = 0;
    for (size_t i = data_offset; i < data_end && result.size() < len; ++i) {
        char c = data[i];
        if (!std::isspace(static_cast<unsigned char>(c))) {
            if (base_idx >= start) {
                result.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
            }
            ++base_idx;
        }
    }

    return result;
}

std::optional<MmapSequence> MmapFastaReader::GetSequence(size_t index) const {
    if (!impl_ || index >= impl_->sequences.size()) {
        return std::nullopt;
    }
    const auto& entry = impl_->sequences[index];
    return MmapSequence{entry.name, entry.data_offset, entry.length};
}

std::optional<MmapSequence> MmapFastaReader::GetSequence(std::string_view name) const {
    if (!impl_) return std::nullopt;

    auto it = impl_->name_to_index.find(std::string(name));
    if (it == impl_->name_to_index.end()) {
        return std::nullopt;
    }
    return GetSequence(it->second);
}

std::string_view MmapFastaReader::GetSequenceRaw(size_t index) const {
    if (!impl_ || index >= impl_->sequences.size()) {
        return {};
    }
    const auto& entry = impl_->sequences[index];
    return std::string_view(impl_->Data() + entry.data_offset,
                            entry.data_end - entry.data_offset);
}

std::string_view MmapFastaReader::GetSequenceRaw(std::string_view name) const {
    if (!impl_) return {};

    auto it = impl_->name_to_index.find(std::string(name));
    if (it == impl_->name_to_index.end()) {
        return {};
    }
    return GetSequenceRaw(it->second);
}

std::string MmapFastaReader::GetSequenceData(size_t index) const {
    if (!impl_ || index >= impl_->sequences.size()) {
        return {};
    }
    const auto& entry = impl_->sequences[index];
    return impl_->ExtractSequence(entry.data_offset, entry.data_end);
}

std::string MmapFastaReader::GetSequenceData(std::string_view name) const {
    if (!impl_) return {};

    auto it = impl_->name_to_index.find(std::string(name));
    if (it == impl_->name_to_index.end()) {
        return {};
    }
    return GetSequenceData(it->second);
}

std::string MmapFastaReader::GetSubsequence(size_t index, size_t start, size_t len) const {
    if (!impl_ || index >= impl_->sequences.size()) {
        return {};
    }
    const auto& entry = impl_->sequences[index];
    return impl_->ExtractSubsequence(entry.data_offset, entry.data_end, start, len);
}

std::string MmapFastaReader::GetSubsequence(std::string_view name, size_t start, size_t len) const {
    if (!impl_) return {};

    auto it = impl_->name_to_index.find(std::string(name));
    if (it == impl_->name_to_index.end()) {
        return {};
    }
    return GetSubsequence(it->second, start, len);
}

MmapStats MmapFastaReader::GetStats() const {
    MmapStats stats;
    if (!impl_ || impl_->mapped == nullptr) {
        return stats;
    }

    stats.file_size = impl_->file_size;
    stats.num_sequences = impl_->sequences.size();

    for (const auto& entry : impl_->sequences) {
        stats.total_bases += entry.length;
    }

#ifdef __linux__
    size_t page_size = GetPageSize();
    size_t num_pages = (impl_->file_size + page_size - 1) / page_size;

    std::vector<unsigned char> vec(num_pages);
    if (mincore(impl_->mapped, impl_->file_size, vec.data()) == 0) {
        for (unsigned char v : vec) {
            if (v & 1) {
                ++stats.resident_pages;
            }
        }
    }
    stats.resident_fraction = num_pages > 0 ?
        static_cast<double>(stats.resident_pages) / num_pages : 0.0;
#endif

    return stats;
}

void MmapFastaReader::AdviseSequential() {
    if (impl_ && impl_->mapped != nullptr) {
        madvise(impl_->mapped, impl_->file_size, MADV_SEQUENTIAL);
    }
}

void MmapFastaReader::AdviseRandom() {
    if (impl_ && impl_->mapped != nullptr) {
        madvise(impl_->mapped, impl_->file_size, MADV_RANDOM);
    }
}

void MmapFastaReader::AdviseWillNeed(size_t index) {
    if (!impl_ || impl_->mapped == nullptr || index >= impl_->sequences.size()) {
        return;
    }
    const auto& entry = impl_->sequences[index];
    size_t page_size = GetPageSize();
    size_t aligned_offset = (entry.data_offset / page_size) * page_size;
    size_t length = entry.data_end - aligned_offset;
    madvise(static_cast<char*>(impl_->mapped) + aligned_offset, length, MADV_WILLNEED);
}

void MmapFastaReader::AdviseDontNeed(size_t index) {
    if (!impl_ || impl_->mapped == nullptr || index >= impl_->sequences.size()) {
        return;
    }
    const auto& entry = impl_->sequences[index];
    size_t page_size = GetPageSize();
    size_t aligned_offset = (entry.data_offset / page_size) * page_size;
    size_t length = entry.data_end - aligned_offset;
    madvise(static_cast<char*>(impl_->mapped) + aligned_offset, length, MADV_DONTNEED);
}

}  // namespace llmap::io
