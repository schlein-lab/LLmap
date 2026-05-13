// LLmap — Memory-mapped FASTA reader core implementation.

#include "io/mmap_fasta.h"
#include "io/mmap_fasta_internal.h"

#include <algorithm>
#include <cctype>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace llmap::io {

size_t GetPageSize() {
    static const size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    return page_size;
}

MmapFastaImpl::~MmapFastaImpl() {
    if (mapped != nullptr && mapped != MAP_FAILED) {
        munmap(mapped, file_size);
    }
    if (fd >= 0) {
        close(fd);
    }
}

const char* MmapFastaImpl::Data() const {
    return static_cast<const char*>(mapped);
}

bool MmapFastaImpl::BuildIndex() {
    if (mapped == nullptr || mapped == MAP_FAILED || file_size == 0) {
        return false;
    }

    const char* data = Data();
    size_t pos = 0;

    while (pos < file_size) {
        while (pos < file_size && std::isspace(static_cast<unsigned char>(data[pos]))) {
            ++pos;
        }
        if (pos >= file_size) break;

        if (data[pos] != '>') {
            last_error = "Expected '>' at position " + std::to_string(pos);
            return false;
        }

        SequenceEntry entry;
        entry.name_offset = pos;
        ++pos;

        size_t name_start = pos;
        while (pos < file_size && !std::isspace(static_cast<unsigned char>(data[pos]))) {
            ++pos;
        }
        entry.name = std::string(data + name_start, pos - name_start);

        while (pos < file_size && data[pos] != '\n') {
            ++pos;
        }
        if (pos < file_size) ++pos;

        entry.data_offset = pos;

        size_t seq_length = 0;
        while (pos < file_size && data[pos] != '>') {
            if (!std::isspace(static_cast<unsigned char>(data[pos]))) {
                ++seq_length;
            }
            ++pos;
        }

        entry.data_end = pos;
        entry.length = seq_length;

        size_t idx = sequences.size();
        name_to_index[entry.name] = idx;
        sequences.push_back(std::move(entry));
    }

    return true;
}

MmapFastaReader::MmapFastaReader(const std::filesystem::path& path,
                                   const MmapFastaConfig& config)
    : path_(path), config_(config), impl_(std::make_unique<MmapFastaImpl>()) {

    impl_->fd = open(path_.c_str(), O_RDONLY);
    if (impl_->fd < 0) {
        impl_->last_error = "Failed to open file: " + path_.string();
        return;
    }

    struct stat st;
    if (fstat(impl_->fd, &st) < 0) {
        impl_->last_error = "Failed to stat file: " + path_.string();
        close(impl_->fd);
        impl_->fd = -1;
        return;
    }

    impl_->file_size = static_cast<size_t>(st.st_size);
    if (impl_->file_size == 0) {
        impl_->last_error = "File is empty: " + path_.string();
        close(impl_->fd);
        impl_->fd = -1;
        return;
    }

    int flags = MAP_PRIVATE;
    if (config_.prefault_pages) {
        flags |= MAP_POPULATE;
    }

    impl_->mapped = mmap(nullptr, impl_->file_size, PROT_READ, flags, impl_->fd, 0);
    if (impl_->mapped == MAP_FAILED) {
        impl_->last_error = "Failed to mmap file: " + path_.string();
        impl_->mapped = nullptr;
        close(impl_->fd);
        impl_->fd = -1;
        return;
    }

    if (config_.lock_pages) {
        mlock(impl_->mapped, impl_->file_size);
    }

    if (config_.read_ahead_bytes > 0) {
        madvise(impl_->mapped, std::min(config_.read_ahead_bytes, impl_->file_size),
                MADV_WILLNEED);
    }

    if (config_.build_index_on_open) {
        if (!impl_->BuildIndex()) {
            munmap(impl_->mapped, impl_->file_size);
            impl_->mapped = nullptr;
            close(impl_->fd);
            impl_->fd = -1;
        }
    }
}

MmapFastaReader::~MmapFastaReader() = default;
MmapFastaReader::MmapFastaReader(MmapFastaReader&&) noexcept = default;
MmapFastaReader& MmapFastaReader::operator=(MmapFastaReader&&) noexcept = default;

bool MmapFastaReader::IsValid() const {
    return impl_ && impl_->mapped != nullptr;
}

std::string MmapFastaReader::LastError() const {
    return impl_ ? impl_->last_error : "Not initialized";
}

size_t MmapFastaReader::NumSequences() const {
    return impl_ ? impl_->sequences.size() : 0;
}

std::vector<std::string_view> MmapFastaReader::SequenceNames() const {
    std::vector<std::string_view> names;
    if (impl_) {
        names.reserve(impl_->sequences.size());
        for (const auto& entry : impl_->sequences) {
            names.push_back(entry.name);
        }
    }
    return names;
}

bool IsFastaFile(const std::filesystem::path& path) {
    auto ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext != ".fa" && ext != ".fasta" && ext != ".fna" && ext != ".fas") {
        return false;
    }

    std::ifstream f(path);
    char c;
    if (f >> c) {
        return c == '>';
    }
    return false;
}

}  // namespace llmap::io
