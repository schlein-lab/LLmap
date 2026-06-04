// LLmap — junction_hunter: persistent k-mer index (mmap reader).

#include "junction_hunter/persistent_index.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace llmap::junction_hunter {

PersistentKmerIndex::~PersistentKmerIndex() { Close(); }

PersistentKmerIndex::PersistentKmerIndex(PersistentKmerIndex&& o) noexcept
    : fd_(o.fd_), mmap_base_(o.mmap_base_), mmap_size_(o.mmap_size_),
      header_(o.header_), entries_(o.entries_),
      last_error_(std::move(o.last_error_)) {
    o.fd_ = -1;
    o.mmap_base_ = nullptr;
    o.mmap_size_ = 0;
    o.header_ = nullptr;
    o.entries_ = nullptr;
}

PersistentKmerIndex& PersistentKmerIndex::operator=(PersistentKmerIndex&& o) noexcept {
    if (this != &o) {
        Close();
        fd_ = o.fd_;
        mmap_base_ = o.mmap_base_;
        mmap_size_ = o.mmap_size_;
        header_ = o.header_;
        entries_ = o.entries_;
        last_error_ = std::move(o.last_error_);
        o.fd_ = -1;
        o.mmap_base_ = nullptr;
        o.mmap_size_ = 0;
        o.header_ = nullptr;
        o.entries_ = nullptr;
    }
    return *this;
}

void PersistentKmerIndex::Close() {
    if (mmap_base_ != nullptr) {
        ::munmap(mmap_base_, mmap_size_);
        mmap_base_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
    mmap_size_ = 0;
    header_ = nullptr;
    entries_ = nullptr;
}

bool PersistentKmerIndex::Open(const std::string& path) {
    Close();
    last_error_.clear();

    fd_ = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) {
        last_error_ = std::string("open: ") + std::strerror(errno);
        return false;
    }
    struct stat st;
    if (::fstat(fd_, &st) != 0) {
        last_error_ = std::string("fstat: ") + std::strerror(errno);
        Close();
        return false;
    }
    mmap_size_ = static_cast<std::size_t>(st.st_size);
    if (mmap_size_ < sizeof(PersistentIndexHeader)) {
        last_error_ = "file smaller than header";
        Close();
        return false;
    }
    mmap_base_ = ::mmap(nullptr, mmap_size_, PROT_READ, MAP_SHARED, fd_, 0);
    if (mmap_base_ == MAP_FAILED) {
        last_error_ = std::string("mmap: ") + std::strerror(errno);
        mmap_base_ = nullptr;
        Close();
        return false;
    }

    header_ = static_cast<const PersistentIndexHeader*>(mmap_base_);
    if (std::memcmp(header_->magic, kPersistentIndexMagic, 16) != 0) {
        last_error_ = "bad magic";
        Close();
        return false;
    }
    if (header_->version != kPersistentIndexVersion) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "version mismatch: file=%u expected=%u",
                      header_->version, kPersistentIndexVersion);
        last_error_ = buf;
        Close();
        return false;
    }
    const std::size_t expected = sizeof(PersistentIndexHeader)
                               + header_->n_entries * sizeof(PersistentIndexEntry);
    if (expected > mmap_size_) {
        last_error_ = "file truncated (entry count exceeds file size)";
        Close();
        return false;
    }

    entries_ = reinterpret_cast<const PersistentIndexEntry*>(
        static_cast<const std::uint8_t*>(mmap_base_) + sizeof(PersistentIndexHeader));
    return true;
}

void PersistentKmerIndex::HintRandomAccess() noexcept {
    if (mmap_base_ && mmap_size_ > 0) {
        ::madvise(mmap_base_, mmap_size_, MADV_RANDOM);
    }
}

std::pair<const PersistentIndexEntry*, const PersistentIndexEntry*>
PersistentKmerIndex::Lookup(std::uint64_t hash) const noexcept {
    if (!entries_ || header_->n_entries == 0) {
        return {nullptr, nullptr};
    }
    const PersistentIndexEntry* begin = entries_;
    const PersistentIndexEntry* end   = entries_ + header_->n_entries;
    auto lo = std::lower_bound(begin, end, hash,
        [](const PersistentIndexEntry& e, std::uint64_t h) {
            return e.hash < h;
        });
    if (lo == end || lo->hash != hash) {
        return {lo, lo};
    }
    auto hi = std::upper_bound(lo, end, hash,
        [](std::uint64_t h, const PersistentIndexEntry& e) {
            return h < e.hash;
        });
    return {lo, hi};
}

}  // namespace llmap::junction_hunter
