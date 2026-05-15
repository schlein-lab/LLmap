#include "annot/variant_priors.h"

#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace llmap::annot {

namespace {

constexpr char kMagic[8] = {'L', 'L', 'M', 'V', 'A', 'R', '0', '1'};
constexpr size_t kHeaderBytes = 4096;

struct OnDiskContigEntry {
    char name[40];
    uint32_t contig_length_bp;
    uint32_t bucket_count;
    uint64_t records_offset;
    char pad[8];
};
static_assert(sizeof(OnDiskContigEntry) == 64,
              "OnDiskContigEntry must be 64 bytes");

struct OnDiskHeader {
    char magic[8];
    uint32_t version;
    uint32_t bucket_bp;
    uint32_t num_contigs;
    uint32_t reserved0;
};
static_assert(sizeof(OnDiskHeader) == 24, "header prefix must be 24 bytes");

}  // namespace

// --- Builder --------------------------------------------------------------

VariantPriorBuilder::VariantPriorBuilder(uint32_t bucket_bp)
    : bucket_bp_(bucket_bp) {}

void VariantPriorBuilder::AddContig(std::string name, uint32_t length_bp) {
    VariantPriorContigInfo info;
    info.name = std::move(name);
    info.contig_length_bp = length_bp;
    info.bucket_count = (length_bp + bucket_bp_ - 1) / bucket_bp_;
    contigs_.push_back(std::move(info));
    records_.emplace_back(contigs_.back().bucket_count);
}

VariantPriorRecord& VariantPriorBuilder::Bucket(uint32_t contig_idx,
                                                uint32_t bucket_idx) {
    return records_[contig_idx][bucket_idx];
}

bool VariantPriorBuilder::Save(const std::filesystem::path& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;

    // header
    char header_block[kHeaderBytes];
    std::memset(header_block, 0, kHeaderBytes);
    OnDiskHeader hdr;
    std::memcpy(hdr.magic, kMagic, 8);
    hdr.version = 1;
    hdr.bucket_bp = bucket_bp_;
    hdr.num_contigs = static_cast<uint32_t>(contigs_.size());
    hdr.reserved0 = 0;
    std::memcpy(header_block, &hdr, sizeof(hdr));

    // contig entries follow the header prefix
    uint64_t cur_offset = kHeaderBytes;
    std::vector<OnDiskContigEntry> entries(contigs_.size());
    for (size_t i = 0; i < contigs_.size(); ++i) {
        std::memset(&entries[i], 0, sizeof(OnDiskContigEntry));
        std::strncpy(entries[i].name, contigs_[i].name.c_str(),
                     sizeof(entries[i].name) - 1);
        entries[i].contig_length_bp = contigs_[i].contig_length_bp;
        entries[i].bucket_count = contigs_[i].bucket_count;
        entries[i].records_offset = cur_offset;
        cur_offset += static_cast<uint64_t>(contigs_[i].bucket_count) *
                      sizeof(VariantPriorRecord);
    }

    size_t entries_size = entries.size() * sizeof(OnDiskContigEntry);
    if (sizeof(OnDiskHeader) + entries_size > kHeaderBytes) {
        std::cerr << "[variant_priors] too many contigs for fixed 4k header\n";
        return false;
    }
    std::memcpy(header_block + sizeof(OnDiskHeader), entries.data(),
                entries_size);

    f.write(header_block, kHeaderBytes);
    if (!f) return false;

    // record arrays in order
    for (size_t i = 0; i < records_.size(); ++i) {
        f.write(reinterpret_cast<const char*>(records_[i].data()),
                static_cast<std::streamsize>(records_[i].size() *
                                             sizeof(VariantPriorRecord)));
        if (!f) return false;
    }
    return true;
}

// --- Reader ---------------------------------------------------------------

std::unique_ptr<VariantPriorStore> VariantPriorStore::Open(
    const std::filesystem::path& path) {

    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        std::cerr << "[variant_priors] open failed: " << path << "\n";
        return nullptr;
    }
    struct stat st;
    if (::fstat(fd, &st) < 0 || st.st_size < static_cast<off_t>(kHeaderBytes)) {
        ::close(fd);
        return nullptr;
    }
    void* base = ::mmap(nullptr, static_cast<size_t>(st.st_size),
                        PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) {
        ::close(fd);
        return nullptr;
    }
    if (std::memcmp(base, kMagic, 8) != 0) {
        ::munmap(base, st.st_size);
        ::close(fd);
        std::cerr << "[variant_priors] bad magic in " << path << "\n";
        return nullptr;
    }

    auto store = std::unique_ptr<VariantPriorStore>(new VariantPriorStore());
    store->fd_ = fd;
    store->mmap_base_ = base;
    store->mmap_size_ = static_cast<size_t>(st.st_size);

    const char* p = static_cast<const char*>(base);
    OnDiskHeader hdr;
    std::memcpy(&hdr, p, sizeof(hdr));
    store->bucket_bp_ = hdr.bucket_bp;

    const OnDiskContigEntry* entries =
        reinterpret_cast<const OnDiskContigEntry*>(p + sizeof(OnDiskHeader));
    for (uint32_t i = 0; i < hdr.num_contigs; ++i) {
        VariantPriorContigInfo c;
        c.name = std::string(entries[i].name, ::strnlen(entries[i].name,
                                                       sizeof(entries[i].name)));
        c.contig_length_bp = entries[i].contig_length_bp;
        c.bucket_count = entries[i].bucket_count;
        c.records_offset = entries[i].records_offset;
        store->contigs_.push_back(c);
        store->contig_name_to_idx_[c.name] = i;
        store->contig_records_.push_back(
            reinterpret_cast<const VariantPriorRecord*>(p + c.records_offset));
    }

    return store;
}

VariantPriorStore::~VariantPriorStore() {
    if (mmap_base_ && mmap_base_ != MAP_FAILED) {
        ::munmap(mmap_base_, mmap_size_);
    }
    if (fd_ >= 0) ::close(fd_);
}

const VariantPriorRecord* VariantPriorStore::LookupBucket(
    uint32_t contig_idx, uint32_t bucket_idx) const {

    if (contig_idx >= contigs_.size()) return nullptr;
    if (bucket_idx >= contigs_[contig_idx].bucket_count) return nullptr;
    return &contig_records_[contig_idx][bucket_idx];
}

const VariantPriorRecord* VariantPriorStore::LookupPosition(
    std::string_view contig_name, uint32_t pos_bp) const {

    auto it = contig_name_to_idx_.find(std::string(contig_name));
    if (it == contig_name_to_idx_.end()) return nullptr;
    uint32_t bucket_idx = pos_bp / bucket_bp_;
    return LookupBucket(it->second, bucket_idx);
}

}  // namespace llmap::annot
