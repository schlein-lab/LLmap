// LLmap — Memory-mapped FASTA reader core implementation.
//
// Transparently handles both plain .fa/.fasta and gzipped .fa.gz inputs.
// Plain inputs use mmap (lazy paging from disk, zero-copy reads).
// Gzipped inputs are decompressed once at open-time into a heap buffer
// of the same memory layout; the rest of the reader operates on a
// `const char*` view and is oblivious to the source.

#include "io/mmap_fasta.h"
#include "io/mmap_fasta_internal.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

namespace llmap::io {

size_t GetPageSize() {
    static const size_t page_size = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    return page_size;
}

MmapFastaImpl::~MmapFastaImpl() {
    if (owned_heap) {
        delete[] heap_buffer;
        heap_buffer = nullptr;
        mapped = nullptr;
    } else if (mapped != nullptr && mapped != MAP_FAILED) {
        munmap(mapped, file_size);
        mapped = nullptr;
    }
    if (fd >= 0) {
        close(fd);
        fd = -1;
    }
}

const char* MmapFastaImpl::Data() const {
    return static_cast<const char*>(mapped);
}

namespace {

// True when path ends with ".gz" (case-insensitive).
bool HasGzExtension(const std::filesystem::path& p) {
    auto ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext == ".gz" || ext == ".bgz";
}

// Decompress an entire gz file into a freshly-allocated heap buffer.
// Returns nullptr on failure with last_error populated.
char* SlurpGzipToHeap(const std::filesystem::path& path,
                     std::size_t& out_size,
                     std::string& last_error) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) {
        last_error = "gzopen failed: " + path.string();
        return nullptr;
    }
    // Larger internal buffer reduces overhead on multi-GB inputs.
    gzbuffer(gz, 1u << 20);  // 1 MiB

    constexpr std::size_t initial_cap = 1ULL << 28;  // 256 MiB
    std::size_t cap  = initial_cap;
    std::size_t size = 0;
    char* buf = new (std::nothrow) char[cap];
    if (!buf) {
        last_error = "heap allocation failed";
        gzclose(gz);
        return nullptr;
    }

    constexpr std::size_t chunk = 1u << 20;  // 1 MiB
    for (;;) {
        if (size + chunk > cap) {
            std::size_t new_cap = cap * 2;
            char* nb = new (std::nothrow) char[new_cap];
            if (!nb) {
                last_error = "heap grow failed at " + std::to_string(new_cap) + " bytes";
                delete[] buf;
                gzclose(gz);
                return nullptr;
            }
            std::copy(buf, buf + size, nb);
            delete[] buf;
            buf = nb;
            cap = new_cap;
        }
        int n = gzread(gz, buf + size, static_cast<unsigned>(chunk));
        if (n < 0) {
            int errnum = 0;
            const char* msg = gzerror(gz, &errnum);
            last_error = std::string("gzread error: ") + (msg ? msg : "(nul)");
            delete[] buf;
            gzclose(gz);
            return nullptr;
        }
        if (n == 0) break;
        size += static_cast<std::size_t>(n);
    }
    gzclose(gz);
    out_size = size;
    return buf;
}

}  // namespace

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

    // ── gz path: decompress once into a heap buffer, then proceed
    //    through the same indexing logic as the mmap path.
    if (HasGzExtension(path_)) {
        std::size_t decompressed = 0;
        char* buf = SlurpGzipToHeap(path_, decompressed, impl_->last_error);
        if (!buf) {
            return;  // last_error set
        }
        impl_->owned_heap  = true;
        impl_->heap_buffer = buf;
        impl_->mapped      = buf;
        impl_->file_size   = decompressed;
        impl_->fd          = -1;  // no fd for the heap path

        if (config_.build_index_on_open) {
            if (!impl_->BuildIndex()) {
                delete[] impl_->heap_buffer;
                impl_->heap_buffer = nullptr;
                impl_->mapped      = nullptr;
                impl_->owned_heap  = false;
            }
        }
        return;
    }

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
    if (ext == ".gz" || ext == ".bgz") {
        // .fa.gz / .fasta.gz: check the secondary extension too.
        auto stem = path.stem();
        auto stem_ext = stem.extension().string();
        std::transform(stem_ext.begin(), stem_ext.end(), stem_ext.begin(), ::tolower);
        if (stem_ext != ".fa" && stem_ext != ".fasta"
            && stem_ext != ".fna" && stem_ext != ".fas") {
            return false;
        }
        // gz path → trust the extension; opening to peek is expensive.
        return true;
    }
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
