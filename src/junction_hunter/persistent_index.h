// LLmap — junction_hunter: persistent k-mer index for the NAHR-pair panel.
//
// Disk format for a single-k index file. Built once via
// `llmap junction-index build` and mmap'd by `llmap junction-hunt
// --index-dir`. The format is intentionally simple: a 256-byte header
// followed by a flat array of entries sorted by (hash, pair_id).
//
// Lookup is a binary search on the sorted array; multiple entries
// per hash live in a contiguous range (returned via Lookup()).
//
// One file per cascade tier:
//   <dir>/tier_k{N}.bin     for each k in CascadeConfig::k_values
//   <dir>/panel.sha256      provenance stamp (matches against runtime panel)
//   <dir>/reference.sha256  provenance stamp (matches against runtime ref)
//   <dir>/MANIFEST.tsv      one row per tier: k, n_entries, byte_size

#pragma once

#include "junction_hunter/pair_kmer_index.h"   // LocusClass

#include <cstdint>
#include <string>
#include <utility>

namespace llmap::junction_hunter {

inline constexpr char kPersistentIndexMagic[16] = {
    'L','L','M','A','P','_','N','A','H','R','_','I','D','X','\0','\0'
};
inline constexpr std::uint32_t kPersistentIndexVersion = 1;

#pragma pack(push, 1)
struct PersistentIndexHeader {
    char         magic[16];          ///< "LLMAP_NAHR_IDX\0\0"
    std::uint32_t version;           ///< == kPersistentIndexVersion
    std::uint8_t  k;                 ///< k for this tier
    std::uint8_t  reserved8[3];      ///< padding to 4-byte boundary
    std::uint64_t n_entries;         ///< number of IndexEntry records that follow
    std::uint64_t n_pairs;           ///< panel size at build time
    std::uint8_t  panel_sha256[32];  ///< sha256 of the pair-panel TSV
    std::uint8_t  ref_sha256[32];    ///< sha256 of the reference FASTA
    char          built_iso[32];     ///< ISO 8601 build timestamp
    char          built_host[32];    ///< hostname at build time
    std::uint8_t  pad[88];           ///< total header size = 256 bytes
};
static_assert(sizeof(PersistentIndexHeader) == 256,
              "PersistentIndexHeader must be exactly 256 bytes");

/// Cache entry. Sorted by (hash, pair_id) on disk. 17 bytes packed.
struct PersistentIndexEntry {
    std::uint64_t hash;       ///< canonical k-mer hash (FNV-1a)
    std::uint32_t pair_id;    ///< 0-based index into the pair-panel TSV
    std::uint32_t offset;     ///< 0-based offset within the cls region
    std::uint8_t  cls;        ///< LocusClass cast to uint8 (LcrUp/LcrDown/Interior/Ambiguous)
};
static_assert(sizeof(PersistentIndexEntry) == 17,
              "PersistentIndexEntry must be exactly 17 bytes (packed)");
#pragma pack(pop)

/// mmap-backed reader for a single tier_k{N}.bin file.
class PersistentKmerIndex {
public:
    PersistentKmerIndex() = default;
    ~PersistentKmerIndex();

    PersistentKmerIndex(const PersistentKmerIndex&) = delete;
    PersistentKmerIndex& operator=(const PersistentKmerIndex&) = delete;
    PersistentKmerIndex(PersistentKmerIndex&&) noexcept;
    PersistentKmerIndex& operator=(PersistentKmerIndex&&) noexcept;

    /// Memory-map an existing tier file. Verifies magic + version.
    /// Returns false on I/O error or format mismatch.
    bool Open(const std::string& path);

    /// Release the mapping.
    void Close();

    bool IsOpen() const noexcept { return fd_ >= 0; }
    const PersistentIndexHeader& Header() const noexcept { return *header_; }
    const PersistentIndexEntry*  Entries() const noexcept { return entries_; }
    std::uint64_t                NumEntries() const noexcept { return header_ ? header_->n_entries : 0; }

    /// Look up all entries with this canonical hash. Returns a contiguous
    /// half-open range [first, last). Empty range if not found.
    std::pair<const PersistentIndexEntry*, const PersistentIndexEntry*>
    Lookup(std::uint64_t hash) const noexcept;

    /// Advise the OS that access pattern is roughly random.
    void HintRandomAccess() noexcept;

    /// Last I/O error message after a failed Open.
    const std::string& LastError() const noexcept { return last_error_; }

private:
    int                          fd_{-1};
    void*                        mmap_base_{nullptr};
    std::size_t                  mmap_size_{0};
    const PersistentIndexHeader* header_{nullptr};
    const PersistentIndexEntry*  entries_{nullptr};
    std::string                  last_error_;
};

}  // namespace llmap::junction_hunter
