// LLmap -- Layer 4 variant-prior store, memory-mapped binary.
//
// Layout (little-endian, fixed-stride records sharded by chromosome):
//
//   header (4096 bytes, page-aligned)
//   per-chromosome records, sorted by bucket index
//
// Header:
//   char     magic[8];           = "LLMVAR01"
//   uint32_t version;            = 1
//   uint32_t bucket_bp;          = 100
//   uint32_t num_contigs;
//   uint32_t reserved0;
//   contig_entry[num_contigs]:
//     char     name[40];          null-padded
//     uint32_t contig_length_bp;
//     uint32_t bucket_count;
//     uint64_t records_offset;   offset within file to this contig's records
//
// Record (24 bytes per bucket):
//   uint32_t n_snv;              # SNV alt alleles seen in DB
//   float    snv_max_af;         max AF across SNVs in bucket
//   uint32_t n_sv;               # SVs overlapping bucket
//   float    sv_max_af;          max AF across SVs in bucket
//   uint16_t sv_type_bits;       bitmask of SV types: DEL=1, DUP=2, INS=4,
//                                INV=8, BND=16, MEI=32, CNV=64
//   uint16_t flags;              bit 0: clinvar_pathogenic_in_bucket
//   uint32_t reserved;
//
// Total: 4096 + 24 * sum(bucket_count). For human GRCh38 with bucket_bp=100,
// ~3.1e9 / 100 = 31 M buckets -> ~720 MB. Mmap-lookup is O(1) by
// contig + bucket index, sub-microsecond.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::annot {

struct VariantPriorRecord {
    uint32_t n_snv = 0;
    float snv_max_af = 0.0f;
    uint32_t n_sv = 0;
    float sv_max_af = 0.0f;
    uint16_t sv_type_bits = 0;
    uint16_t flags = 0;
    uint32_t reserved = 0;
};

static_assert(sizeof(VariantPriorRecord) == 24,
              "VariantPriorRecord on-disk layout must be 24 bytes");

enum VariantSvType : uint16_t {
    kSvDel = 1u << 0,
    kSvDup = 1u << 1,
    kSvIns = 1u << 2,
    kSvInv = 1u << 3,
    kSvBnd = 1u << 4,
    kSvMei = 1u << 5,
    kSvCnv = 1u << 6,
};

enum VariantFlag : uint16_t {
    kClinVarPathogenic = 1u << 0,
};

struct VariantPriorContigInfo {
    std::string name;
    uint32_t contig_length_bp = 0;
    uint32_t bucket_count = 0;
    uint64_t records_offset = 0;  // byte offset in file
};

// Open a .priors file via mmap and provide O(1) bucket lookup.
class VariantPriorStore {
public:
    static std::unique_ptr<VariantPriorStore> Open(
        const std::filesystem::path& path);

    ~VariantPriorStore();
    VariantPriorStore(const VariantPriorStore&) = delete;
    VariantPriorStore& operator=(const VariantPriorStore&) = delete;

    uint32_t BucketBp() const { return bucket_bp_; }
    const std::vector<VariantPriorContigInfo>& Contigs() const { return contigs_; }

    // O(1) by contig_idx + bucket. Returns nullptr if out of range.
    const VariantPriorRecord* LookupBucket(uint32_t contig_idx,
                                           uint32_t bucket_idx) const;

    // Convenience: lookup by chromosome name + bp position.
    const VariantPriorRecord* LookupPosition(std::string_view contig_name,
                                             uint32_t pos_bp) const;

private:
    VariantPriorStore() = default;

    int fd_ = -1;
    void* mmap_base_ = nullptr;
    size_t mmap_size_ = 0;

    uint32_t bucket_bp_ = 100;
    std::vector<VariantPriorContigInfo> contigs_;
    std::unordered_map<std::string, uint32_t> contig_name_to_idx_;
    // Pointer to the start of each contig's record array (inside the mmap).
    std::vector<const VariantPriorRecord*> contig_records_;
};

// Build a .priors file from an in-memory bucket array. The caller fills
// records[contig_idx][bucket] for every contig+bucket, then calls Save().
class VariantPriorBuilder {
public:
    VariantPriorBuilder(uint32_t bucket_bp = 100);
    void AddContig(std::string name, uint32_t length_bp);
    VariantPriorRecord& Bucket(uint32_t contig_idx, uint32_t bucket_idx);
    bool Save(const std::filesystem::path& path) const;

    uint32_t BucketBp() const { return bucket_bp_; }
    size_t NumContigs() const { return contigs_.size(); }

private:
    uint32_t bucket_bp_;
    std::vector<VariantPriorContigInfo> contigs_;
    std::vector<std::vector<VariantPriorRecord>> records_;
};

}  // namespace llmap::annot
