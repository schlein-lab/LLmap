// LLmap — BAM/SAM output writer.
//
// Outputs alignment results in SAM/BAM format for drop-in compatibility
// with samtools/bcftools/IGV. Handles CIGAR generation from alignment hits.
//
// Uses htslib if available (LLMAP_HAS_HTSLIB), otherwise writes plain SAM.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "core/alignment_record.h"

namespace llmap::output {

// SAM/BAM output format
enum class BamOutputFormat {
    SAM,  // Plain text (always available)
    BAM,  // Binary (requires htslib)
};

// Reference sequence for header generation
struct ReferenceSequence {
    std::string name;
    std::uint64_t length;
};

// Statistics from writing alignments
struct BamWriterStats {
    std::size_t records_written = 0;
    std::size_t mapped_written = 0;
    std::size_t unmapped_written = 0;
    std::size_t tentative_written = 0;
    std::size_t alternatives_written = 0;
    float write_time_ms = 0.0f;
};

// Configuration for the BAM writer
struct BamWriterConfig {
    BamOutputFormat format = BamOutputFormat::SAM;

    // Header fields
    std::string sample_name = "sample";
    std::string read_group_id = "1";
    std::string platform = "PACBIO";
    std::string program_name = "llmap";
    std::string program_version = "1.0.0";

    // Output options
    bool include_alternatives = true;     // XA tag for alternative alignments
    bool include_tentative = true;        // Write tentative as unmapped with tags
    bool include_wavecollapse_tags = true; // XC:cluster, XL:level, XI:iteration
    bool include_paralog_tags = true;      // XP:paralog probability
    bool write_unmapped = true;            // Include unmapped reads

    // SAM-specific
    bool write_header = true;

    // Sequences (for CIGAR validation)
    bool include_sequences = false;  // If true, query sequences are written
};

// Forward declaration
class BamWriterImpl;

// BAM/SAM writer for AlignmentRecords
class BamWriter {
public:
    // Factory: create a writer for the specified output path
    static std::unique_ptr<BamWriter> Create(
        const std::filesystem::path& path,
        std::span<const ReferenceSequence> references,
        const BamWriterConfig& config = {});

    ~BamWriter();

    // Non-copyable, movable
    BamWriter(const BamWriter&) = delete;
    BamWriter& operator=(const BamWriter&) = delete;
    BamWriter(BamWriter&&) noexcept;
    BamWriter& operator=(BamWriter&&) noexcept;

    // Write a single alignment record
    // query_seq is optional; only needed if config.include_sequences=true
    bool Write(const AlignmentRecord& record,
               std::string_view query_seq = "");

    // Write a batch of alignment records
    bool WriteBatch(std::span<const AlignmentRecord> records);

    // Finalize and close the file
    bool Close();

    // Get accumulated statistics
    BamWriterStats GetStats() const;

    // Get last error
    std::string LastError() const;

    // Check if htslib is available for BAM output
    static bool HtslibAvailable();

private:
    explicit BamWriter(const std::filesystem::path& path,
                       std::span<const ReferenceSequence> references,
                       const BamWriterConfig& config);
    bool Initialize();

    std::filesystem::path path_;
    BamWriterConfig config_;
    std::vector<ReferenceSequence> references_;
    std::unique_ptr<BamWriterImpl> impl_;
};

// Convenience: write all records to a file in one call
bool WriteAlignments(
    const std::filesystem::path& path,
    std::span<const AlignmentRecord> records,
    std::span<const ReferenceSequence> references,
    const BamWriterConfig& config = {});

// CIGAR manipulation utilities
namespace cigar {

// Generate CIGAR string from alignment positions and edit distance
// Simple heuristic when full traceback is unavailable
std::string GenerateSimpleCigar(std::uint64_t query_len,
                                 std::uint64_t ref_start,
                                 std::uint64_t ref_end,
                                 std::uint32_t nm);

// Parse CIGAR ops string into length and reference consumed
std::pair<std::uint64_t, std::uint64_t> CigarStats(std::string_view cigar);

// Validate CIGAR string format
bool ValidateCigar(std::string_view cigar);

}  // namespace cigar

}  // namespace llmap::output
