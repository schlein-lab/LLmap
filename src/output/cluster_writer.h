// LLmap — Cluster assignment writer.
//
// Outputs Stage 1 (self-interference) cluster assignments. Supports TSV format
// in V1.0; Parquet support added in Phase 6 when Arrow is integrated.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace llmap::output {

// A read with its cluster assignment (output record)
struct ClusterAssignment {
    std::string read_id;
    uint32_t cluster_id;
    float confidence;
    bool is_representative;
    size_t cluster_size;

    // Original sequence length
    size_t read_length;
};

// Statistics from writing cluster assignments
struct ClusterWriterStats {
    size_t records_written = 0;
    size_t clusters_written = 0;
    size_t representatives_written = 0;
    float write_time_ms = 0.0f;
};

// Output format
enum class ClusterOutputFormat {
    TSV,      // Tab-separated values (always available)
    Parquet,  // Apache Parquet (requires Arrow integration, Phase 6+)
};

// Configuration for the cluster writer
struct ClusterWriterConfig {
    ClusterOutputFormat format = ClusterOutputFormat::TSV;
    bool include_header = true;
    bool include_sequence_length = true;
    std::string delimiter = "\t";  // For TSV format
};

// Forward declaration
class ClusterWriterImpl;

// Cluster assignment writer
class ClusterWriter {
public:
    // Factory: create a writer for the specified output path
    static std::unique_ptr<ClusterWriter> Create(
        const std::filesystem::path& path,
        const ClusterWriterConfig& config = {});

    ~ClusterWriter();

    // Non-copyable, movable
    ClusterWriter(const ClusterWriter&) = delete;
    ClusterWriter& operator=(const ClusterWriter&) = delete;
    ClusterWriter(ClusterWriter&&) noexcept;
    ClusterWriter& operator=(ClusterWriter&&) noexcept;

    // Write a single assignment
    bool Write(const ClusterAssignment& assignment);

    // Write a batch of assignments
    bool WriteBatch(std::span<const ClusterAssignment> assignments);

    // Finalize and close the file
    bool Close();

    // Get accumulated statistics
    ClusterWriterStats GetStats() const;

    // Get last error
    std::string LastError() const;

private:
    explicit ClusterWriter(const std::filesystem::path& path,
                           const ClusterWriterConfig& config);
    bool Initialize();

    std::filesystem::path path_;
    ClusterWriterConfig config_;
    std::unique_ptr<ClusterWriterImpl> impl_;
};

// Convenience: write all assignments to a file in one call
bool WriteClusterAssignments(
    const std::filesystem::path& path,
    std::span<const ClusterAssignment> assignments,
    const ClusterWriterConfig& config = {});

}  // namespace llmap::output
