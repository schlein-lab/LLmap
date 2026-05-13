// LLmap — Cluster assignment writer implementation.

#include "output/cluster_writer.h"

#include <chrono>
#include <fstream>
#include <set>
#include <sstream>

namespace llmap::output {

class ClusterWriterImpl {
public:
    std::ofstream file;
    std::string last_error;
    ClusterWriterStats stats;
    std::set<uint32_t> seen_clusters;
    bool header_written = false;
};

ClusterWriter::ClusterWriter(const std::filesystem::path& path,
                             const ClusterWriterConfig& config)
    : path_(path), config_(config), impl_(std::make_unique<ClusterWriterImpl>()) {}

ClusterWriter::~ClusterWriter() {
    Close();
}

ClusterWriter::ClusterWriter(ClusterWriter&&) noexcept = default;
ClusterWriter& ClusterWriter::operator=(ClusterWriter&&) noexcept = default;

std::unique_ptr<ClusterWriter> ClusterWriter::Create(
    const std::filesystem::path& path,
    const ClusterWriterConfig& config) {

    if (config.format == ClusterOutputFormat::Parquet) {
        // Parquet not yet supported
        return nullptr;
    }

    auto writer = std::unique_ptr<ClusterWriter>(new ClusterWriter(path, config));
    if (!writer->Initialize()) {
        return nullptr;
    }
    return writer;
}

bool ClusterWriter::Initialize() {
    impl_->file.open(path_, std::ios::out | std::ios::trunc);
    if (!impl_->file) {
        impl_->last_error = "Failed to open output file: " + path_.string();
        return false;
    }

    // Write header if requested
    if (config_.include_header) {
        impl_->file << "read_id" << config_.delimiter
                    << "cluster_id" << config_.delimiter
                    << "confidence" << config_.delimiter
                    << "is_representative" << config_.delimiter
                    << "cluster_size";
        if (config_.include_sequence_length) {
            impl_->file << config_.delimiter << "read_length";
        }
        impl_->file << "\n";
        impl_->header_written = true;
    }

    return true;
}

bool ClusterWriter::Write(const ClusterAssignment& assignment) {
    if (!impl_->file.is_open()) {
        impl_->last_error = "Writer not initialized or already closed";
        return false;
    }

    auto start = std::chrono::steady_clock::now();

    impl_->file << assignment.read_id << config_.delimiter
                << assignment.cluster_id << config_.delimiter
                << assignment.confidence << config_.delimiter
                << (assignment.is_representative ? "true" : "false") << config_.delimiter
                << assignment.cluster_size;

    if (config_.include_sequence_length) {
        impl_->file << config_.delimiter << assignment.read_length;
    }
    impl_->file << "\n";

    // Update stats
    impl_->stats.records_written++;
    if (impl_->seen_clusters.insert(assignment.cluster_id).second) {
        impl_->stats.clusters_written++;
    }
    if (assignment.is_representative) {
        impl_->stats.representatives_written++;
    }

    auto end = std::chrono::steady_clock::now();
    impl_->stats.write_time_ms += std::chrono::duration<float, std::milli>(end - start).count();

    return true;
}

bool ClusterWriter::WriteBatch(std::span<const ClusterAssignment> assignments) {
    for (const auto& assignment : assignments) {
        if (!Write(assignment)) {
            return false;
        }
    }
    return true;
}

bool ClusterWriter::Close() {
    if (impl_->file.is_open()) {
        impl_->file.close();
        return true;
    }
    return false;
}

ClusterWriterStats ClusterWriter::GetStats() const {
    return impl_->stats;
}

std::string ClusterWriter::LastError() const {
    return impl_->last_error;
}

// ========== Convenience functions ==========

bool WriteClusterAssignments(
    const std::filesystem::path& path,
    std::span<const ClusterAssignment> assignments,
    const ClusterWriterConfig& config) {

    auto writer = ClusterWriter::Create(path, config);
    if (!writer) {
        return false;
    }

    if (!writer->WriteBatch(assignments)) {
        return false;
    }

    return writer->Close();
}

}  // namespace llmap::output
