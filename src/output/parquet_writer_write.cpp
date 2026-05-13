// LLmap — Parquet writer: Write/Batch/Close methods.

#include "output/parquet_writer.h"
#include "output/parquet_writer_impl.h"

#include <chrono>

namespace llmap::output {

bool ParquetWriter::Write(const AlignmentRecord& record) {
    if (impl_->closed) {
        impl_->last_error = "Writer already closed";
        return false;
    }

    auto start = std::chrono::steady_clock::now();

    // Convert record to probability entries
    auto entries = RecordToEntries(record, config_.min_probability);

    // Filter based on config
    if (record.status == AlignmentStatus::Unmapped && !config_.include_unmapped) {
        impl_->stats.records_written++;
        impl_->stats.unmapped_records++;
        return true;
    }
    if (record.status == AlignmentStatus::Tentative && !config_.include_tentative) {
        impl_->stats.records_written++;
        impl_->stats.tentative_records++;
        return true;
    }

    // Add to buffer
    for (auto& entry : entries) {
        impl_->buffer.push_back(std::move(entry));
    }

    // Update stats
    impl_->stats.records_written++;
    switch (record.status) {
        case AlignmentStatus::Mapped:
            impl_->stats.mapped_records++;
            break;
        case AlignmentStatus::Tentative:
            impl_->stats.tentative_records++;
            break;
        case AlignmentStatus::Unmapped:
            impl_->stats.unmapped_records++;
            break;
    }

    // Flush if buffer is large enough
    if (impl_->buffer.size() >= config_.row_group_size) {
        if (!impl_->FlushBuffer(config_)) {
            return false;
        }
    }

    auto end = std::chrono::steady_clock::now();
    impl_->stats.write_time_ms +=
        std::chrono::duration<float, std::milli>(end - start).count();

    return true;
}

bool ParquetWriter::WriteEntries(std::span<const ProbabilityEntry> entries) {
    auto start = std::chrono::steady_clock::now();

    for (const auto& entry : entries) {
        if (entry.probability >= config_.min_probability) {
            impl_->buffer.push_back(entry);
        }
    }

    if (impl_->buffer.size() >= config_.row_group_size) {
        if (!impl_->FlushBuffer(config_)) {
            return false;
        }
    }

    auto end = std::chrono::steady_clock::now();
    impl_->stats.write_time_ms +=
        std::chrono::duration<float, std::milli>(end - start).count();

    return true;
}

bool ParquetWriter::WriteBatch(std::span<const AlignmentRecord> records) {
    for (const auto& record : records) {
        if (!Write(record)) {
            return false;
        }
    }
    return true;
}

bool ParquetWriter::Close() {
    if (impl_->closed) {
        return false;
    }

    // Flush remaining buffer
    if (!impl_->buffer.empty()) {
        if (!impl_->FlushBuffer(config_)) {
            impl_->closed = true;
            return false;
        }
    }

    impl_->closed = true;

#ifdef LLMAP_HAS_ARROW
    if (!impl_->is_csv_mode && impl_->parquet_writer) {
        auto status = impl_->parquet_writer->Close();
        if (!status.ok()) {
            impl_->last_error = "Failed to close Parquet writer: " + status.ToString();
            return false;
        }
        impl_->parquet_writer.reset();

        if (impl_->arrow_file) {
            auto close_status = impl_->arrow_file->Close();
            if (!close_status.ok()) {
                impl_->last_error = "Failed to close Parquet file: " + close_status.ToString();
                return false;
            }
        }
        return true;
    }
#endif

    if (impl_->csv_file.is_open()) {
        impl_->csv_file.close();
        return true;
    }

    return false;
}

ParquetWriterStats ParquetWriter::GetStats() const {
    return impl_->stats;
}

std::string ParquetWriter::LastError() const {
    return impl_->last_error;
}

bool ParquetWriter::ArrowAvailable() {
#ifdef LLMAP_HAS_ARROW
    return true;
#else
    return false;
#endif
}

// ========== Convenience functions ==========

bool WriteParquet(
    const std::filesystem::path& path,
    std::span<const AlignmentRecord> records,
    const ParquetWriterConfig& config) {

    auto writer = ParquetWriter::Create(path, config);
    if (!writer) {
        return false;
    }

    if (!writer->WriteBatch(records)) {
        return false;
    }

    return writer->Close();
}

}  // namespace llmap::output
