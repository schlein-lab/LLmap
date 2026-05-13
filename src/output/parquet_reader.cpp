// LLmap — Parquet probabilistic reader: core implementation.

#include "output/parquet_reader.h"
#include "output/parquet_reader_impl.h"

#include <filesystem>

namespace llmap::output {

ParquetReader::ParquetReader(const std::filesystem::path& path,
                             const ParquetReaderConfig& config)
    : path_(path), config_(config),
      impl_(std::make_unique<ParquetReaderImpl>()) {}

ParquetReader::~ParquetReader() = default;

ParquetReader::ParquetReader(ParquetReader&&) noexcept = default;
ParquetReader& ParquetReader::operator=(ParquetReader&&) noexcept = default;

std::unique_ptr<ParquetReader> ParquetReader::Open(
    const std::filesystem::path& path,
    const ParquetReaderConfig& config) {

    auto reader = std::unique_ptr<ParquetReader>(
        new ParquetReader(path, config));
    if (!reader->Initialize()) {
        return nullptr;
    }
    return reader;
}

bool ParquetReader::Initialize() {
    if (!std::filesystem::exists(path_)) {
        impl_->last_error = "File not found: " + path_.string();
        return false;
    }

#ifdef LLMAP_HAS_ARROW
    // Try to open as Parquet first
    if (path_.extension() == ".parquet") {
        auto result = arrow::io::ReadableFile::Open(path_.string());
        if (result.ok()) {
            auto arrow_file = *result;
            parquet::arrow::FileReaderBuilder builder;
            auto status = builder.Open(arrow_file);
            if (status.ok()) {
                status = builder.Build(&impl_->parquet_reader);
                if (status.ok()) {
                    impl_->is_csv_mode = false;
                    impl_->detected_format = ParquetOutputFormat::Parquet;

                    // Read entire table for simplicity
                    status = impl_->parquet_reader->ReadTable(&impl_->table);
                    if (status.ok()) {
                        impl_->total_rows = impl_->table->num_rows();
                        return true;
                    }
                }
            }
        }
        // Fall through to CSV if Parquet reading fails
    }
#endif

    // Open as CSV
    auto csv_path = path_;
    if (path_.extension() == ".parquet") {
        csv_path.replace_extension(".csv");
        if (!std::filesystem::exists(csv_path)) {
            impl_->last_error = "Neither Parquet nor CSV file found";
            return false;
        }
    }

    impl_->csv_file.open(csv_path);
    if (!impl_->csv_file) {
        impl_->last_error = "Failed to open CSV file: " + csv_path.string();
        return false;
    }

    impl_->is_csv_mode = true;
    impl_->detected_format = ParquetOutputFormat::CSV;

    // Skip header line
    std::string header;
    if (!std::getline(impl_->csv_file, header)) {
        impl_->last_error = "CSV file is empty";
        return false;
    }
    impl_->header_read = true;

    return true;
}

bool ParquetReader::HasMore() const {
    return !impl_->at_eof;
}

ParquetReaderStats ParquetReader::GetStats() const {
    return impl_->stats;
}

std::string ParquetReader::LastError() const {
    return impl_->last_error;
}

ParquetOutputFormat ParquetReader::DetectedFormat() const {
    return impl_->detected_format;
}

}  // namespace llmap::output
