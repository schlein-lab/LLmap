// LLmap — Parquet probabilistic output writer: core implementation.

#include "output/parquet_writer.h"
#include "output/parquet_writer_impl.h"

#include <iomanip>

namespace llmap::output {

// ========== Impl methods ==========

bool ParquetWriterImpl::FlushBuffer(const ParquetWriterConfig& config) {
    if (buffer.empty()) return true;

#ifdef LLMAP_HAS_ARROW
    if (!is_csv_mode) {
        return WriteBufferToParquet(config);
    }
#endif
    return WriteBufferToCSV();
}

bool ParquetWriterImpl::WriteBufferToCSV() {
    if (!csv_file.is_open()) {
        last_error = "CSV file not open";
        return false;
    }

    for (const auto& entry : buffer) {
        csv_file << entry.read_id << ","
                 << entry.bucket_id << ","
                 << std::fixed << std::setprecision(6) << entry.probability << ","
                 << std::fixed << std::setprecision(6) << entry.confidence << ","
                 << static_cast<int>(entry.level) << ","
                 << entry.iteration << ","
                 << (entry.is_collapsed ? "1" : "0") << "\n";

        stats.entries_written++;
    }

    buffer.clear();
    return true;
}

#ifdef LLMAP_HAS_ARROW
bool ParquetWriterImpl::WriteBufferToParquet(const ParquetWriterConfig& config) {
    if (!parquet_writer) {
        last_error = "Parquet writer not initialized";
        return false;
    }

    auto n = static_cast<int64_t>(buffer.size());

    // Build arrays
    arrow::StringBuilder read_id_builder;
    arrow::StringBuilder bucket_id_builder;
    arrow::FloatBuilder prob_builder;
    arrow::FloatBuilder conf_builder;
    arrow::UInt8Builder level_builder;
    arrow::UInt32Builder iter_builder;
    arrow::BooleanBuilder collapsed_builder;

    ARROW_RETURN_NOT_OK(read_id_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(bucket_id_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(prob_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(conf_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(level_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(iter_builder.Reserve(n));
    ARROW_RETURN_NOT_OK(collapsed_builder.Reserve(n));

    for (const auto& entry : buffer) {
        ARROW_RETURN_NOT_OK(read_id_builder.Append(entry.read_id));
        ARROW_RETURN_NOT_OK(bucket_id_builder.Append(entry.bucket_id));
        ARROW_RETURN_NOT_OK(prob_builder.Append(entry.probability));
        ARROW_RETURN_NOT_OK(conf_builder.Append(entry.confidence));
        ARROW_RETURN_NOT_OK(level_builder.Append(entry.level));
        ARROW_RETURN_NOT_OK(iter_builder.Append(entry.iteration));
        ARROW_RETURN_NOT_OK(collapsed_builder.Append(entry.is_collapsed));
    }

    std::shared_ptr<arrow::Array> read_ids, bucket_ids, probs, confs, levels, iters, collapsed;
    ARROW_RETURN_NOT_OK(read_id_builder.Finish(&read_ids));
    ARROW_RETURN_NOT_OK(bucket_id_builder.Finish(&bucket_ids));
    ARROW_RETURN_NOT_OK(prob_builder.Finish(&probs));
    ARROW_RETURN_NOT_OK(conf_builder.Finish(&confs));
    ARROW_RETURN_NOT_OK(level_builder.Finish(&levels));
    ARROW_RETURN_NOT_OK(iter_builder.Finish(&iters));
    ARROW_RETURN_NOT_OK(collapsed_builder.Finish(&collapsed));

    auto table = arrow::Table::Make(
        schema,
        {read_ids, bucket_ids, probs, confs, levels, iters, collapsed});

    auto status = parquet_writer->WriteTable(*table, n);
    if (!status.ok()) {
        last_error = "Failed to write Parquet table: " + status.ToString();
        return false;
    }

    stats.entries_written += buffer.size();
    buffer.clear();
    return true;
}
#endif

// ========== ParquetWriter ==========

ParquetWriter::ParquetWriter(const std::filesystem::path& path,
                             const ParquetWriterConfig& config)
    : path_(path), config_(config),
      impl_(std::make_unique<ParquetWriterImpl>()) {}

ParquetWriter::~ParquetWriter() {
    Close();
}

ParquetWriter::ParquetWriter(ParquetWriter&&) noexcept = default;
ParquetWriter& ParquetWriter::operator=(ParquetWriter&&) noexcept = default;

std::unique_ptr<ParquetWriter> ParquetWriter::Create(
    const std::filesystem::path& path,
    const ParquetWriterConfig& config) {

    auto writer = std::unique_ptr<ParquetWriter>(
        new ParquetWriter(path, config));
    if (!writer->Initialize()) {
        return nullptr;
    }
    return writer;
}

bool ParquetWriter::Initialize() {
#ifdef LLMAP_HAS_ARROW
    // Determine output format
    bool use_arrow = ArrowAvailable() &&
                     config_.format == ParquetOutputFormat::Parquet;

    if (use_arrow) {
        impl_->is_csv_mode = false;

        // Build schema
        impl_->schema = arrow::schema({
            arrow::field("read_id", arrow::utf8()),
            arrow::field("bucket_id", arrow::utf8()),
            arrow::field("probability", arrow::float32()),
            arrow::field("confidence", arrow::float32()),
            arrow::field("level", arrow::uint8()),
            arrow::field("iteration", arrow::uint32()),
            arrow::field("is_collapsed", arrow::boolean())
        });

        // Open output file
        auto result = arrow::io::FileOutputStream::Open(path_.string());
        if (!result.ok()) {
            impl_->last_error = "Failed to open Parquet file: " + result.status().ToString();
            return false;
        }
        impl_->arrow_file = *result;

        // Create Parquet writer
        parquet::WriterProperties::Builder props_builder;
        if (config_.compress) {
            if (config_.compression == "snappy") {
                props_builder.compression(parquet::Compression::SNAPPY);
            } else if (config_.compression == "gzip") {
                props_builder.compression(parquet::Compression::GZIP);
            } else if (config_.compression == "zstd") {
                props_builder.compression(parquet::Compression::ZSTD);
            } else {
                props_builder.compression(parquet::Compression::UNCOMPRESSED);
            }
        } else {
            props_builder.compression(parquet::Compression::UNCOMPRESSED);
        }
        auto props = props_builder.build();

        auto arrow_props = parquet::ArrowWriterProperties::Builder()
            .store_schema()
            ->build();

        auto status = parquet::arrow::FileWriter::Open(
            *impl_->schema,
            arrow::default_memory_pool(),
            impl_->arrow_file,
            props,
            arrow_props,
            &impl_->parquet_writer);

        if (!status.ok()) {
            impl_->last_error = "Failed to create Parquet writer: " + status.ToString();
            return false;
        }

        return true;
    }
#endif

    // CSV fallback
    impl_->is_csv_mode = true;

    // For CSV, change extension if path ends with .parquet
    auto csv_path = path_;
    if (path_.extension() == ".parquet") {
        csv_path.replace_extension(".csv");
    }

    impl_->csv_file.open(csv_path, std::ios::out | std::ios::trunc);
    if (!impl_->csv_file) {
        impl_->last_error = "Failed to open CSV file: " + csv_path.string();
        return false;
    }

    // Write CSV header
    impl_->csv_file << kCsvHeader << "\n";
    impl_->header_written = true;

    return true;
}

}  // namespace llmap::output
