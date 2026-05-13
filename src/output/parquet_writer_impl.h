// LLmap — Parquet writer internal implementation header.
// Shared across parquet_writer*.cpp files.

#pragma once

#include "output/parquet_writer.h"

#include <fstream>
#include <vector>

#ifdef LLMAP_HAS_ARROW
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

namespace llmap::output {

// CSV header for fallback format
inline constexpr std::string_view kCsvHeader =
    "read_id,bucket_id,probability,confidence,level,iteration,is_collapsed";

class ParquetWriterImpl {
public:
    std::vector<ProbabilityEntry> buffer;
    std::ofstream csv_file;
    std::string last_error;
    ParquetWriterStats stats;
    bool is_csv_mode = false;
    bool header_written = false;
    bool closed = false;

#ifdef LLMAP_HAS_ARROW
    std::shared_ptr<arrow::io::FileOutputStream> arrow_file;
    std::unique_ptr<parquet::arrow::FileWriter> parquet_writer;
    std::shared_ptr<arrow::Schema> schema;
#endif

    bool FlushBuffer(const ParquetWriterConfig& config);
    bool WriteBufferToCSV();
#ifdef LLMAP_HAS_ARROW
    bool WriteBufferToParquet(const ParquetWriterConfig& config);
#endif
};

}  // namespace llmap::output
