// LLmap — Parquet reader internal implementation header.
// Shared between split parquet_reader_*.cpp files.

#pragma once

#include "output/parquet_reader.h"

#include <fstream>
#include <set>
#include <string>

#ifdef LLMAP_HAS_ARROW
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#endif

namespace llmap::output {

// Internal implementation state
class ParquetReaderImpl {
public:
    std::ifstream csv_file;
    std::string last_error;
    ParquetReaderStats stats;
    bool is_csv_mode = true;
    bool header_read = false;
    bool at_eof = false;
    ParquetOutputFormat detected_format = ParquetOutputFormat::CSV;

#ifdef LLMAP_HAS_ARROW
    std::unique_ptr<parquet::arrow::FileReader> parquet_reader;
    std::shared_ptr<arrow::Table> table;
    int64_t current_row = 0;
    int64_t total_rows = 0;
#endif

    std::set<std::string> unique_reads;
    std::set<std::string> unique_buckets;

    void UpdateStats(const ProbabilityEntry& entry) {
        stats.entries_read++;
        unique_reads.insert(entry.read_id);
        unique_buckets.insert(entry.bucket_id);
        if (entry.is_collapsed) stats.collapsed_entries++;
        if (entry.bucket_id == "*") stats.unmapped_entries++;
    }

    void FinalizeStats() {
        stats.unique_reads = unique_reads.size();
        stats.unique_buckets = unique_buckets.size();
    }
};

// CSV parsing utilities (defined in parquet_reader_util.cpp)
namespace detail {
std::string_view Trim(std::string_view s);
std::vector<std::string> SplitCSV(std::string_view line);
}  // namespace detail

}  // namespace llmap::output
