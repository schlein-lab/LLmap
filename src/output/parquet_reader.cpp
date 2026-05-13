// LLmap — Parquet probabilistic reader: core implementation.

#include "output/parquet_reader.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <unordered_map>

#ifdef LLMAP_HAS_ARROW
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#endif

namespace llmap::output {

namespace {

// Trim whitespace from both ends
std::string_view Trim(std::string_view s) {
    while (!s.empty() && std::isspace(s.front())) s.remove_prefix(1);
    while (!s.empty() && std::isspace(s.back())) s.remove_suffix(1);
    return s;
}

// Split a CSV line respecting basic quoting
std::vector<std::string> SplitCSV(std::string_view line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quote = false;

    for (char c : line) {
        if (c == '"') {
            in_quote = !in_quote;
        } else if (c == ',' && !in_quote) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

}  // namespace

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

std::vector<ProbabilityEntry> ParquetReader::ReadAll() {
    auto start = std::chrono::steady_clock::now();
    std::vector<ProbabilityEntry> entries;

#ifdef LLMAP_HAS_ARROW
    if (!impl_->is_csv_mode && impl_->table) {
        auto n = impl_->table->num_rows();

        auto read_ids = std::static_pointer_cast<arrow::StringArray>(
            impl_->table->column(0)->chunk(0));
        auto bucket_ids = std::static_pointer_cast<arrow::StringArray>(
            impl_->table->column(1)->chunk(0));
        auto probs = std::static_pointer_cast<arrow::FloatArray>(
            impl_->table->column(2)->chunk(0));
        auto confs = std::static_pointer_cast<arrow::FloatArray>(
            impl_->table->column(3)->chunk(0));
        auto levels = std::static_pointer_cast<arrow::UInt8Array>(
            impl_->table->column(4)->chunk(0));
        auto iters = std::static_pointer_cast<arrow::UInt32Array>(
            impl_->table->column(5)->chunk(0));
        auto collapsed = std::static_pointer_cast<arrow::BooleanArray>(
            impl_->table->column(6)->chunk(0));

        entries.reserve(n);
        for (int64_t i = 0; i < n; ++i) {
            ProbabilityEntry entry{
                .read_id = read_ids->GetString(i),
                .bucket_id = bucket_ids->GetString(i),
                .probability = probs->Value(i),
                .confidence = confs->Value(i),
                .level = levels->Value(i),
                .iteration = iters->Value(i),
                .is_collapsed = collapsed->Value(i),
            };

            if (entry.probability >= config_.min_probability) {
                if (!config_.skip_unmapped || entry.bucket_id != "*") {
                    impl_->UpdateStats(entry);
                    entries.push_back(std::move(entry));
                }
            }
        }
        impl_->at_eof = true;
        impl_->FinalizeStats();

        auto end = std::chrono::steady_clock::now();
        impl_->stats.read_time_ms =
            std::chrono::duration<float, std::milli>(end - start).count();
        return entries;
    }
#endif

    // CSV mode
    std::string line;
    while (std::getline(impl_->csv_file, line)) {
        auto entry_opt = ParseCSVLine(line);
        if (entry_opt) {
            auto& entry = *entry_opt;
            if (entry.probability >= config_.min_probability) {
                if (!config_.skip_unmapped || entry.bucket_id != "*") {
                    impl_->UpdateStats(entry);
                    entries.push_back(std::move(entry));
                }
            }
        }
    }
    impl_->at_eof = true;
    impl_->FinalizeStats();

    auto end = std::chrono::steady_clock::now();
    impl_->stats.read_time_ms =
        std::chrono::duration<float, std::milli>(end - start).count();

    return entries;
}

std::vector<ProbabilityEntry> ParquetReader::ReadBatch(std::size_t batch_size) {
    std::vector<ProbabilityEntry> entries;
    entries.reserve(batch_size);

#ifdef LLMAP_HAS_ARROW
    if (!impl_->is_csv_mode && impl_->table) {
        auto n = impl_->table->num_rows();
        if (impl_->current_row >= n) {
            impl_->at_eof = true;
            return entries;
        }

        auto read_ids = std::static_pointer_cast<arrow::StringArray>(
            impl_->table->column(0)->chunk(0));
        auto bucket_ids = std::static_pointer_cast<arrow::StringArray>(
            impl_->table->column(1)->chunk(0));
        auto probs = std::static_pointer_cast<arrow::FloatArray>(
            impl_->table->column(2)->chunk(0));
        auto confs = std::static_pointer_cast<arrow::FloatArray>(
            impl_->table->column(3)->chunk(0));
        auto levels = std::static_pointer_cast<arrow::UInt8Array>(
            impl_->table->column(4)->chunk(0));
        auto iters = std::static_pointer_cast<arrow::UInt32Array>(
            impl_->table->column(5)->chunk(0));
        auto collapsed = std::static_pointer_cast<arrow::BooleanArray>(
            impl_->table->column(6)->chunk(0));

        int64_t end_row = std::min(impl_->current_row + static_cast<int64_t>(batch_size), n);
        for (int64_t i = impl_->current_row; i < end_row; ++i) {
            ProbabilityEntry entry{
                .read_id = read_ids->GetString(i),
                .bucket_id = bucket_ids->GetString(i),
                .probability = probs->Value(i),
                .confidence = confs->Value(i),
                .level = levels->Value(i),
                .iteration = iters->Value(i),
                .is_collapsed = collapsed->Value(i),
            };

            if (entry.probability >= config_.min_probability) {
                if (!config_.skip_unmapped || entry.bucket_id != "*") {
                    impl_->UpdateStats(entry);
                    entries.push_back(std::move(entry));
                }
            }
        }
        impl_->current_row = end_row;
        if (impl_->current_row >= n) {
            impl_->at_eof = true;
            impl_->FinalizeStats();
        }
        return entries;
    }
#endif

    // CSV mode
    std::string line;
    while (entries.size() < batch_size && std::getline(impl_->csv_file, line)) {
        auto entry_opt = ParseCSVLine(line);
        if (entry_opt) {
            auto& entry = *entry_opt;
            if (entry.probability >= config_.min_probability) {
                if (!config_.skip_unmapped || entry.bucket_id != "*") {
                    impl_->UpdateStats(entry);
                    entries.push_back(std::move(entry));
                }
            }
        }
    }

    if (!impl_->csv_file) {
        impl_->at_eof = true;
        impl_->FinalizeStats();
    }

    return entries;
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

// ========== Free functions ==========

std::vector<ProbabilityEntry> ReadParquet(
    const std::filesystem::path& path,
    const ParquetReaderConfig& config) {

    auto reader = ParquetReader::Open(path, config);
    if (!reader) {
        return {};
    }
    return reader->ReadAll();
}

std::optional<ProbabilityEntry> ParseCSVLine(std::string_view line) {
    line = Trim(line);
    if (line.empty()) return std::nullopt;

    auto fields = SplitCSV(line);
    if (fields.size() < 7) return std::nullopt;

    try {
        ProbabilityEntry entry{
            .read_id = fields[0],
            .bucket_id = fields[1],
            .probability = std::stof(fields[2]),
            .confidence = std::stof(fields[3]),
            .level = static_cast<std::uint8_t>(std::stoul(fields[4])),
            .iteration = static_cast<std::uint32_t>(std::stoul(fields[5])),
            .is_collapsed = fields[6] == "1" || fields[6] == "true",
        };
        return entry;
    } catch (...) {
        return std::nullopt;
    }
}

std::vector<std::vector<ProbabilityEntry>> GroupByReadId(
    const std::vector<ProbabilityEntry>& entries) {

    std::unordered_map<std::string, std::vector<ProbabilityEntry>> groups;
    for (const auto& entry : entries) {
        groups[entry.read_id].push_back(entry);
    }

    std::vector<std::vector<ProbabilityEntry>> result;
    result.reserve(groups.size());
    for (auto& [read_id, group] : groups) {
        result.push_back(std::move(group));
    }
    return result;
}

RoundTripResult ValidateRoundTrip(
    const std::vector<ProbabilityEntry>& original,
    const std::vector<ProbabilityEntry>& reread,
    float probability_tolerance) {

    RoundTripResult result;
    result.entries_original = original.size();
    result.entries_reread = reread.size();

    if (original.size() != reread.size()) {
        result.error = "Entry count mismatch: " +
            std::to_string(original.size()) + " vs " +
            std::to_string(reread.size());
        return result;
    }

    // Sort both by (read_id, bucket_id) for comparison
    auto sorted_orig = original;
    auto sorted_reread = reread;

    auto cmp = [](const ProbabilityEntry& a, const ProbabilityEntry& b) {
        if (a.read_id != b.read_id) return a.read_id < b.read_id;
        return a.bucket_id < b.bucket_id;
    };

    std::sort(sorted_orig.begin(), sorted_orig.end(), cmp);
    std::sort(sorted_reread.begin(), sorted_reread.end(), cmp);

    for (std::size_t i = 0; i < sorted_orig.size(); ++i) {
        const auto& orig = sorted_orig[i];
        const auto& re = sorted_reread[i];

        bool match = (orig.read_id == re.read_id) &&
                     (orig.bucket_id == re.bucket_id) &&
                     (std::fabs(orig.probability - re.probability) <= probability_tolerance) &&
                     (std::fabs(orig.confidence - re.confidence) <= probability_tolerance) &&
                     (orig.level == re.level) &&
                     (orig.iteration == re.iteration) &&
                     (orig.is_collapsed == re.is_collapsed);

        if (!match) {
            result.mismatches++;
        }
    }

    result.success = (result.mismatches == 0);
    if (!result.success && result.error.empty()) {
        result.error = std::to_string(result.mismatches) + " entry mismatches";
    }

    return result;
}

}  // namespace llmap::output
