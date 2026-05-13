// LLmap — Parquet reader: ReadAll and ReadBatch implementations.

#include "output/parquet_reader.h"
#include "output/parquet_reader_impl.h"

#include <algorithm>
#include <chrono>

namespace llmap::output {

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

        int64_t end_row = std::min(
            impl_->current_row + static_cast<int64_t>(batch_size), n);
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

}  // namespace llmap::output
