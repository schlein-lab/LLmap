// LLmap — FASTQ reader implementation.

#include "io/fastq_reader.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <sstream>

namespace llmap::io {

namespace {

constexpr size_t kDefaultBufferSize = 1 << 20;  // 1 MB

// Trim whitespace from end of string in place
void TrimRight(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' ||
                          s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
}

}  // namespace

class FastqReaderImpl {
public:
    std::ifstream file;
    bool at_eof = false;
    std::string last_error;
    FastqReaderStats stats;
    bool is_compressed = false;

    std::string line_buffer;
};

FastqReader::FastqReader(const std::filesystem::path& path,
                         const FastqReaderConfig& config)
    : path_(path), config_(config), impl_(std::make_unique<FastqReaderImpl>()) {}

FastqReader::~FastqReader() = default;

FastqReader::FastqReader(FastqReader&&) noexcept = default;
FastqReader& FastqReader::operator=(FastqReader&&) noexcept = default;

std::unique_ptr<FastqReader> FastqReader::Open(
    const std::filesystem::path& path,
    const FastqReaderConfig& config) {

    if (!std::filesystem::exists(path)) {
        return nullptr;
    }

    auto reader = std::unique_ptr<FastqReader>(new FastqReader(path, config));
    if (!reader->Initialize()) {
        return nullptr;
    }
    return reader;
}

bool FastqReader::Initialize() {
    impl_->is_compressed = IsGzipFile(path_);

    if (impl_->is_compressed) {
        impl_->last_error = "Gzip compression not yet supported in V1.0";
        return false;
    }

    impl_->file.open(path_, std::ios::binary);
    if (!impl_->file) {
        impl_->last_error = "Failed to open file: " + path_.string();
        return false;
    }

    return true;
}

std::optional<FastqRecord> FastqReader::Next() {
    if (impl_->at_eof) {
        return std::nullopt;
    }

    if (config_.max_records > 0 && impl_->stats.total_records >= config_.max_records) {
        impl_->at_eof = true;
        return std::nullopt;
    }

    FastqRecord record;

    // Line 1: @ + read ID
    if (!std::getline(impl_->file, impl_->line_buffer)) {
        impl_->at_eof = true;
        return std::nullopt;
    }
    TrimRight(impl_->line_buffer);

    if (impl_->line_buffer.empty() || impl_->line_buffer[0] != '@') {
        impl_->last_error = "Expected '@' at start of FASTQ record";
        impl_->stats.invalid_records++;
        return std::nullopt;
    }
    record.id = impl_->line_buffer.substr(1);

    // Line 2: sequence
    if (!std::getline(impl_->file, record.sequence)) {
        impl_->last_error = "Unexpected end of file reading sequence";
        impl_->at_eof = true;
        impl_->stats.invalid_records++;
        return std::nullopt;
    }
    TrimRight(record.sequence);

    if (config_.uppercase_sequence) {
        std::transform(record.sequence.begin(), record.sequence.end(),
                       record.sequence.begin(), ::toupper);
    }

    // Line 3: + separator
    if (!std::getline(impl_->file, impl_->line_buffer)) {
        impl_->last_error = "Unexpected end of file reading separator";
        impl_->at_eof = true;
        impl_->stats.invalid_records++;
        return std::nullopt;
    }
    TrimRight(impl_->line_buffer);

    if (impl_->line_buffer.empty() || impl_->line_buffer[0] != '+') {
        impl_->last_error = "Expected '+' separator line";
        impl_->stats.invalid_records++;
        return std::nullopt;
    }

    // Line 4: quality
    if (!std::getline(impl_->file, record.quality)) {
        impl_->last_error = "Unexpected end of file reading quality";
        impl_->at_eof = true;
        impl_->stats.invalid_records++;
        return std::nullopt;
    }
    TrimRight(record.quality);

    // Validate length match
    if (config_.validate_quality && record.quality.size() != record.sequence.size()) {
        impl_->last_error = "Quality string length mismatch";
        impl_->stats.invalid_records++;
        return std::nullopt;
    }

    // Apply length filters
    const size_t len = record.sequence.size();
    if (config_.min_length > 0 && len < config_.min_length) {
        return Next();  // Skip and try next
    }
    if (config_.max_length > 0 && len > config_.max_length) {
        return Next();  // Skip and try next
    }

    // Update stats
    impl_->stats.total_records++;
    impl_->stats.total_bases += len;
    if (impl_->stats.total_records == 1) {
        impl_->stats.min_length = len;
        impl_->stats.max_length = len;
    } else {
        impl_->stats.min_length = std::min(impl_->stats.min_length, len);
        impl_->stats.max_length = std::max(impl_->stats.max_length, len);
    }
    impl_->stats.avg_length = static_cast<float>(impl_->stats.total_bases) /
                               static_cast<float>(impl_->stats.total_records);

    return record;
}

std::vector<FastqRecord> FastqReader::NextBatch(size_t batch_size) {
    std::vector<FastqRecord> batch;
    batch.reserve(batch_size);

    for (size_t i = 0; i < batch_size; ++i) {
        auto record = Next();
        if (!record) break;
        batch.push_back(std::move(*record));
    }

    return batch;
}

std::vector<FastqRecord> FastqReader::ReadAll() {
    std::vector<FastqRecord> records;

    while (auto record = Next()) {
        records.push_back(std::move(*record));
    }

    return records;
}

bool FastqReader::HasMore() const {
    if (impl_->at_eof) {
        return false;
    }
    // Peek to check if there's more data
    if (impl_->file.peek() == std::ifstream::traits_type::eof()) {
        return false;
    }
    return true;
}

bool FastqReader::IsCompressed() const {
    return impl_->is_compressed;
}

FastqReaderStats FastqReader::GetStats() const {
    return impl_->stats;
}

std::string FastqReader::LastError() const {
    return impl_->last_error;
}

// ========== Convenience functions ==========

std::vector<FastqRecord> ReadFastq(
    const std::filesystem::path& path,
    size_t max_records) {

    FastqReaderConfig config;
    config.max_records = max_records;

    auto reader = FastqReader::Open(path, config);
    if (!reader) {
        return {};
    }

    return reader->ReadAll();
}

size_t CountFastqRecords(const std::filesystem::path& path) {
    auto reader = FastqReader::Open(path);
    if (!reader) {
        return 0;
    }

    size_t count = 0;
    while (reader->Next()) {
        ++count;
    }

    return count;
}

bool IsGzipFile(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;

    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);

    // Gzip magic: 0x1f 0x8b
    return file.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

bool ValidateFastqRecord(const FastqRecord& record) {
    if (record.id.empty()) return false;
    if (record.sequence.empty()) return false;
    if (record.quality.size() != record.sequence.size()) return false;

    // Validate sequence characters
    for (char c : record.sequence) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            return false;
        }
    }

    return true;
}

}  // namespace llmap::io
