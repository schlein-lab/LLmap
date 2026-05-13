// LLmap — FASTA reader implementation.

#include "io/fasta_reader.h"

#include <algorithm>

namespace llmap::io {

namespace {

void TrimRight(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' ||
                          s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
}

}  // namespace

class FastaReaderImpl {
public:
    std::ifstream file;
    bool at_eof = false;
    std::string last_error;
    size_t records_read = 0;
    std::string pending_header;  // Header line read ahead
    std::string line_buffer;
};

FastaReader::FastaReader(const std::filesystem::path& path,
                         const FastaReaderConfig& config)
    : path_(path), config_(config), impl_(std::make_unique<FastaReaderImpl>()) {

    impl_->file.open(path_, std::ios::in);
    if (!impl_->file) {
        impl_->last_error = "Failed to open file: " + path_.string();
        impl_->at_eof = true;
    }
}

FastaReader::~FastaReader() = default;

FastaReader::FastaReader(FastaReader&&) noexcept = default;
FastaReader& FastaReader::operator=(FastaReader&&) noexcept = default;

FastaRecord FastaReader::Next() {
    if (impl_->at_eof) {
        return {};
    }

    if (config_.max_records > 0 && impl_->records_read >= config_.max_records) {
        impl_->at_eof = true;
        return {};
    }

    FastaRecord record;

    // Read header line (either pending or next line)
    std::string header;
    if (!impl_->pending_header.empty()) {
        header = std::move(impl_->pending_header);
        impl_->pending_header.clear();
    } else {
        // Skip any leading empty lines
        while (std::getline(impl_->file, header)) {
            TrimRight(header);
            if (!header.empty()) break;
        }
        if (header.empty()) {
            impl_->at_eof = true;
            return {};
        }
    }

    if (header.empty() || header[0] != '>') {
        impl_->last_error = "Expected '>' at start of FASTA record";
        impl_->at_eof = true;
        return {};
    }

    // Extract name (first word after >)
    size_t name_end = header.find_first_of(" \t", 1);
    if (name_end == std::string::npos) {
        record.name = header.substr(1);
    } else {
        record.name = header.substr(1, name_end - 1);
    }

    // Read sequence lines until next header or EOF
    while (std::getline(impl_->file, impl_->line_buffer)) {
        TrimRight(impl_->line_buffer);

        if (impl_->line_buffer.empty()) {
            continue;  // Skip empty lines within sequence
        }

        if (impl_->line_buffer[0] == '>') {
            // Next record header - save it for later
            impl_->pending_header = impl_->line_buffer;
            break;
        }

        record.sequence += impl_->line_buffer;
    }

    // Check for EOF
    if (impl_->pending_header.empty() && !impl_->file) {
        impl_->at_eof = true;
    }

    // Process sequence
    if (config_.uppercase_sequence) {
        std::transform(record.sequence.begin(), record.sequence.end(),
                       record.sequence.begin(), ::toupper);
    }

    // Skip N-only sequences if configured
    if (config_.skip_N_only && record.IsValid()) {
        bool all_n = std::all_of(record.sequence.begin(), record.sequence.end(),
                                 [](char c) { return c == 'N' || c == 'n'; });
        if (all_n) {
            return Next();  // Skip and try next
        }
    }

    if (record.IsValid()) {
        ++impl_->records_read;
    }

    return record;
}

bool FastaReader::HasMore() const {
    if (impl_->at_eof) {
        return false;
    }
    if (!impl_->pending_header.empty()) {
        return true;
    }
    return impl_->file.peek() != std::ifstream::traits_type::eof();
}

size_t FastaReader::RecordsRead() const {
    return impl_->records_read;
}

std::string FastaReader::LastError() const {
    return impl_->last_error;
}

std::vector<FastaRecord> ReadFasta(const std::filesystem::path& path) {
    std::vector<FastaRecord> records;
    FastaReader reader(path);

    while (reader.HasMore()) {
        auto record = reader.Next();
        if (record.IsValid()) {
            records.push_back(std::move(record));
        }
    }

    return records;
}

}  // namespace llmap::io
