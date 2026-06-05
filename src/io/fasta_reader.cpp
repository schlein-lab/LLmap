// LLmap — FASTA reader implementation.
//
// Transparently reads plain .fa/.fasta or gzipped .fa.gz inputs via zlib.

#include "io/fasta_reader.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <memory>
#include <zlib.h>

namespace llmap::io {

namespace {

void TrimRight(std::string& s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' ||
                          s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
}

bool HasGzExtension(const std::filesystem::path& p) {
    auto ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext == ".gz" || ext == ".bgz";
}

}  // namespace

// Thin wrapper that hides whether we read from ifstream or gzFile.
class FastaInput {
public:
    virtual ~FastaInput() = default;
    virtual bool GetLine(std::string& out) = 0;
    virtual bool Eof() const = 0;
    virtual std::string LastError() const { return {}; }
};

namespace {

class IfstreamInput : public FastaInput {
public:
    explicit IfstreamInput(const std::filesystem::path& p) {
        file.open(p, std::ios::in);
        if (!file) err = "Failed to open: " + p.string();
    }
    bool GetLine(std::string& out) override {
        if (!file) return false;
        if (!std::getline(file, out)) return false;
        return true;
    }
    bool Eof() const override { return !file || file.eof(); }
    std::string LastError() const override { return err; }
private:
    std::ifstream file;
    std::string err;
};

class GzInput : public FastaInput {
public:
    explicit GzInput(const std::filesystem::path& p) {
        gz = gzopen(p.c_str(), "rb");
        if (!gz) { err = "gzopen failed: " + p.string(); return; }
        gzbuffer(gz, 1u << 20);
    }
    ~GzInput() override { if (gz) gzclose(gz); }
    bool GetLine(std::string& out) override {
        if (!gz) return false;
        out.clear();
        char buf[8192];
        for (;;) {
            if (!gzgets(gz, buf, sizeof(buf))) {
                eof = gzeof(gz);
                if (!eof) {
                    int errnum = 0;
                    const char* msg = gzerror(gz, &errnum);
                    err = std::string("gzgets error: ") + (msg ? msg : "(nul)");
                }
                return !out.empty();
            }
            out.append(buf);
            // Stop once we've consumed a full line (newline at end) or EOF.
            if (!out.empty() && out.back() == '\n') return true;
            // If the buffer didn't contain a newline, loop to append more.
        }
    }
    bool Eof() const override { return eof || !gz; }
    std::string LastError() const override { return err; }
private:
    gzFile gz = nullptr;
    bool eof = false;
    std::string err;
};

}  // namespace

class FastaReaderImpl {
public:
    std::unique_ptr<FastaInput> input;
    bool at_eof = false;
    std::string last_error;
    size_t records_read = 0;
    std::string pending_header;
    std::string line_buffer;
};

FastaReader::FastaReader(const std::filesystem::path& path,
                         const FastaReaderConfig& config)
    : path_(path), config_(config), impl_(std::make_unique<FastaReaderImpl>()) {

    if (HasGzExtension(path_)) {
        impl_->input = std::make_unique<GzInput>(path_);
    } else {
        impl_->input = std::make_unique<IfstreamInput>(path_);
    }
    impl_->last_error = impl_->input->LastError();
    if (!impl_->last_error.empty()) {
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

    std::string header;
    if (!impl_->pending_header.empty()) {
        header = std::move(impl_->pending_header);
        impl_->pending_header.clear();
    } else {
        while (impl_->input->GetLine(header)) {
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

    size_t name_end = header.find_first_of(" \t", 1);
    if (name_end == std::string::npos) {
        record.name = header.substr(1);
    } else {
        record.name = header.substr(1, name_end - 1);
    }

    while (impl_->input->GetLine(impl_->line_buffer)) {
        TrimRight(impl_->line_buffer);

        if (impl_->line_buffer.empty()) {
            continue;
        }

        if (impl_->line_buffer[0] == '>') {
            impl_->pending_header = impl_->line_buffer;
            break;
        }

        record.sequence += impl_->line_buffer;
    }

    if (impl_->pending_header.empty() && impl_->input->Eof()) {
        impl_->at_eof = true;
    }

    if (config_.uppercase_sequence) {
        std::transform(record.sequence.begin(), record.sequence.end(),
                       record.sequence.begin(), ::toupper);
    }

    if (config_.skip_N_only && record.IsValid()) {
        bool all_n = std::all_of(record.sequence.begin(), record.sequence.end(),
                                 [](char c) { return c == 'N' || c == 'n'; });
        if (all_n) {
            return Next();
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
    return !impl_->input->Eof();
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
