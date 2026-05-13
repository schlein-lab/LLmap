// LLmap — Parsing functions for real reference validation.
//
// Handles BED ground truth parsing and BAM baseline extraction.

#include "validation/real_reference.h"
#include "validation/real_reference_internal.h"

#include <array>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>

namespace llmap::validation {

namespace internal {

std::vector<std::string_view> SplitLine(std::string_view line, char delim) {
    std::vector<std::string_view> fields;
    size_t start = 0;
    while (start < line.size()) {
        auto pos = line.find(delim, start);
        if (pos == std::string_view::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

uint64_t ParseUint64(std::string_view s) {
    uint64_t val = 0;
    for (char c : s) {
        if (c >= '0' && c <= '9') {
            val = val * 10 + (c - '0');
        }
    }
    return val;
}

std::string TrimWhitespace(std::string_view s) {
    size_t start = 0;
    while (start < s.size() && (s[start] == ' ' || s[start] == '\t' ||
                                 s[start] == '\n' || s[start] == '\r')) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && (s[end-1] == ' ' || s[end-1] == '\t' ||
                           s[end-1] == '\n' || s[end-1] == '\r')) {
        --end;
    }
    return std::string(s.substr(start, end - start));
}

}  // namespace internal

std::vector<RealGroundTruth> ParseGroundTruthBed(
    const std::filesystem::path& bed_path) {

    std::vector<RealGroundTruth> truths;
    std::ifstream in(bed_path);
    if (!in.is_open()) {
        return truths;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        auto fields = internal::SplitLine(line, '\t');
        if (fields.size() < 4) continue;

        RealGroundTruth truth;
        truth.chrom = std::string(fields[0]);
        truth.start = internal::ParseUint64(fields[1]);
        truth.end = internal::ParseUint64(fields[2]);
        truth.read_id = std::string(fields[3]);

        if (fields.size() > 4) {
            truth.strand = std::string(fields[5]);
        }
        if (fields.size() > 6) {
            truth.gene = std::string(fields[6]);
        }

        truths.push_back(std::move(truth));
    }

    return truths;
}

std::unordered_map<std::string, bool> ParseMinimap2Bam(
    const std::filesystem::path& bam_path) {

    std::unordered_map<std::string, bool> mapped;

    std::string cmd = "samtools view " + bam_path.string() +
                      " 2>/dev/null | cut -f1,2";

    std::array<char, 4096> buffer;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(
        popen(cmd.c_str(), "r"), pclose);

    if (!pipe) {
        return mapped;
    }

    std::string output;
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output += buffer.data();
    }

    std::istringstream iss(output);
    std::string line;
    while (std::getline(iss, line)) {
        auto fields = internal::SplitLine(line, '\t');
        if (fields.size() < 2) continue;

        std::string read_name(fields[0]);
        uint64_t flag = internal::ParseUint64(fields[1]);

        bool is_unmapped = (flag & 0x4) != 0;
        mapped[read_name] = !is_unmapped;
    }

    return mapped;
}

}  // namespace llmap::validation
