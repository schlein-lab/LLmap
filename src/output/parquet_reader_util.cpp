// LLmap — Parquet reader: utility and free functions.

#include "output/parquet_reader.h"
#include "output/parquet_reader_impl.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace llmap::output {

namespace detail {

std::string_view Trim(std::string_view s) {
    while (!s.empty() && std::isspace(s.front())) s.remove_prefix(1);
    while (!s.empty() && std::isspace(s.back())) s.remove_suffix(1);
    return s;
}

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

}  // namespace detail

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
    line = detail::Trim(line);
    if (line.empty()) return std::nullopt;

    auto fields = detail::SplitCSV(line);
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
                     (std::fabs(orig.probability - re.probability) <=
                      probability_tolerance) &&
                     (std::fabs(orig.confidence - re.confidence) <=
                      probability_tolerance) &&
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
