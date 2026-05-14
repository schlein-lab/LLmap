// LLmap — Configuration file support — TOML parsing.

#include "core/config.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace llmap {

namespace {

std::string_view Trim(std::string_view s) noexcept {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
        s.remove_prefix(1);
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.remove_suffix(1);
    }
    return s;
}

std::string Unquote(std::string_view s) {
    s = Trim(s);
    if (s.size() >= 2) {
        if ((s.front() == '"' && s.back() == '"') ||
            (s.front() == '\'' && s.back() == '\'')) {
            s = s.substr(1, s.size() - 2);
        }
    }
    return std::string(s);
}

bool IsBoolTrue(std::string_view s) noexcept {
    s = Trim(s);
    return s == "true" || s == "yes" || s == "on" || s == "1";
}

bool IsBoolFalse(std::string_view s) noexcept {
    s = Trim(s);
    return s == "false" || s == "no" || s == "off" || s == "0";
}

}  // namespace

bool ConfigValue::AsBool(bool default_val) const noexcept {
    if (IsBoolTrue(value)) return true;
    if (IsBoolFalse(value)) return false;
    return default_val;
}

int64_t ConfigValue::AsInt(int64_t default_val) const noexcept {
    auto s = Trim(value);
    int64_t result = default_val;
    std::from_chars(s.data(), s.data() + s.size(), result);
    return result;
}

double ConfigValue::AsDouble(double default_val) const noexcept {
    auto s = Trim(value);
    double result = default_val;
    std::from_chars(s.data(), s.data() + s.size(), result);
    return result;
}

LLmapConfig LLmapConfig::Defaults() {
    return LLmapConfig{};
}

LLmapError ConfigParser::MakeParseError(std::string_view what, size_t line) const {
    return ParseError(ErrorCode::kParseInvalidFormat, what, line);
}

Result<void, LLmapError> ConfigParser::ParseLine(std::string_view line, size_t line_num) {
    line = Trim(line);

    if (line.empty() || line.front() == '#') {
        return MakeOk();
    }

    if (line.front() == '[') {
        auto end = line.find(']');
        if (end == std::string_view::npos) {
            return MakeErr<void>(MakeParseError("unclosed section bracket", line_num));
        }
        current_section_ = std::string(Trim(line.substr(1, end - 1)));
        return MakeOk();
    }

    auto eq_pos = line.find('=');
    if (eq_pos == std::string_view::npos) {
        return MakeErr<void>(MakeParseError("missing '=' in key-value pair", line_num));
    }

    auto key = Trim(line.substr(0, eq_pos));
    auto val = Trim(line.substr(eq_pos + 1));

    if (key.empty()) {
        return MakeErr<void>(MakeParseError("empty key", line_num));
    }

    ConfigValue cv;
    cv.key = std::string(key);
    cv.value = Unquote(val);
    cv.section = current_section_;
    cv.line = line_num;

    values_.push_back(std::move(cv));
    return MakeOk();
}

Result<LLmapConfig, LLmapError> ConfigParser::Parse(std::string_view content,
                                                     std::string_view filename) {
    filename_ = std::string(filename);
    values_.clear();
    current_section_.clear();

    size_t line_num = 0;
    size_t pos = 0;

    while (pos < content.size()) {
        ++line_num;
        auto end = content.find('\n', pos);
        if (end == std::string_view::npos) {
            end = content.size();
        }

        auto line = content.substr(pos, end - pos);
        auto result = ParseLine(line, line_num);
        if (!result.ok()) {
            return MakeErr<LLmapConfig>(std::move(result.error()));
        }

        pos = end + 1;
    }

    LLmapConfig config = LLmapConfig::Defaults();
    config.source_file = filename_;

    for (const auto& val : values_) {
        if (val.section == "align") {
            if (val.key == "threads") config.align.threads = static_cast<uint32_t>(val.AsInt(0));
            else if (val.key == "kmer_size") config.align.kmer_size = static_cast<uint32_t>(val.AsInt(15));
            else if (val.key == "window_size") config.align.window_size = static_cast<uint32_t>(val.AsInt(10));
            else if (val.key == "min_chain_score") config.align.min_chain_score = static_cast<uint32_t>(val.AsInt(50));
            else if (val.key == "max_gaps_in_chain") config.align.max_gaps_in_chain = static_cast<uint32_t>(val.AsInt(5));
            else if (val.key == "min_identity") config.align.min_identity = val.AsDouble(0.9);
            else if (val.key == "min_mapping_quality") config.align.min_mapping_quality = val.AsDouble(10.0);
            else if (val.key == "output_format") config.align.output_format = std::string(val.AsString());
            else if (val.key == "include_unmapped") config.align.include_unmapped = val.AsBool(false);
            else if (val.key == "include_secondary") config.align.include_secondary = val.AsBool(true);
        } else if (val.section == "llm") {
            if (val.key == "api_key") config.llm.api_key = std::string(val.AsString());
            else if (val.key == "model") config.llm.model = std::string(val.AsString());
            else if (val.key == "temperature") config.llm.temperature = val.AsDouble(0.7);
            else if (val.key == "max_tokens") config.llm.max_tokens = static_cast<uint32_t>(val.AsInt(4096));
            else if (val.key == "stall_threshold") config.llm.stall_threshold = val.AsDouble(0.1);
            else if (val.key == "work_dir") config.llm.work_dir = std::string(val.AsString());
            else if (val.key == "enabled") config.llm.enabled = val.AsBool(false);
        } else if (val.section == "singlecell") {
            if (val.key == "cb_tag") config.singlecell.cb_tag = std::string(val.AsString());
            else if (val.key == "umi_tag") config.singlecell.umi_tag = std::string(val.AsString());
            else if (val.key == "cb_pattern") config.singlecell.cb_pattern = std::string(val.AsString());
            else if (val.key == "cb_whitelist_file") config.singlecell.cb_whitelist_file = std::string(val.AsString());
            else if (val.key == "platform") config.singlecell.platform = std::string(val.AsString());
            else if (val.key == "min_reads_per_cell") config.singlecell.min_reads_per_cell = static_cast<uint32_t>(val.AsInt(100));
            else if (val.key == "min_confidence") config.singlecell.min_confidence = val.AsDouble(0.5);
        } else if (val.section == "psv") {
            if (val.key == "catalog_file") config.psv.catalog_file = std::string(val.AsString());
            else if (val.key == "weight") config.psv.weight = val.AsDouble(0.5);
            else if (val.key == "min_posterior") config.psv.min_posterior = val.AsDouble(0.9);
            else if (val.key == "enabled") config.psv.enabled = val.AsBool(false);
        } else if (val.section == "logging") {
            if (val.key == "level") config.logging.level = ParseLogLevel(val.AsString());
            else if (val.key == "format") {
                if (val.AsString() == "json") config.logging.format = LogFormat::kJson;
                else config.logging.format = LogFormat::kText;
            }
            else if (val.key == "file") config.logging.file = std::string(val.AsString());
        } else if (val.section.empty() || val.section == "general") {
            if (val.key == "reference") config.reference_path = std::string(val.AsString());
            else if (val.key == "index") config.index_path = std::string(val.AsString());
            else if (val.key == "temp_dir") config.temp_dir = std::string(val.AsString());
        }
    }

    return MakeOk(std::move(config));
}

Result<LLmapConfig, LLmapError> ConfigParser::ParseFile(const std::filesystem::path& path) {
    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) {
        return MakeErr<LLmapConfig>(IoError(ErrorCode::kIoFileNotFound, path.string()));
    }

    std::ifstream file(path);
    if (!file) {
        return MakeErr<LLmapConfig>(IoError(ErrorCode::kIoReadError, path.string()));
    }

    std::ostringstream oss;
    oss << file.rdbuf();
    return Parse(oss.str(), path.string());
}

std::optional<ConfigValue> ConfigParser::Get(std::string_view section,
                                              std::string_view key) const noexcept {
    for (const auto& v : values_) {
        if (v.section == section && v.key == key) {
            return v;
        }
    }
    return std::nullopt;
}

}  // namespace llmap
