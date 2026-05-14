// LLmap — Configuration file support.
//
// TOML-based configuration with cascading search paths:
//   1. ./llmap.toml (project local)
//   2. ~/.config/llmap/config.toml (user)
//   3. /etc/llmap/config.toml (system)
// CLI flags override config file values.

#pragma once

#include "core/error.h"
#include "core/logging.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <optional>
#include <vector>
#include <filesystem>
#include <unordered_map>

namespace llmap {

struct AlignConfig {
    uint32_t threads{0};
    uint32_t kmer_size{15};
    uint32_t window_size{10};
    uint32_t min_chain_score{50};
    uint32_t max_gaps_in_chain{5};
    double min_identity{0.9};
    double min_mapping_quality{10.0};
    std::string output_format{"bam"};
    bool include_unmapped{false};
    bool include_secondary{true};
};

struct LlmConfig {
    std::string api_key;
    std::string model{"claude-sonnet-4-20250514"};
    double temperature{0.7};
    uint32_t max_tokens{4096};
    double stall_threshold{0.1};
    std::string work_dir{"/tmp/llmap_llm"};
    bool enabled{false};
};

struct SingleCellConfig {
    std::string cb_tag{"CB"};
    std::string umi_tag{"UB"};
    std::string cb_pattern;
    std::string cb_whitelist_file;
    std::string platform{"10x-v3"};
    uint32_t min_reads_per_cell{100};
    double min_confidence{0.5};
};

struct PsvConfig {
    std::string catalog_file;
    double weight{0.5};
    double min_posterior{0.9};
    bool enabled{false};
};

struct LoggingConfig {
    LogLevel level{LogLevel::kInfo};
    LogFormat format{LogFormat::kText};
    std::string file;
};

struct LLmapConfig {
    AlignConfig align;
    LlmConfig llm;
    SingleCellConfig singlecell;
    PsvConfig psv;
    LoggingConfig logging;

    std::string reference_path;
    std::string index_path;
    std::string temp_dir{"/tmp/llmap"};

    std::filesystem::path source_file;

    static LLmapConfig Defaults();
};

struct ConfigValue {
    std::string key;
    std::string value;
    std::string section;
    size_t line{0};

    bool IsEmpty() const noexcept { return value.empty(); }
    bool AsBool(bool default_val = false) const noexcept;
    int64_t AsInt(int64_t default_val = 0) const noexcept;
    double AsDouble(double default_val = 0.0) const noexcept;
    std::string_view AsString() const noexcept { return value; }
};

class ConfigParser {
public:
    Result<LLmapConfig, LLmapError> Parse(std::string_view content,
                                          std::string_view filename = "<string>");

    Result<LLmapConfig, LLmapError> ParseFile(const std::filesystem::path& path);

    const std::vector<ConfigValue>& Values() const noexcept { return values_; }
    std::optional<ConfigValue> Get(std::string_view section,
                                   std::string_view key) const noexcept;

private:
    Result<void, LLmapError> ParseLine(std::string_view line, size_t line_num);
    LLmapError MakeParseError(std::string_view what, size_t line) const;

    std::vector<ConfigValue> values_;
    std::string current_section_;
    std::string filename_;
};

Result<std::filesystem::path, LLmapError> FindConfigFile();
std::vector<std::filesystem::path> GetConfigSearchPaths();

Result<LLmapConfig, LLmapError> LoadConfig();
Result<LLmapConfig, LLmapError> LoadConfigFromFile(const std::filesystem::path& path);

void ApplyEnvironmentOverrides(LLmapConfig& config);

struct ConfigOverride {
    std::string key;
    std::string value;
};

void ApplyOverrides(LLmapConfig& config, const std::vector<ConfigOverride>& overrides);

std::string ConfigToToml(const LLmapConfig& config);

Result<void, LLmapError> ValidateConfig(const LLmapConfig& config);

}  // namespace llmap
