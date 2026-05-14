// LLmap — Configuration file support — utilities and validation.

#include "core/config.h"

#include <cstdlib>
#include <sstream>

namespace llmap {

namespace {

const char* GetEnvOrNull(const char* name) {
    return std::getenv(name);
}

std::filesystem::path GetHomeDir() {
    if (auto* home = GetEnvOrNull("HOME")) {
        return std::filesystem::path(home);
    }
    if (auto* userprofile = GetEnvOrNull("USERPROFILE")) {
        return std::filesystem::path(userprofile);
    }
    return {};
}

bool IsBoolTrue(std::string_view s) noexcept {
    return s == "true" || s == "yes" || s == "on" || s == "1";
}

}  // namespace

std::vector<std::filesystem::path> GetConfigSearchPaths() {
    std::vector<std::filesystem::path> paths;

    paths.emplace_back("./llmap.toml");
    paths.emplace_back("./.llmap.toml");

    auto home = GetHomeDir();
    if (!home.empty()) {
        paths.push_back(home / ".config" / "llmap" / "config.toml");
        paths.push_back(home / ".llmap.toml");
    }

    paths.emplace_back("/etc/llmap/config.toml");

    return paths;
}

Result<std::filesystem::path, LLmapError> FindConfigFile() {
    std::error_code ec;
    for (const auto& path : GetConfigSearchPaths()) {
        if (std::filesystem::exists(path, ec) && std::filesystem::is_regular_file(path, ec)) {
            return MakeOk(path);
        }
    }
    return MakeErr<std::filesystem::path>(
        IoError(ErrorCode::kIoFileNotFound, "no config file found"));
}

Result<LLmapConfig, LLmapError> LoadConfigFromFile(const std::filesystem::path& path) {
    ConfigParser parser;
    return parser.ParseFile(path);
}

Result<LLmapConfig, LLmapError> LoadConfig() {
    auto path_result = FindConfigFile();
    if (!path_result.ok()) {
        return MakeOk(LLmapConfig::Defaults());
    }
    return LoadConfigFromFile(path_result.value());
}

void ApplyEnvironmentOverrides(LLmapConfig& config) {
    if (auto* val = GetEnvOrNull("LLMAP_REFERENCE")) {
        config.reference_path = val;
    }
    if (auto* val = GetEnvOrNull("LLMAP_INDEX")) {
        config.index_path = val;
    }
    if (auto* val = GetEnvOrNull("LLMAP_TEMP_DIR")) {
        config.temp_dir = val;
    }
    if (auto* val = GetEnvOrNull("LLMAP_THREADS")) {
        config.align.threads = static_cast<uint32_t>(std::atoi(val));
    }
    if (auto* val = GetEnvOrNull("ANTHROPIC_API_KEY")) {
        config.llm.api_key = val;
    }
    if (auto* val = GetEnvOrNull("LLMAP_LLM_MODEL")) {
        config.llm.model = val;
    }
    if (auto* val = GetEnvOrNull("LLMAP_LOG_LEVEL")) {
        config.logging.level = ParseLogLevel(val);
    }
    if (auto* val = GetEnvOrNull("LLMAP_LOG_FORMAT")) {
        if (std::string_view(val) == "json") {
            config.logging.format = LogFormat::kJson;
        }
    }
}

void ApplyOverrides(LLmapConfig& config, const std::vector<ConfigOverride>& overrides) {
    for (const auto& ov : overrides) {
        if (ov.key == "threads") config.align.threads = static_cast<uint32_t>(std::atoi(ov.value.c_str()));
        else if (ov.key == "kmer_size") config.align.kmer_size = static_cast<uint32_t>(std::atoi(ov.value.c_str()));
        else if (ov.key == "window_size") config.align.window_size = static_cast<uint32_t>(std::atoi(ov.value.c_str()));
        else if (ov.key == "min_identity") config.align.min_identity = std::atof(ov.value.c_str());
        else if (ov.key == "output_format") config.align.output_format = ov.value;
        else if (ov.key == "reference") config.reference_path = ov.value;
        else if (ov.key == "index") config.index_path = ov.value;
        else if (ov.key == "temp_dir") config.temp_dir = ov.value;
        else if (ov.key == "llm.api_key") config.llm.api_key = ov.value;
        else if (ov.key == "llm.model") config.llm.model = ov.value;
        else if (ov.key == "llm.enabled") config.llm.enabled = IsBoolTrue(ov.value);
        else if (ov.key == "log_level") config.logging.level = ParseLogLevel(ov.value);
    }
}

Result<void, LLmapError> ValidateConfig(const LLmapConfig& config) {
    ErrorList errors;

    if (config.align.kmer_size < 5 || config.align.kmer_size > 31) {
        errors.Add(ConfigError(ErrorCode::kConfigInvalidValue, "align.kmer_size must be 5-31"));
    }

    if (config.align.min_identity < 0.0 || config.align.min_identity > 1.0) {
        errors.Add(ConfigError(ErrorCode::kConfigInvalidValue, "align.min_identity must be 0.0-1.0"));
    }

    if (config.llm.temperature < 0.0 || config.llm.temperature > 2.0) {
        errors.Add(ConfigError(ErrorCode::kConfigInvalidValue, "llm.temperature must be 0.0-2.0"));
    }

    if (config.llm.enabled && config.llm.api_key.empty()) {
        errors.Add(ConfigError(ErrorCode::kConfigMissingRequired, "llm.api_key required when llm.enabled=true"));
    }

    if (config.psv.weight < 0.0 || config.psv.weight > 1.0) {
        errors.Add(ConfigError(ErrorCode::kConfigInvalidValue, "psv.weight must be 0.0-1.0"));
    }

    if (config.singlecell.min_confidence < 0.0 || config.singlecell.min_confidence > 1.0) {
        errors.Add(ConfigError(ErrorCode::kConfigInvalidValue, "singlecell.min_confidence must be 0.0-1.0"));
    }

    return errors.ToResult();
}

std::string ConfigToToml(const LLmapConfig& config) {
    std::ostringstream out;

    out << "# LLmap configuration file\n\n";

    if (!config.reference_path.empty()) {
        out << "reference = \"" << config.reference_path << "\"\n";
    }
    if (!config.index_path.empty()) {
        out << "index = \"" << config.index_path << "\"\n";
    }
    out << "temp_dir = \"" << config.temp_dir << "\"\n";

    out << "\n[align]\n";
    out << "threads = " << config.align.threads << "\n";
    out << "kmer_size = " << config.align.kmer_size << "\n";
    out << "window_size = " << config.align.window_size << "\n";
    out << "min_chain_score = " << config.align.min_chain_score << "\n";
    out << "max_gaps_in_chain = " << config.align.max_gaps_in_chain << "\n";
    out << "min_identity = " << config.align.min_identity << "\n";
    out << "min_mapping_quality = " << config.align.min_mapping_quality << "\n";
    out << "output_format = \"" << config.align.output_format << "\"\n";
    out << "include_unmapped = " << (config.align.include_unmapped ? "true" : "false") << "\n";
    out << "include_secondary = " << (config.align.include_secondary ? "true" : "false") << "\n";

    out << "\n[llm]\n";
    out << "enabled = " << (config.llm.enabled ? "true" : "false") << "\n";
    out << "model = \"" << config.llm.model << "\"\n";
    out << "temperature = " << config.llm.temperature << "\n";
    out << "max_tokens = " << config.llm.max_tokens << "\n";
    out << "stall_threshold = " << config.llm.stall_threshold << "\n";
    out << "work_dir = \"" << config.llm.work_dir << "\"\n";

    out << "\n[singlecell]\n";
    out << "cb_tag = \"" << config.singlecell.cb_tag << "\"\n";
    out << "umi_tag = \"" << config.singlecell.umi_tag << "\"\n";
    out << "platform = \"" << config.singlecell.platform << "\"\n";
    out << "min_reads_per_cell = " << config.singlecell.min_reads_per_cell << "\n";
    out << "min_confidence = " << config.singlecell.min_confidence << "\n";

    out << "\n[psv]\n";
    out << "enabled = " << (config.psv.enabled ? "true" : "false") << "\n";
    out << "weight = " << config.psv.weight << "\n";
    out << "min_posterior = " << config.psv.min_posterior << "\n";

    out << "\n[logging]\n";
    out << "level = \"" << LogLevelNameLower(config.logging.level) << "\"\n";
    out << "format = \"" << (config.logging.format == LogFormat::kJson ? "json" : "text") << "\"\n";

    return out.str();
}

}  // namespace llmap
