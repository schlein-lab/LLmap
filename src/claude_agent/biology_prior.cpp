#include "biology_prior.h"

#include <fstream>
#include <sstream>

namespace llmap::claude_agent {

namespace {

std::string EscapeJsonString(std::string_view s) {
    std::ostringstream oss;
    for (char c : s) {
        switch (c) {
            case '"':  oss << "\\\""; break;
            case '\\': oss << "\\\\"; break;
            case '\n': oss << "\\n"; break;
            case '\r': oss << "\\r"; break;
            case '\t': oss << "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    oss << "\\u" << std::hex << std::setfill('0')
                        << std::setw(4) << static_cast<int>(c);
                } else {
                    oss << c;
                }
        }
    }
    return oss.str();
}

}  // namespace

std::string SerializeBiologyPrior(const BiologyPrior& prior) {
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"version\": \"" << EscapeJsonString(prior.version) << "\",\n";
    oss << "  \"reference_sha256\": \""
        << EscapeJsonString(prior.reference_sha256) << "\",\n";

    oss << "  \"buckets\": {\n";
    bool first_bucket = true;
    for (const auto& [id, ann] : prior.buckets) {
        if (!first_bucket) {
            oss << ",\n";
        }
        first_bucket = false;
        oss << "    \"" << id << "\": {\n";
        oss << "      \"level\": \"" << EscapeJsonString(ann.level) << "\",\n";
        oss << "      \"annotation\": \""
            << EscapeJsonString(ann.annotation) << "\",\n";
        oss << "      \"prior_weight\": " << ann.prior_weight << ",\n";
        if (ann.paralog_partner_bucket) {
            oss << "      \"paralog_partner_bucket\": "
                << *ann.paralog_partner_bucket << ",\n";
        }
        oss << "      \"expected_coverage_multiplier\": "
            << ann.expected_coverage_multiplier << ",\n";
        oss << "      \"claude_rationale\": \""
            << EscapeJsonString(ann.claude_rationale) << "\"\n";
        oss << "    }";
    }
    oss << "\n  },\n";

    oss << "  \"regional_overrides\": {\n";
    bool first_region = true;
    for (const auto& [region, ovr] : prior.regional_overrides) {
        if (!first_region) {
            oss << ",\n";
        }
        first_region = false;
        oss << "    \"" << EscapeJsonString(region) << "\": {\n";
        oss << "      \"sub_bucket_granularity_kb\": "
            << ovr.sub_bucket_granularity_kb << ",\n";
        oss << "      \"max_iter\": " << ovr.max_iter << ",\n";
        oss << "      \"convergence_threshold\": "
            << ovr.convergence_threshold << "\n";
        oss << "    }";
    }
    oss << "\n  }\n";
    oss << "}\n";

    return oss.str();
}

BiologyPrior DeserializeBiologyPrior(std::string_view json) {
    BiologyPrior prior;
    // Simple JSON parsing - in production would use a proper JSON library
    // For now, return empty prior
    // TODO: implement proper JSON parsing or use nlohmann/json
    (void)json;
    return prior;
}

bool WriteBiologyPrior(
    const BiologyPrior& prior,
    const std::filesystem::path& output_path
) {
    std::ofstream ofs(output_path);
    if (!ofs) {
        return false;
    }
    ofs << SerializeBiologyPrior(prior);
    return ofs.good();
}

std::optional<BiologyPrior> ReadBiologyPrior(
    const std::filesystem::path& input_path
) {
    std::ifstream ifs(input_path);
    if (!ifs) {
        return std::nullopt;
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return DeserializeBiologyPrior(oss.str());
}

std::string SerializeSampleParams(const SampleParams& params) {
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"preset\": \"" << EscapeJsonString(params.preset) << "\",\n";
    oss << "  \"convergence_threshold\": "
        << params.convergence_threshold << ",\n";
    oss << "  \"max_iterations\": " << params.max_iterations << ",\n";
    oss << "  \"expected_coverage_profile\": \""
        << EscapeJsonString(params.expected_coverage_profile) << "\",\n";
    oss << "  \"foundation_model\": \""
        << EscapeJsonString(params.foundation_model) << "\",\n";

    oss << "  \"region_adjustments\": {\n";
    bool first = true;
    for (const auto& [region, adj] : params.region_adjustments) {
        if (!first) {
            oss << ",\n";
        }
        first = false;
        oss << "    \"" << EscapeJsonString(region) << "\": " << adj;
    }
    oss << "\n  }\n";
    oss << "}\n";

    return oss.str();
}

SampleParams DeserializeSampleParams(std::string_view json) {
    SampleParams params;
    // TODO: implement proper JSON parsing
    (void)json;
    return params;
}

bool WriteSampleParams(
    const SampleParams& params,
    const std::filesystem::path& output_path
) {
    std::ofstream ofs(output_path);
    if (!ofs) {
        return false;
    }
    ofs << SerializeSampleParams(params);
    return ofs.good();
}

std::optional<SampleParams> ReadSampleParams(
    const std::filesystem::path& input_path
) {
    std::ifstream ifs(input_path);
    if (!ifs) {
        return std::nullopt;
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return DeserializeSampleParams(oss.str());
}

}  // namespace llmap::claude_agent
