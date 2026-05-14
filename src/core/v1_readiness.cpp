#include "v1_readiness.h"

#include <chrono>
#include <sstream>
#include <iomanip>

#include "core/version.h"

namespace llmap::core {

std::string CategoryToString(CapabilityCategory category) {
    switch (category) {
        case CapabilityCategory::Core: return "Core";
        case CapabilityCategory::Foundation: return "Foundation";
        case CapabilityCategory::SelfInterference: return "SelfInterference";
        case CapabilityCategory::ReferenceCollapse: return "ReferenceCollapse";
        case CapabilityCategory::Classical: return "Classical";
        case CapabilityCategory::Validation: return "Validation";
        case CapabilityCategory::Output: return "Output";
        case CapabilityCategory::Agent: return "Agent";
        case CapabilityCategory::Performance: return "Performance";
        case CapabilityCategory::SingleCell: return "SingleCell";
        case CapabilityCategory::Production: return "Production";
    }
    return "Unknown";
}

CapabilityCategory StringToCategory(std::string_view str) {
    if (str == "Core") return CapabilityCategory::Core;
    if (str == "Foundation") return CapabilityCategory::Foundation;
    if (str == "SelfInterference") return CapabilityCategory::SelfInterference;
    if (str == "ReferenceCollapse") return CapabilityCategory::ReferenceCollapse;
    if (str == "Classical") return CapabilityCategory::Classical;
    if (str == "Validation") return CapabilityCategory::Validation;
    if (str == "Output") return CapabilityCategory::Output;
    if (str == "Agent") return CapabilityCategory::Agent;
    if (str == "Performance") return CapabilityCategory::Performance;
    if (str == "SingleCell") return CapabilityCategory::SingleCell;
    if (str == "Production") return CapabilityCategory::Production;
    return CapabilityCategory::Core;
}

CategoryResult RunCategoryCheck(CapabilityCategory category,
                                const ReadinessConfig& config) {
    CategoryResult result;
    result.category = category;
    result.all_passed = true;
    result.total_duration_ms = 0.0;

    auto add_check = [&](CapabilityCheck check) {
        if (!check.passed) result.all_passed = false;
        result.total_duration_ms += check.duration_ms;
        result.checks.push_back(std::move(check));
    };

    switch (category) {
        case CapabilityCategory::Core:
            if (config.check_core) {
                add_check(checks::CheckAlignmentRecord());
                add_check(checks::CheckBucketPyramid());
                add_check(checks::CheckWaveState());
                add_check(checks::CheckFastaReader());
            }
            break;
        case CapabilityCategory::Foundation:
            if (config.check_foundation) {
                add_check(checks::CheckFoundationEmbedder());
                add_check(checks::CheckBucketEmbedder());
            }
            break;
        case CapabilityCategory::SelfInterference:
            if (config.check_self_interference) {
                add_check(checks::CheckFaissWrapper());
                add_check(checks::CheckSimilarityGraph());
                add_check(checks::CheckLeidenClustering());
                add_check(checks::CheckSelfWaveCollapse());
                add_check(checks::CheckClusterRep());
            }
            break;
        case CapabilityCategory::ReferenceCollapse:
            if (config.check_reference_collapse) {
                add_check(checks::CheckReferenceIndex());
                add_check(checks::CheckEmIterator());
                add_check(checks::CheckCollapseCheck());
                add_check(checks::CheckRefinement());
                add_check(checks::CheckMemberPropagation());
            }
            break;
        case CapabilityCategory::Classical:
            if (config.check_classical) {
                add_check(checks::CheckMinimizerIndex());
                add_check(checks::CheckChainExtraction());
                add_check(checks::CheckWfa2Aligner());
                add_check(checks::CheckClassicalPipeline());
            }
            break;
        case CapabilityCategory::Validation:
            if (config.check_validation) {
                add_check(checks::CheckKillSwitch());
                add_check(checks::CheckEndToEnd());
                add_check(checks::CheckRealReference());
            }
            break;
        case CapabilityCategory::Output:
            if (config.check_output) {
                add_check(checks::CheckBamWriter());
                add_check(checks::CheckParquetWriter());
                add_check(checks::CheckParquetReader());
            }
            break;
        case CapabilityCategory::Agent:
            if (config.check_agent) {
                add_check(checks::CheckAgentTypes());
                add_check(checks::CheckAnthropicClient());
                add_check(checks::CheckAgentSession());
                add_check(checks::CheckCudaSandbox());
                add_check(checks::CheckPipelineAgent());
            }
            break;
        case CapabilityCategory::Performance:
            if (config.check_performance) {
                add_check(checks::CheckProfiler());
                add_check(checks::CheckArena());
                add_check(checks::CheckSimd());
                add_check(checks::CheckMmapFasta());
                add_check(checks::CheckThreadPool());
                add_check(checks::CheckCacheLayout());
            }
            break;
        case CapabilityCategory::SingleCell:
            if (config.check_single_cell) {
                add_check(checks::CheckCbPreservation());
                add_check(checks::CheckPerCellParalog());
                add_check(checks::CheckPsvCatalog());
                add_check(checks::CheckPsvAssigner());
                add_check(checks::CheckScParalogQc());
            }
            break;
        case CapabilityCategory::Production:
            if (config.check_production) {
                add_check(checks::CheckLogging());
                add_check(checks::CheckErrorHandling());
                add_check(checks::CheckConfig());
                add_check(checks::CheckVersion());
            }
            break;
    }

    return result;
}

ReadinessReport RunReadinessCheck(const ReadinessConfig& config) {
    ReadinessReport report;
    report.version = std::string(kVersion);
    report.build_date = __DATE__;

#if defined(__linux__)
    report.platform = "Linux";
#elif defined(__APPLE__)
    report.platform = "macOS";
#elif defined(_WIN32)
    report.platform = "Windows";
#else
    report.platform = "Unknown";
#endif

    report.all_passed = true;
    report.total_checks = 0;
    report.passed_checks = 0;
    report.failed_checks = 0;
    report.total_duration_ms = 0.0;

    std::vector<CapabilityCategory> categories = {
        CapabilityCategory::Core,
        CapabilityCategory::Foundation,
        CapabilityCategory::SelfInterference,
        CapabilityCategory::ReferenceCollapse,
        CapabilityCategory::Classical,
        CapabilityCategory::Validation,
        CapabilityCategory::Output,
        CapabilityCategory::Agent,
        CapabilityCategory::Performance,
        CapabilityCategory::SingleCell,
        CapabilityCategory::Production
    };

    for (auto cat : categories) {
        auto result = RunCategoryCheck(cat, config);
        for (const auto& check : result.checks) {
            report.total_checks++;
            if (check.passed) {
                report.passed_checks++;
            } else {
                report.failed_checks++;
                report.all_passed = false;
            }
        }
        report.total_duration_ms += result.total_duration_ms;
        report.categories.push_back(std::move(result));
    }

    return report;
}

std::string FormatReport(const ReadinessReport& report, bool verbose) {
    std::ostringstream ss;

    ss << "LLmap V1.0 Readiness Check\n";
    ss << "==========================\n\n";
    ss << "Version:    " << report.version << "\n";
    ss << "Build Date: " << report.build_date << "\n";
    ss << "Platform:   " << report.platform << "\n\n";

    for (const auto& cat : report.categories) {
        if (cat.checks.empty()) continue;

        ss << "[" << CategoryToString(cat.category) << "] ";
        ss << (cat.all_passed ? "PASS" : "FAIL");
        ss << " (" << std::fixed << std::setprecision(1);
        ss << cat.total_duration_ms << "ms)\n";

        if (verbose || !cat.all_passed) {
            for (const auto& check : cat.checks) {
                ss << "  " << (check.passed ? "[+]" : "[-]");
                ss << " " << check.name;
                if (!check.passed && check.error_message) {
                    ss << ": " << *check.error_message;
                }
                ss << "\n";
            }
        }
    }

    ss << "\nSummary: " << report.passed_checks << "/" << report.total_checks;
    ss << " checks passed";
    if (report.failed_checks > 0) {
        ss << " (" << report.failed_checks << " failed)";
    }
    ss << "\n";
    ss << "Total time: " << std::fixed << std::setprecision(1);
    ss << report.total_duration_ms << "ms\n";
    ss << "\nResult: " << (report.all_passed ? "READY" : "NOT READY") << "\n";

    return ss.str();
}

std::string FormatReportJson(const ReadinessReport& report) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2);

    ss << "{\n";
    ss << "  \"version\": \"" << report.version << "\",\n";
    ss << "  \"build_date\": \"" << report.build_date << "\",\n";
    ss << "  \"platform\": \"" << report.platform << "\",\n";
    ss << "  \"all_passed\": " << (report.all_passed ? "true" : "false") << ",\n";
    ss << "  \"total_checks\": " << report.total_checks << ",\n";
    ss << "  \"passed_checks\": " << report.passed_checks << ",\n";
    ss << "  \"failed_checks\": " << report.failed_checks << ",\n";
    ss << "  \"total_duration_ms\": " << report.total_duration_ms << ",\n";
    ss << "  \"categories\": [\n";

    bool first_cat = true;
    for (const auto& cat : report.categories) {
        if (cat.checks.empty()) continue;
        if (!first_cat) ss << ",\n";
        first_cat = false;

        ss << "    {\n";
        ss << "      \"name\": \"" << CategoryToString(cat.category) << "\",\n";
        ss << "      \"all_passed\": " << (cat.all_passed ? "true" : "false") << ",\n";
        ss << "      \"duration_ms\": " << cat.total_duration_ms << ",\n";
        ss << "      \"checks\": [\n";

        bool first_check = true;
        for (const auto& check : cat.checks) {
            if (!first_check) ss << ",\n";
            first_check = false;

            ss << "        {";
            ss << "\"name\": \"" << check.name << "\", ";
            ss << "\"passed\": " << (check.passed ? "true" : "false") << ", ";
            ss << "\"duration_ms\": " << check.duration_ms;
            if (check.error_message) {
                ss << ", \"error\": \"" << *check.error_message << "\"";
            }
            ss << "}";
        }
        ss << "\n      ]\n    }";
    }
    ss << "\n  ]\n}\n";

    return ss.str();
}

}  // namespace llmap::core
