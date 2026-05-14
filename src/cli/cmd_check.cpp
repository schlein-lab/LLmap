#include "commands.h"

#include <cstring>
#include <iostream>
#include <string>

#include "../core/v1_readiness.h"

namespace llmap::cli {

namespace {

void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " check [options]\n\n"
              << "Run V1.0 readiness check to verify all capabilities.\n\n"
              << "Options:\n"
              << "  -v, --verbose     Show all checks, not just failures\n"
              << "  -j, --json        Output as JSON\n"
              << "  -c, --category    Check only specific category\n"
              << "                    (Core, Foundation, SelfInterference,\n"
              << "                     ReferenceCollapse, Classical, Validation,\n"
              << "                     Output, Agent, Performance, SingleCell,\n"
              << "                     Production)\n"
              << "  -h, --help        Show this help message\n\n"
              << "Examples:\n"
              << "  " << prog << " check              # Run all checks\n"
              << "  " << prog << " check -v           # Verbose output\n"
              << "  " << prog << " check -j           # JSON output\n"
              << "  " << prog << " check -c Core      # Check only Core category\n";
}

}  // namespace

int run_check(int argc, char** argv) {
    core::ReadinessConfig config;
    bool verbose = false;
    bool json_output = false;
    std::string category_filter;

    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "-h") == 0 ||
            std::strcmp(argv[i], "--help") == 0) {
            PrintUsage(argv[0]);
            return 0;
        }
        if (std::strcmp(argv[i], "-v") == 0 ||
            std::strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
            config.verbose = true;
            continue;
        }
        if (std::strcmp(argv[i], "-j") == 0 ||
            std::strcmp(argv[i], "--json") == 0) {
            json_output = true;
            continue;
        }
        if ((std::strcmp(argv[i], "-c") == 0 ||
             std::strcmp(argv[i], "--category") == 0) && i + 1 < argc) {
            category_filter = argv[++i];
            // Disable all, then enable only the requested category
            config.check_core = false;
            config.check_foundation = false;
            config.check_self_interference = false;
            config.check_reference_collapse = false;
            config.check_classical = false;
            config.check_validation = false;
            config.check_output = false;
            config.check_agent = false;
            config.check_performance = false;
            config.check_single_cell = false;
            config.check_production = false;

            if (category_filter == "Core") config.check_core = true;
            else if (category_filter == "Foundation") config.check_foundation = true;
            else if (category_filter == "SelfInterference") config.check_self_interference = true;
            else if (category_filter == "ReferenceCollapse") config.check_reference_collapse = true;
            else if (category_filter == "Classical") config.check_classical = true;
            else if (category_filter == "Validation") config.check_validation = true;
            else if (category_filter == "Output") config.check_output = true;
            else if (category_filter == "Agent") config.check_agent = true;
            else if (category_filter == "Performance") config.check_performance = true;
            else if (category_filter == "SingleCell") config.check_single_cell = true;
            else if (category_filter == "Production") config.check_production = true;
            else {
                std::cerr << "Error: Unknown category '" << category_filter << "'\n";
                return 1;
            }
            continue;
        }
        std::cerr << "Error: Unknown option '" << argv[i] << "'\n";
        PrintUsage(argv[0]);
        return 1;
    }

    auto report = core::RunReadinessCheck(config);

    if (json_output) {
        std::cout << core::FormatReportJson(report);
    } else {
        std::cout << core::FormatReport(report, verbose);
    }

    return report.all_passed ? 0 : 1;
}

}  // namespace llmap::cli
