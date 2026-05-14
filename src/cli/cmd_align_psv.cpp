// LLmap — cmd_align PSV integration.
// PSV-based paralog disambiguation for the align pipeline.

#include "cli/cmd_align_internal.h"

#include <cstdio>
#include <filesystem>
#include <optional>

#include "psv/psv_assigner.h"
#include "psv/psv_catalog.h"

namespace llmap::cli::align_internal {

std::optional<psv::PsvCatalog> LoadPsvCatalog(const AlignArgs& args) {
    if (args.psv_catalog.empty()) {
        return std::nullopt;
    }

    if (!std::filesystem::exists(args.psv_catalog)) {
        std::fprintf(stderr, "Error: PSV catalog not found: %s\n",
                     args.psv_catalog.c_str());
        return std::nullopt;
    }

    std::filesystem::path path(args.psv_catalog);
    std::string ext = path.extension().string();

    std::optional<psv::PsvCatalog> catalog;

    if (ext == ".vcf" || ext == ".vcf.gz") {
        catalog = psv::LoadPsvCatalogFromVcf(path);
    } else {
        catalog = psv::LoadPsvCatalogFromBed(path);
    }

    if (!catalog) {
        std::fprintf(stderr, "Error: failed to load PSV catalog: %s\n",
                     args.psv_catalog.c_str());
        return std::nullopt;
    }

    catalog->BuildIndex();
    return catalog;
}

psv::PsvAssignmentConfig CreatePsvConfig(const AlignArgs& args) {
    psv::PsvAssignmentConfig config;
    config.confidence_threshold = args.psv_min_posterior;
    config.error_rate = 0.01f;
    return config;
}

void ApplyPsvAssignments(
    const psv::PsvCatalog& catalog,
    const AlignArgs& args,
    std::vector<AlignmentRecord>& records,
    const std::vector<std::string>& read_sequences,
    bool verbose)
{
    psv::PsvAssigner assigner(catalog, CreatePsvConfig(args));

    std::size_t n_assigned = 0;
    std::size_t n_confident = 0;

    for (std::size_t i = 0; i < records.size(); ++i) {
        auto& record = records[i];

        if (!record.primary) {
            continue;
        }

        if (i >= read_sequences.size()) {
            continue;
        }

        auto observations = assigner.ExtractObservations(*record.primary, read_sequences[i]);

        if (observations.empty()) {
            continue;
        }

        auto result = assigner.Assign(record.read_id, observations);

        if (result.likelihoods.empty()) {
            continue;
        }

        n_assigned++;

        auto psv_call = psv::ResultToParalogCall(result);

        if (args.psv_only || !record.paralog_assignment) {
            record.paralog_assignment = psv_call;
        } else {
            record.paralog_assignment = psv::MergeAssignments(*record.paralog_assignment, result, args.psv_weight);
        }

        if (result.best_posterior >= args.psv_min_posterior) {
            n_confident++;
        }
    }

    if (verbose) {
        const auto& stats = assigner.GetStats();
        std::fprintf(stderr, "PSV assignment:\n");
        std::fprintf(stderr, "  Reads with PSV observations: %zu\n", n_assigned);
        std::fprintf(stderr, "  Confident paralog calls:     %zu\n", n_confident);
        std::fprintf(stderr, "  Total PSV sites queried:     %zu\n", stats.total_psv_observations);
        std::fprintf(stderr, "  Avg observations per read:   %.1f\n",
                     n_assigned > 0 ? static_cast<float>(stats.total_psv_observations) / n_assigned : 0.0f);
    }
}

}  // namespace llmap::cli::align_internal
