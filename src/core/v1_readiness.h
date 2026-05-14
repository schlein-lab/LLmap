#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace llmap::core {

// Individual capability check result
struct CapabilityCheck {
    std::string name;
    std::string description;
    bool passed;
    std::optional<std::string> error_message;
    double duration_ms;
};

// Capability categories matching the 10-phase architecture
enum class CapabilityCategory {
    Core,          // Phase 0: AlignmentRecord, BucketPyramid, WaveState
    Foundation,    // Phase 1: Foundation model embeddings
    SelfInterference,  // Phase 2: FAISS, Leiden, Self-WaveCollapse
    ReferenceCollapse, // Phase 3: Reference index, EM iteration
    Classical,     // Phase 4: Minimizer index, chain, WFA2
    Validation,    // Phase 5: Kill-switch validation
    Output,        // Phase 6: BAM/SAM, Parquet
    Agent,         // Phase 7: Claude agent integration
    Performance,   // Phase 8: Profiling, arena, SIMD, mmap
    SingleCell,    // Phase 9: Cell barcode, PSV, QC
    Production     // Phase 10: Logging, error handling, config
};

std::string CategoryToString(CapabilityCategory category);
CapabilityCategory StringToCategory(std::string_view str);

// Results for a category of capabilities
struct CategoryResult {
    CapabilityCategory category;
    std::vector<CapabilityCheck> checks;
    bool all_passed;
    double total_duration_ms;
};

// Overall readiness report
struct ReadinessReport {
    std::string version;
    std::string build_date;
    std::string platform;
    std::vector<CategoryResult> categories;
    bool all_passed;
    uint32_t total_checks;
    uint32_t passed_checks;
    uint32_t failed_checks;
    double total_duration_ms;
};

// Configuration for readiness check
struct ReadinessConfig {
    bool check_core = true;
    bool check_foundation = true;
    bool check_self_interference = true;
    bool check_reference_collapse = true;
    bool check_classical = true;
    bool check_validation = true;
    bool check_output = true;
    bool check_agent = true;
    bool check_performance = true;
    bool check_single_cell = true;
    bool check_production = true;
    bool verbose = false;
    bool skip_gpu = true;  // Default to CPU-only checks
};

// Run the full V1.0 readiness check
ReadinessReport RunReadinessCheck(const ReadinessConfig& config = {});

// Run checks for a specific category
CategoryResult RunCategoryCheck(CapabilityCategory category,
                                const ReadinessConfig& config);

// Format report for display
std::string FormatReport(const ReadinessReport& report, bool verbose = false);

// Format report as JSON
std::string FormatReportJson(const ReadinessReport& report);

// Individual capability checks (used internally and for testing)
namespace checks {

// Phase 0: Core
CapabilityCheck CheckAlignmentRecord();
CapabilityCheck CheckBucketPyramid();
CapabilityCheck CheckWaveState();
CapabilityCheck CheckFastaReader();

// Phase 1: Foundation
CapabilityCheck CheckFoundationEmbedder();
CapabilityCheck CheckBucketEmbedder();

// Phase 2: Self-Interference
CapabilityCheck CheckFaissWrapper();
CapabilityCheck CheckSimilarityGraph();
CapabilityCheck CheckLeidenClustering();
CapabilityCheck CheckSelfWaveCollapse();
CapabilityCheck CheckClusterRep();

// Phase 3: Reference Collapse
CapabilityCheck CheckReferenceIndex();
CapabilityCheck CheckEmIterator();
CapabilityCheck CheckCollapseCheck();
CapabilityCheck CheckRefinement();
CapabilityCheck CheckMemberPropagation();

// Phase 4: Classical
CapabilityCheck CheckMinimizerIndex();
CapabilityCheck CheckChainExtraction();
CapabilityCheck CheckWfa2Aligner();
CapabilityCheck CheckClassicalPipeline();

// Phase 5: Validation
CapabilityCheck CheckKillSwitch();
CapabilityCheck CheckEndToEnd();
CapabilityCheck CheckRealReference();

// Phase 6: Output
CapabilityCheck CheckBamWriter();
CapabilityCheck CheckParquetWriter();
CapabilityCheck CheckParquetReader();

// Phase 7: Agent
CapabilityCheck CheckAgentTypes();
CapabilityCheck CheckAnthropicClient();
CapabilityCheck CheckAgentSession();
CapabilityCheck CheckCudaSandbox();
CapabilityCheck CheckPipelineAgent();

// Phase 8: Performance
CapabilityCheck CheckProfiler();
CapabilityCheck CheckArena();
CapabilityCheck CheckSimd();
CapabilityCheck CheckMmapFasta();
CapabilityCheck CheckThreadPool();
CapabilityCheck CheckCacheLayout();

// Phase 9: Single-Cell
CapabilityCheck CheckCbPreservation();
CapabilityCheck CheckPerCellParalog();
CapabilityCheck CheckPsvCatalog();
CapabilityCheck CheckPsvAssigner();
CapabilityCheck CheckScParalogQc();

// Phase 10: Production
CapabilityCheck CheckLogging();
CapabilityCheck CheckErrorHandling();
CapabilityCheck CheckConfig();
CapabilityCheck CheckVersion();

}  // namespace checks

}  // namespace llmap::core
