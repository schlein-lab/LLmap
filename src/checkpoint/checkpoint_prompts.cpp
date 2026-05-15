// LLmap -- prompt templates and response parsing for checkpoints.

#include "checkpoint/checkpoint_prompts.h"
#include "checkpoint/checkpoint_cache.h"  // reuse DeserializeDecision

#include <sstream>

namespace llmap::checkpoint {

namespace {

void AppendCommonContext(std::ostringstream& oss, const CheckpointContext& ctx) {
    oss << "read_id: " << ctx.read_id << "\n"
        << "ref_id: " << ctx.ref_id << "\n";
    if (!ctx.region_name.empty()) {
        oss << "region: " << ctx.region_name << "\n";
    }
    if (!ctx.user_note.empty()) {
        oss << "note: " << ctx.user_note << "\n";
    }
    if (!ctx.candidate_positions.empty()) {
        oss << "candidates (pos, chain_score):\n";
        for (const auto& [p, s] : ctx.candidate_positions) {
            oss << "  - " << p << " " << s << "\n";
        }
    }
    if (!ctx.read_seq.empty()) {
        // Truncate at 1 kb for prompt economy.
        size_t n = std::min<size_t>(ctx.read_seq.size(), 1024);
        oss << "read_seq (first " << n << " bp): "
            << ctx.read_seq.substr(0, n) << "\n";
    }
}

std::string AmbiguousChainPrompt(const CheckpointContext& ctx) {
    std::ostringstream oss;
    oss << "CHECKPOINT: ambiguous-chain\n"
        << "Multiple chains are within 5% of the best chain score. The\n"
        << "classical mapper cannot pick one with confidence.\n\n";
    AppendCommonContext(oss, ctx);
    oss << "\nTask: identify which candidate is most plausible. If the read\n"
        << "spans a paralog family or tandem repeat, emit a wave of multiple\n"
        << "positions with amplitudes summing to 1.0 instead of forcing a\n"
        << "single answer. Explain your reasoning.\n";
    return oss.str();
}

std::string UnknownRegionPrompt(const CheckpointContext& ctx) {
    std::ostringstream oss;
    oss << "CHECKPOINT: unknown-region\n"
        << "The reference window where this read lands has no taxonomy or\n"
        << "specific-locus annotation. The mapper is using default params.\n\n";
    AppendCommonContext(oss, ctx);
    oss << "\nTask: hypothesise what this region is (centromere? novel SD?\n"
        << "contamination?) and propose a ParamOverride if defaults are\n"
        << "likely wrong. Use region_lookup or local_grep if useful.\n";
    return oss.str();
}

std::string ParalogDisambigPrompt(const CheckpointContext& ctx) {
    std::ostringstream oss;
    oss << "CHECKPOINT: paralog-disambiguation\n"
        << "This region is flagged require_psv_disambig. The classical\n"
        << "scorer cannot tell paralogs apart -- they are ~99% identical.\n\n";
    AppendCommonContext(oss, ctx);
    oss << "\nTask: consult psv_check for known paralog-specific variants in\n"
        << "this region, decide which paralog the read belongs to, or emit\n"
        << "a multi-position wave if the read is genuinely ambiguous.\n";
    return oss.str();
}

std::string SDExpansionPrompt(const CheckpointContext& ctx) {
    std::ostringstream oss;
    oss << "CHECKPOINT: sd-expansion\n"
        << "The read is more than 10% longer than the reference window.\n"
        << "This suggests a segmental duplication expansion not present in\n"
        << "the reference.\n\n";
    AppendCommonContext(oss, ctx);
    oss << "\nTask: is this a known SD-expansion locus (region_lookup)?\n"
        << "If so, what is the expected copy number? Emit a 'special_finding'\n"
        << "describing the SD class if you can identify it.\n";
    return oss.str();
}

std::string NovelInsertionPrompt(const CheckpointContext& ctx) {
    std::ostringstream oss;
    oss << "CHECKPOINT: novel-insertion\n"
        << "A large soft-clip remains after extension. The read may carry a\n"
        << "novel insertion (transposable element, SV, contamination).\n\n";
    AppendCommonContext(oss, ctx);
    oss << "\nTask: classify the soft-clipped tail. Use local_grep or\n"
        << "web_fetch_ucsc to compare against known repeat elements. Emit a\n"
        << "'special_finding' describing the insertion class.\n";
    return oss.str();
}

}  // namespace

std::string SystemPrompt() {
    return
        "You are LLmap's Layer 3 active-mapping consultant. The chain DP /\n"
        "extension stage will hand you cases its deterministic scorer cannot\n"
        "resolve. Available tools: bash, region_lookup, psv_check,\n"
        "local_grep, web_fetch_ucsc.\n"
        "\n"
        "Reply with a single JSON object of the shape:\n"
        "{\n"
        "  \"schema_version\": 1,\n"
        "  \"consulted\": true,\n"
        "  \"fallback_used\": false,\n"
        "  \"reasoning\": \"<one sentence>\",\n"
        "  \"special_finding\": \"<short label or empty>\",\n"
        "  \"wave\": [[pos, amplitude], ...],\n"
        "  \"override\": {<optional ParamOverride fields>}\n"
        "}\n"
        "\n"
        "If you cannot decide, set fallback_used=true and leave wave empty.";
}

std::string BuildPrompt(CheckpointType type, const CheckpointContext& ctx) {
    switch (type) {
        case CheckpointType::AmbiguousChain:        return AmbiguousChainPrompt(ctx);
        case CheckpointType::UnknownRegion:         return UnknownRegionPrompt(ctx);
        case CheckpointType::ParalogDisambiguation: return ParalogDisambigPrompt(ctx);
        case CheckpointType::SDExpansion:           return SDExpansionPrompt(ctx);
        case CheckpointType::NovelInsertion:        return NovelInsertionPrompt(ctx);
    }
    return AmbiguousChainPrompt(ctx);
}

AgentDecision ParseAgentResponse(const std::string& response) {
    // Extract the first { ... } block; the agent may surround it with prose.
    size_t a = response.find('{');
    size_t b = response.rfind('}');
    if (a == std::string::npos || b == std::string::npos || b <= a) {
        AgentDecision d;
        d.fallback_used = true;
        d.reasoning = "no JSON in response";
        return d;
    }
    auto parsed = DeserializeDecision(response.substr(a, b - a + 1));
    if (parsed) return *parsed;

    AgentDecision d;
    d.fallback_used = true;
    d.reasoning = "JSON parse failed";
    return d;
}

}  // namespace llmap::checkpoint
