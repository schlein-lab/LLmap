// LLmap -- checkpoint tool registry, with stub implementations.

#include "checkpoint/checkpoint_tools.h"

namespace llmap::checkpoint {

namespace tools {

ToolResult BashStub(const std::map<std::string, std::string>&) {
    // Future: execute a whitelisted command. For now, refuse safely.
    return {"not_implemented",
            "bash: tool stub -- whitelist not yet defined"};
}

ToolResult RegionLookupStub(const std::map<std::string, std::string>& args) {
    auto it = args.find("region");
    std::string region = it != args.end() ? it->second : "<missing>";
    return {"not_implemented",
            "region_lookup: would query specific_loci/ for " + region};
}

ToolResult PsvCheckStub(const std::map<std::string, std::string>& args) {
    auto it = args.find("region");
    std::string region = it != args.end() ? it->second : "<missing>";
    return {"not_implemented",
            "psv_check: would query loaded PSV catalog for " + region};
}

ToolResult LocalGrepStub(const std::map<std::string, std::string>& args) {
    auto it = args.find("sequence");
    std::string seq = it != args.end() ? it->second : "<missing>";
    return {"not_implemented",
            "local_grep: would grep reference for sequence (len="
            + std::to_string(seq.size()) + ")"};
}

ToolResult WebFetchUcscStub(const std::map<std::string, std::string>& args) {
    auto it = args.find("track");
    std::string track = it != args.end() ? it->second : "<missing>";
    return {"not_implemented",
            "web_fetch_ucsc: would HTTP GET UCSC track " + track};
}

}  // namespace tools

CheckpointToolRegistry::CheckpointToolRegistry() {
    Register({"bash",
              "Execute a whitelisted bash command",
              {"command"},
              tools::BashStub});
    Register({"region_lookup",
              "Look up a coordinate in the specific_loci/ database",
              {"region"},
              tools::RegionLookupStub});
    Register({"psv_check",
              "Query a PSV catalog for paralog-specific variants",
              {"region"},
              tools::PsvCheckStub});
    Register({"local_grep",
              "Grep a chromosome for a sequence",
              {"sequence"},
              tools::LocalGrepStub});
    Register({"web_fetch_ucsc",
              "HTTP GET a UCSC track with caching",
              {"track"},
              tools::WebFetchUcscStub});
}

void CheckpointToolRegistry::Register(ToolDescriptor desc) {
    auto it = tools_.find(desc.name);
    bool was_stub = false;
    if (it != tools_.end()) {
        // crude stub-detection: stubs return "not_implemented" for empty args
        ToolResult probe = it->second.fn({});
        was_stub = probe.status == "not_implemented";
    }
    // probe new tool to update non_stub_count_
    ToolResult new_probe = desc.fn({});
    bool new_is_stub = new_probe.status == "not_implemented";

    if (it != tools_.end()) {
        if (was_stub && !new_is_stub) ++non_stub_count_;
        else if (!was_stub && new_is_stub) --non_stub_count_;
    } else {
        if (!new_is_stub) ++non_stub_count_;
    }
    tools_[desc.name] = std::move(desc);
}

ToolResult CheckpointToolRegistry::Invoke(
    const std::string& name,
    const std::map<std::string, std::string>& args) const {

    auto it = tools_.find(name);
    if (it == tools_.end()) {
        return {"error", "unknown tool: " + name};
    }
    for (const auto& req : it->second.required_args) {
        if (args.find(req) == args.end()) {
            return {"error", "missing required arg '" + req +
                             "' for tool '" + name + "'"};
        }
    }
    return it->second.fn(args);
}

std::vector<ToolDescriptor> CheckpointToolRegistry::List() const {
    std::vector<ToolDescriptor> out;
    out.reserve(tools_.size());
    for (const auto& [_, d] : tools_) out.push_back(d);
    return out;
}

}  // namespace llmap::checkpoint
