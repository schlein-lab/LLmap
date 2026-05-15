// LLmap -- tools the active-mapping agent can call.
//
// Each tool is a stub for now: it advertises its name + schema and returns
// "not yet implemented" when invoked. What matters at this stage is the
// REGISTRATION mechanism: future code can wire up actual tool execution
// without changing call sites.
//
// The dispatcher passes a CheckpointToolRegistry to the agent session and
// the agent invokes tools through it. Tool names mirror the names used in
// SystemPrompt() in checkpoint_prompts.cpp.

#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace llmap::checkpoint {

// One tool call result. status="ok" or "error" or "not_implemented".
struct ToolResult {
    std::string status;
    std::string output;
};

// Tool signature: receives a flat key/value argument map, returns ToolResult.
using ToolFn = std::function<ToolResult(const std::map<std::string, std::string>&)>;

// A registered tool: name, JSON-schema-ish description, callable.
struct ToolDescriptor {
    std::string name;
    std::string description;
    std::vector<std::string> required_args;
    ToolFn fn;
};

class CheckpointToolRegistry {
public:
    CheckpointToolRegistry();

    // Register a tool, replacing any existing one with the same name.
    void Register(ToolDescriptor desc);

    // Invoke a tool by name. Returns status="error" if the tool is unknown
    // and status="not_implemented" if the tool's stub has not been replaced.
    ToolResult Invoke(const std::string& name,
                      const std::map<std::string, std::string>& args) const;

    // Tool inventory, e.g. for emitting a tool catalog to the agent.
    std::vector<ToolDescriptor> List() const;

    // Convenience: are any non-stub tools registered?
    bool HasNonStubTools() const { return non_stub_count_ > 0; }

private:
    std::map<std::string, ToolDescriptor> tools_;
    size_t non_stub_count_ = 0;
};

// Pre-defined stub tools used at construction time. Exposed for tests that
// want to verify the inventory without instantiating a full agent.
namespace tools {
ToolResult BashStub(const std::map<std::string, std::string>& args);
ToolResult RegionLookupStub(const std::map<std::string, std::string>& args);
ToolResult PsvCheckStub(const std::map<std::string, std::string>& args);
ToolResult LocalGrepStub(const std::map<std::string, std::string>& args);
ToolResult WebFetchUcscStub(const std::map<std::string, std::string>& args);
}  // namespace tools

}  // namespace llmap::checkpoint
