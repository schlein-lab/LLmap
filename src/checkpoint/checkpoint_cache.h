// LLmap -- file-backed cache for agent checkpoint decisions.
//
// The cache key is the SHA-256 of (schema_version, checkpoint type, read_id,
// region_name, candidate_positions sorted by (pos, score)). The cache value
// is the JSON-serialised AgentDecision. Cache files are sharded by the first
// byte of the hash to avoid huge flat directories:
//
//   ~/.cache/llmap/checkpoints/<hash[0:2]>/<hash>.json
//
// Bumping kSchemaVersion invalidates the whole cache (old hashes stop
// matching, eventually swept by maintenance). A missing cache file is not
// an error; we just report std::nullopt.

#pragma once

#include "checkpoint/checkpoint_types.h"

#include <filesystem>
#include <optional>
#include <string>

namespace llmap::checkpoint {

// Bump this whenever AgentDecision serialisation changes shape.
constexpr int kSchemaVersion = 1;

// Hex SHA-256 over a deterministic flattening of context. Pure function;
// no IO. The schema_version is folded in so any format change invalidates
// previous keys without manual cleanup.
std::string MakeCacheKey(CheckpointType type, const CheckpointContext& ctx);

// JSON serialisation. We use a tiny hand-rolled writer/parser; no nlohmann.
std::string SerializeDecision(const AgentDecision& decision);
std::optional<AgentDecision> DeserializeDecision(const std::string& json);

class CheckpointCache {
public:
    // root_dir defaults to ~/.cache/llmap/checkpoints if empty.
    explicit CheckpointCache(std::filesystem::path root_dir = {});

    // Read a previously-stored decision. nullopt on miss / unreadable file.
    std::optional<AgentDecision> Lookup(const std::string& key) const;

    // Atomic-ish write: writes to tmp + renames. Best effort; returns false
    // on filesystem error but does not throw.
    bool Store(const std::string& key, const AgentDecision& decision);

    // Where the cache lives, e.g. for tests / inspection.
    const std::filesystem::path& RootDir() const { return root_dir_; }

private:
    std::filesystem::path PathForKey(const std::string& key) const;

    std::filesystem::path root_dir_;
};

}  // namespace llmap::checkpoint
