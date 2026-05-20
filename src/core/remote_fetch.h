// LLmap — remote reference/anchor fetching.
//
// Lets any reader accept an http(s):// or s3:// URL instead of a local path.
// The file is fetched once into a local cache and the cached path is returned,
// so existing FASTA/anchor readers work unchanged on remote inputs. Designed
// for HPC compute nodes that have outbound internet (verified on Hummel) and
// for the planned LLM-in-LLmap consult layer.
//
// Implementation is dependency-light: it shells out (via fork/exec, no shell,
// no injection) to `curl` then `wget`. No compile-time libcurl linkage, so the
// feature degrades gracefully — if neither downloader exists, FetchToCache
// returns nullopt and the caller reports a normal "could not open" error.

#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace llmap::core {

// True for http://, https://, or s3:// URLs.
[[nodiscard]] bool IsRemotePath(std::string_view path);

// Translate s3://bucket/key to https://bucket.s3.amazonaws.com/key.
// Non-s3 URLs are returned unchanged.
[[nodiscard]] std::string S3ToHttps(std::string_view url);

// Fetch `url` into a local cache and return the cached file path. Idempotent:
// a non-empty cached copy is reused. `cache_dir` defaults to
// $LLMAP_REMOTE_CACHE, else $TMPDIR/llmap_remote_cache, else /tmp/...
// If `also_index` is true, sibling .fai and .gzi are fetched too (best-effort)
// so the cached file can back samtools-style region queries.
// Returns nullopt on download failure or if no downloader is available.
[[nodiscard]] std::optional<std::string> FetchToCache(
    std::string_view url, std::string_view cache_dir = "",
    bool also_index = false);

}  // namespace llmap::core
