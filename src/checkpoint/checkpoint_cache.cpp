// LLmap -- file-backed checkpoint cache implementation.
//
// Includes:
//   - a small self-contained SHA-256 (FIPS 180-4) so we don't depend on
//     OpenSSL just for the cache key.
//   - a small JSON writer/parser specialised for AgentDecision; no external
//     JSON dependency is needed at the checkpoint layer.

#include "checkpoint/checkpoint_cache.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <system_error>

namespace llmap::checkpoint {

namespace {

// ---- SHA-256 (FIPS 180-4) ---------------------------------------------

constexpr std::array<uint32_t, 64> kSha256K = {
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,
    0x923f82a4u,0xab1c5ed5u,0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
    0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,0xe49b69c1u,0xefbe4786u,
    0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,
    0x06ca6351u,0x14292967u,0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
    0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,0xa2bfe8a1u,0xa81a664bu,
    0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,
    0x5b9cca4fu,0x682e6ff3u,0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
    0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

inline uint32_t RotR(uint32_t x, uint32_t n) { return (x >> n) | (x << (32 - n)); }

std::string Sha256Hex(const std::string& data) {
    uint32_t h[8] = {0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
                     0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u};

    // Pad: data || 0x80 || 0x00s || 64-bit length.
    std::vector<uint8_t> msg(data.begin(), data.end());
    const uint64_t bit_len = static_cast<uint64_t>(data.size()) * 8;
    msg.push_back(0x80u);
    while (msg.size() % 64 != 56) msg.push_back(0x00u);
    for (int i = 7; i >= 0; --i) {
        msg.push_back(static_cast<uint8_t>((bit_len >> (i * 8)) & 0xFF));
    }

    for (size_t base = 0; base < msg.size(); base += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (uint32_t(msg[base + i*4]) << 24) |
                   (uint32_t(msg[base + i*4 + 1]) << 16) |
                   (uint32_t(msg[base + i*4 + 2]) << 8)  |
                   (uint32_t(msg[base + i*4 + 3]));
        }
        for (int i = 16; i < 64; ++i) {
            uint32_t s0 = RotR(w[i-15], 7) ^ RotR(w[i-15], 18) ^ (w[i-15] >> 3);
            uint32_t s1 = RotR(w[i-2], 17) ^ RotR(w[i-2], 19) ^ (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }
        uint32_t a=h[0], b=h[1], c=h[2], d=h[3], e=h[4], f=h[5], g=h[6], hh=h[7];
        for (int i = 0; i < 64; ++i) {
            uint32_t S1 = RotR(e, 6) ^ RotR(e, 11) ^ RotR(e, 25);
            uint32_t ch = (e & f) ^ (~e & g);
            uint32_t t1 = hh + S1 + ch + kSha256K[i] + w[i];
            uint32_t S0 = RotR(a, 2) ^ RotR(a, 13) ^ RotR(a, 22);
            uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = S0 + mj;
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    std::ostringstream oss;
    oss << std::hex;
    for (uint32_t v : h) {
        char buf[9];
        std::snprintf(buf, sizeof(buf), "%08x", v);
        oss << buf;
    }
    return oss.str();
}

// ---- Default cache root -----------------------------------------------

std::filesystem::path DefaultRootDir() {
    const char* xdg = std::getenv("XDG_CACHE_HOME");
    if (xdg && *xdg) {
        return std::filesystem::path(xdg) / "llmap" / "checkpoints";
    }
    const char* home = std::getenv("HOME");
    if (home && *home) {
        return std::filesystem::path(home) / ".cache" / "llmap" / "checkpoints";
    }
    return std::filesystem::temp_directory_path() / "llmap_checkpoints";
}

// ---- Tiny JSON helpers -----------------------------------------------

void EscapeString(std::ostringstream& oss, const std::string& s) {
    oss << '"';
    for (char c : s) {
        switch (c) {
            case '"':  oss << "\\\""; break;
            case '\\': oss << "\\\\"; break;
            case '\n': oss << "\\n";  break;
            case '\r': oss << "\\r";  break;
            case '\t': oss << "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    oss << buf;
                } else {
                    oss << c;
                }
        }
    }
    oss << '"';
}

// Look up `"key":<value>` in a flat JSON blob and return the substring of
// <value>. Caller decides how to interpret it. Returns empty optional on
// "not found". This is intentionally tolerant: it does not validate the
// surrounding JSON.
std::optional<std::string> FindRawField(const std::string& json,
                                        const std::string& key) {
    std::string needle = "\"" + key + "\"";
    size_t p = json.find(needle);
    if (p == std::string::npos) return std::nullopt;
    p = json.find(':', p);
    if (p == std::string::npos) return std::nullopt;
    ++p;
    while (p < json.size() && (json[p] == ' ' || json[p] == '\t')) ++p;
    if (p >= json.size()) return std::nullopt;

    if (json[p] == '"') {
        // string
        size_t q = p + 1;
        std::string out;
        while (q < json.size() && json[q] != '"') {
            if (json[q] == '\\' && q + 1 < json.size()) {
                char c = json[q + 1];
                if (c == 'n') out += '\n';
                else if (c == 't') out += '\t';
                else if (c == 'r') out += '\r';
                else out += c;
                q += 2;
            } else {
                out += json[q++];
            }
        }
        return out;
    }

    // number, bool, null, array, or object: read until matching delimiter
    int depth = 0;
    size_t q = p;
    while (q < json.size()) {
        char c = json[q];
        if (c == '[' || c == '{') ++depth;
        else if (c == ']' || c == '}') {
            if (depth == 0) break;
            --depth;
        } else if ((c == ',' || c == '\n') && depth == 0) {
            break;
        }
        ++q;
    }
    return json.substr(p, q - p);
}

}  // namespace

// ---- MakeCacheKey -----------------------------------------------------

std::string MakeCacheKey(CheckpointType type, const CheckpointContext& ctx) {
    auto candidates = ctx.candidate_positions;
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });

    std::ostringstream oss;
    oss << "v" << kSchemaVersion
        << '|' << static_cast<int>(type)
        << '|' << ctx.read_id
        << '|' << ctx.ref_id
        << '|' << ctx.region_name
        << '|';
    for (const auto& [pos, score] : candidates) {
        oss << pos << ':' << score << ';';
    }
    return Sha256Hex(oss.str());
}

// ---- Serialise / parse AgentDecision ---------------------------------

std::string SerializeDecision(const AgentDecision& d) {
    std::ostringstream oss;
    oss << "{\"schema_version\":" << kSchemaVersion
        << ",\"consulted\":" << (d.consulted ? "true" : "false")
        << ",\"fallback_used\":" << (d.fallback_used ? "true" : "false")
        << ",\"reasoning\":";
    EscapeString(oss, d.reasoning);
    oss << ",\"special_finding\":";
    EscapeString(oss, d.special_finding);

    oss << ",\"wave\":[";
    for (size_t i = 0; i < d.wave.size(); ++i) {
        if (i) oss << ',';
        oss << "[" << d.wave[i].first << "," << d.wave[i].second << "]";
    }
    oss << "]";

    // ParamOverride: write only the fields that are set, as floats/ints
    oss << ",\"override\":{";
    bool first = true;
    auto emit = [&](const char* k, auto v) {
        if (!first) oss << ',';
        first = false;
        oss << "\"" << k << "\":" << v;
    };
    if (d.override.k)                emit("k", static_cast<int>(*d.override.k));
    if (d.override.w)                emit("w", static_cast<int>(*d.override.w));
    if (d.override.max_occ)          emit("max_occ", *d.override.max_occ);
    if (d.override.lambda_scale)     emit("lambda_scale", *d.override.lambda_scale);
    if (d.override.identity_threshold) emit("identity_threshold", *d.override.identity_threshold);
    if (d.override.anchor_weight_scale) emit("anchor_weight_scale", *d.override.anchor_weight_scale);
    if (d.override.report_multi_position) emit("report_multi_position", *d.override.report_multi_position ? "true":"false");
    if (d.override.require_psv_disambig)  emit("require_psv_disambig",  *d.override.require_psv_disambig ? "true":"false");
    if (d.override.allow_high_mismatch)   emit("allow_high_mismatch",   *d.override.allow_high_mismatch ? "true":"false");
    if (d.override.require_llm_at_runtime) emit("require_llm_at_runtime", *d.override.require_llm_at_runtime ? "true":"false");
    oss << "}}";
    return oss.str();
}

std::optional<AgentDecision> DeserializeDecision(const std::string& json) {
    auto schema = FindRawField(json, "schema_version");
    if (!schema || std::atoi(schema->c_str()) != kSchemaVersion) {
        return std::nullopt;
    }
    AgentDecision d;
    if (auto v = FindRawField(json, "consulted"))      d.consulted = (*v == "true");
    if (auto v = FindRawField(json, "fallback_used"))  d.fallback_used = (*v == "true");
    if (auto v = FindRawField(json, "reasoning"))      d.reasoning = *v;
    if (auto v = FindRawField(json, "special_finding")) d.special_finding = *v;

    if (auto raw = FindRawField(json, "wave")) {
        // Parse "[[pos,amp],[pos,amp],...]"
        const std::string& s = *raw;
        size_t i = 0;
        while (i < s.size()) {
            while (i < s.size() && s[i] != '[') ++i;
            if (i >= s.size()) break;
            ++i;
            char* end = nullptr;
            uint32_t pos = static_cast<uint32_t>(std::strtoul(s.c_str() + i, &end, 10));
            if (!end || end == s.c_str() + i) break;
            i = end - s.c_str();
            while (i < s.size() && s[i] == ',') ++i;
            float amp = std::strtof(s.c_str() + i, &end);
            if (!end) break;
            i = end - s.c_str();
            while (i < s.size() && s[i] != ']') ++i;
            d.wave.emplace_back(pos, amp);
            ++i;
        }
    }
    return d;
}

// ---- CheckpointCache --------------------------------------------------

CheckpointCache::CheckpointCache(std::filesystem::path root_dir)
    : root_dir_(root_dir.empty() ? DefaultRootDir() : std::move(root_dir)) {}

std::filesystem::path CheckpointCache::PathForKey(const std::string& key) const {
    std::string shard = key.size() >= 2 ? key.substr(0, 2) : "xx";
    return root_dir_ / shard / (key + ".json");
}

std::optional<AgentDecision> CheckpointCache::Lookup(const std::string& key) const {
    auto path = PathForKey(key);
    std::error_code ec;
    if (!std::filesystem::exists(path, ec) || ec) return std::nullopt;

    std::ifstream in(path);
    if (!in.good()) return std::nullopt;
    std::ostringstream buf;
    buf << in.rdbuf();
    return DeserializeDecision(buf.str());
}

bool CheckpointCache::Store(const std::string& key, const AgentDecision& decision) {
    auto path = PathForKey(key);
    std::error_code ec;
    std::filesystem::create_directories(path.parent_path(), ec);
    if (ec) return false;

    auto tmp = path;
    tmp += ".tmp";
    {
        std::ofstream out(tmp);
        if (!out.good()) return false;
        out << SerializeDecision(decision);
    }
    std::filesystem::rename(tmp, path, ec);
    return !ec;
}

}  // namespace llmap::checkpoint
