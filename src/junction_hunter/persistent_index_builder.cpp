// LLmap — junction_hunter: persistent k-mer index builder (impl).

#include "junction_hunter/persistent_index_builder.h"

#include "io/mmap_fasta.h"
#include "junction_hunter/pair_kmer_index.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace llmap::junction_hunter {

namespace {

inline char Comp(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}
inline bool IsAcgt(char c) noexcept {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T'
        || c == 'a' || c == 'c' || c == 'g' || c == 't';
}

/// Canonical FNV-1a hash matching HashKmer in cascade_pair_index.cpp.
/// MUST stay byte-identical to the runtime hash or the cache breaks.
std::uint64_t HashKmer(std::string_view k) noexcept {
    std::uint64_t hf = 1469598103934665603ULL;
    std::uint64_t hr = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    for (std::size_t i = 0; i < k.size(); ++i) {
        char cf = k[i];
        char cr = Comp(k[k.size() - 1 - i]);
        hf = (hf ^ static_cast<std::uint64_t>(cf)) * fnv_prime;
        hr = (hr ^ static_cast<std::uint64_t>(cr)) * fnv_prime;
    }
    return hf < hr ? hf : hr;
}

/// Cheap content-stamp: FNV-1a over file bytes + length, written into
/// the 32-byte sha256 field as 8 raw bytes (rest zero-padded). Not
/// cryptographic — only used to detect "cache built against a different
/// panel/reference than now in use". Two distinct files giving the same
/// 64-bit stamp is astronomically unlikely for our inputs.
bool ComputeContentStamp(const std::string& path, std::uint8_t out[32], std::string& err) {
    std::memset(out, 0, 32);
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open for stamp: " + path; return false; }
    constexpr std::uint64_t fnv_init  = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime = 1099511628211ULL;
    std::uint64_t h = fnv_init;
    std::uint64_t total = 0;
    char buf[1 << 16];
    while (f) {
        f.read(buf, sizeof(buf));
        std::streamsize got = f.gcount();
        for (std::streamsize i = 0; i < got; ++i) {
            h = (h ^ static_cast<std::uint8_t>(buf[i])) * fnv_prime;
        }
        total += static_cast<std::uint64_t>(got);
    }
    // Mix file length so length-changes are guaranteed to alter the stamp.
    h = (h ^ total) * fnv_prime;
    std::memcpy(out, &h, sizeof(h));
    std::memcpy(out + 8, &total, sizeof(total));
    return true;
}

void EnumerateKmers(std::string_view seq,
                    std::uint8_t k,
                    LocusClass cls,
                    std::uint32_t pair_id,
                    std::vector<PersistentIndexEntry>& sink) {
    if (seq.size() < k) return;
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        bool clean = true;
        for (std::size_t j = 0; j < k; ++j) {
            if (!IsAcgt(seq[i + j])) { clean = false; break; }
        }
        if (!clean) continue;
        PersistentIndexEntry e{};
        e.hash    = HashKmer(seq.substr(i, k));
        e.pair_id = pair_id;
        e.offset  = static_cast<std::uint32_t>(i);
        e.cls     = static_cast<std::uint8_t>(cls);
        sink.push_back(e);
    }
}

/// After sorting by (hash, pair_id), collapse entries with the SAME
/// (hash, pair_id):
///   - cls differs across collapsed entries → upgrade to Ambiguous
///   - cls identical → keep first (lowest offset)
/// Result mirrors the runtime KmerClassMap behaviour byte-for-byte.
std::uint64_t CollapseInPlace(std::vector<PersistentIndexEntry>& v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end(),
              [](const PersistentIndexEntry& a, const PersistentIndexEntry& b) {
                  if (a.hash != b.hash) return a.hash < b.hash;
                  if (a.pair_id != b.pair_id) return a.pair_id < b.pair_id;
                  return a.offset < b.offset;
              });
    std::size_t w = 0;
    for (std::size_t r = 0; r < v.size(); ) {
        std::size_t s = r;
        while (s < v.size()
               && v[s].hash == v[r].hash
               && v[s].pair_id == v[r].pair_id) {
            ++s;
        }
        PersistentIndexEntry merged = v[r];
        for (std::size_t k = r + 1; k < s; ++k) {
            if (v[k].cls != merged.cls
                && v[k].cls != static_cast<std::uint8_t>(LocusClass::Ambiguous)) {
                merged.cls = static_cast<std::uint8_t>(LocusClass::Ambiguous);
            }
        }
        v[w++] = merged;
        r = s;
    }
    v.resize(w);
    return v.size();
}

bool WriteTier(const std::string& path,
               std::uint8_t k,
               std::uint64_t n_pairs,
               const std::uint8_t panel_sha[32],
               const std::uint8_t ref_sha[32],
               const std::vector<PersistentIndexEntry>& entries,
               std::string& err) {
    PersistentIndexHeader hdr{};
    std::memcpy(hdr.magic, kPersistentIndexMagic, 16);
    hdr.version = kPersistentIndexVersion;
    hdr.k = k;
    hdr.n_entries = entries.size();
    hdr.n_pairs   = n_pairs;
    std::memcpy(hdr.panel_sha256, panel_sha, 32);
    std::memcpy(hdr.ref_sha256,   ref_sha,   32);
    {
        std::time_t t = std::time(nullptr);
        std::tm tm{};
        gmtime_r(&t, &tm);
        std::strftime(hdr.built_iso, sizeof(hdr.built_iso),
                      "%Y-%m-%dT%H:%M:%SZ", &tm);
    }
    {
        char host[64] = {0};
        ::gethostname(host, sizeof(host) - 1);
        std::strncpy(hdr.built_host, host, sizeof(hdr.built_host) - 1);
    }

    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) { err = "open for write: " + path; return false; }
    if (std::fwrite(&hdr, sizeof(hdr), 1, fp) != 1) {
        err = "write header"; std::fclose(fp); return false;
    }
    if (!entries.empty()) {
        if (std::fwrite(entries.data(), sizeof(PersistentIndexEntry),
                        entries.size(), fp) != entries.size()) {
            err = "write entries"; std::fclose(fp); return false;
        }
    }
    if (std::fclose(fp) != 0) { err = "close: " + path; return false; }
    return true;
}

}  // namespace

bool BuildPersistentIndex(const std::vector<NahrPair>& pairs,
                          const std::string& reference_path,
                          const BuildOptions& opts,
                          BuildStats& stats,
                          std::string& err) {
    auto t_total = std::chrono::steady_clock::now();
    err.clear();

    std::filesystem::create_directories(opts.out_dir);

    std::uint8_t panel_stamp[32] = {0};
    std::uint8_t ref_stamp[32]   = {0};
    if (!ComputeContentStamp(opts.panel_path, panel_stamp, err)) return false;
    if (!ComputeContentStamp(reference_path,  ref_stamp,   err)) return false;

    llmap::io::MmapFastaReader ref(reference_path);
    if (!ref.IsValid()) {
        err = "mmap reference failed: " + ref.LastError();
        return false;
    }

    struct PairSeqs { std::string up, down, inn; };
    std::vector<PairSeqs> pair_seqs(pairs.size());
    if (opts.verbose) std::fprintf(stderr, "[build] extracting %zu pair sequences\n", pairs.size());
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& p = pairs[i];
        const auto up_len   = p.lcr_up_end   > p.lcr_up_start   ? p.lcr_up_end   - p.lcr_up_start   : 0;
        const auto down_len = p.lcr_down_end > p.lcr_down_start ? p.lcr_down_end - p.lcr_down_start : 0;
        const auto in_len   = p.interior_end > p.interior_start ? p.interior_end - p.interior_start : 0;
        pair_seqs[i].up   = ref.GetSubsequence(p.chrom, p.lcr_up_start,   up_len);
        pair_seqs[i].down = ref.GetSubsequence(p.chrom, p.lcr_down_start, down_len);
        pair_seqs[i].inn  = ref.GetSubsequence(p.chrom, p.interior_start, in_len);
    }

    stats.entries_per_tier.assign(opts.k_values.size(), 0);
    stats.bytes_per_tier.assign(opts.k_values.size(), 0);

    for (std::size_t t = 0; t < opts.k_values.size(); ++t) {
        auto t_tier = std::chrono::steady_clock::now();
        const std::uint8_t k = opts.k_values[t];
        if (opts.verbose) std::fprintf(stderr, "[build] tier %zu (k=%u): enumerating\n", t, k);

        std::vector<PersistentIndexEntry> entries;
        for (std::uint32_t pid = 0; pid < pairs.size(); ++pid) {
            EnumerateKmers(pair_seqs[pid].up,   k, LocusClass::LcrUp,    pid, entries);
            EnumerateKmers(pair_seqs[pid].down, k, LocusClass::LcrDown,  pid, entries);
            EnumerateKmers(pair_seqs[pid].inn,  k, LocusClass::Interior, pid, entries);
        }
        if (opts.verbose)
            std::fprintf(stderr, "[build]   raw=%zu (%.1f GB), sorting+collapsing\n",
                          entries.size(),
                          entries.size() * sizeof(PersistentIndexEntry) / 1e9);

        std::uint64_t n_final = CollapseInPlace(entries);

        char tier_path[1024];
        std::snprintf(tier_path, sizeof(tier_path), "%s/tier_k%u.bin",
                      opts.out_dir.c_str(), static_cast<unsigned>(k));
        if (!WriteTier(tier_path, k, pairs.size(), panel_stamp, ref_stamp, entries, err)) {
            return false;
        }
        stats.entries_per_tier[t] = n_final;
        stats.bytes_per_tier[t]   = sizeof(PersistentIndexHeader)
                                  + n_final * sizeof(PersistentIndexEntry);

        // Free the entry vector before next tier.
        std::vector<PersistentIndexEntry>().swap(entries);

        auto dt = std::chrono::steady_clock::now() - t_tier;
        long secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
        if (opts.verbose)
            std::fprintf(stderr, "[build]   tier %u: %llu entries, %.2f GB, %lds\n",
                          static_cast<unsigned>(k),
                          static_cast<unsigned long long>(n_final),
                          stats.bytes_per_tier[t] / 1e9,
                          secs);
    }

    {
        std::string mpath = opts.out_dir + "/MANIFEST.tsv";
        std::ofstream m(mpath);
        if (!m) { err = "open manifest"; return false; }
        m << "k\tn_entries\tbytes\n";
        for (std::size_t t = 0; t < opts.k_values.size(); ++t) {
            m << static_cast<unsigned>(opts.k_values[t]) << '\t'
              << stats.entries_per_tier[t] << '\t'
              << stats.bytes_per_tier[t] << '\n';
        }
    }
    auto write_hex = [](const std::string& path, const std::uint8_t s[32]) {
        std::ofstream f(path);
        for (int i = 0; i < 32; ++i) {
            char buf[3]; std::snprintf(buf, sizeof(buf), "%02x", s[i]);
            f << buf;
        }
        f << '\n';
    };
    write_hex(opts.out_dir + "/panel.stamp",     panel_stamp);
    write_hex(opts.out_dir + "/reference.stamp", ref_stamp);

    auto dt_total = std::chrono::steady_clock::now() - t_total;
    stats.seconds_total =
        std::chrono::duration_cast<std::chrono::duration<double>>(dt_total).count();
    return true;
}

}  // namespace llmap::junction_hunter
