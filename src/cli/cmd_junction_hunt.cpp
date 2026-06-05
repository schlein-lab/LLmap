// LLmap — `llmap junction-hunt` CLI command (Mode-5 entry point).
//
// Two operating modes:
//
//  legacy (default):  Builds per-pair k-mer indices from the reference
//                     FASTA on every invocation. Same code path as the
//                     original Mode-5 prototype. Used when --index-dir
//                     is not supplied.
//
//  cached:            Loads pre-built tier_k{N}.bin files from
//                     --index-dir, mmap-maps them, and runs an inverted
//                     read loop (tile each read once, look up the
//                     k-mer hash in the global index, distribute hits
//                     to per-pair accumulators). Memory peak is roughly
//                     proportional to the OS-managed page-cache working
//                     set rather than the panel size. Build the cache
//                     once with `llmap junction-index build`.

#include "cli/commands.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/nahr_pair_loader.h"
#include "junction_hunter/pair_kmer_index.h"
#include "junction_hunter/read_tiler.h"
#include "junction_hunter/consensus_caller.h"
#include "junction_hunter/cascade_pair_index.h"
#include "junction_hunter/cascade_caller.h"
#include "junction_hunter/persistent_index.h"
#include "junction_hunter/reference_gate.h"
#include "io/fasta_reader.h"
#include "io/mmap_fasta.h"

#include <array>

namespace llmap::cli {

namespace {

struct JhArgs {
    std::string pairs_tsv;
    std::string reference;
    std::string reads;
    std::string output;
    std::string index_dir;
    std::string reference_index;   ///< .llmi for routed mode
    std::string nahr_bed;          ///< NAHR-arm BED for routed mode
    std::string routed_out;        ///< optional SKIP/UNMAPPED audit TSV
    int max_pairs{0};
    bool verbose{false};
    bool help{false};
};

void PrintUsage() {
    std::puts(
        "Usage: llmap junction-hunt [options]\n"
        "\n"
        "Detect NAHR breakpoints in reads via multi-k consensus —\n"
        "alignment-free, no minimap2. Outputs one record per (read, pair).\n"
        "\n"
        "Required:\n"
        "  --pairs FILE          NAHR-pair TSV (linter output)\n"
        "  --reads FILE          Input reads (FASTA/FASTQ)\n"
        "  -o, --output FILE     Per-(read,pair) record TSV\n"
        "\n"
        "Mode A (legacy, builds index in-RAM each call):\n"
        "  --reference FILE      Reference FASTA (GRCh38)\n"
        "\n"
        "Mode B (cached, pre-built panel index — much smaller RSS):\n"
        "  --index-dir DIR       Directory holding tier_k*.bin files\n"
        "                        (see `llmap junction-index build`)\n"
        "\n"
        "Mode C (routed, reference-anchored gate + cache fallback):\n"
        "  --reference-index F   .llmi minimizer index of GRCh38/CHM13\n"
        "                        (see `llmap index`); enables routed mode\n"
        "  --nahr-bed FILE       NAHR-arm BED (chrom\\tstart\\tend\\tpair_id\\tarm)\n"
        "  --routed-out FILE     Optional audit TSV for SKIP/UNMAPPED reads\n"
        "  --index-dir DIR       Same persistent cache as Mode B\n"
        "\n"
        "Optional:\n"
        "  --max-pairs N         Limit panel to first N pairs (debug)\n"
        "  -v, --verbose         Verbose progress\n"
        "  -h, --help            This help\n"
    );
}

bool ParseArgs(int argc, char** argv, JhArgs& a) {
    for (int i = 0; i < argc; ++i) {
        std::string x = argv[i];
        auto take = [&](const char* n) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "error: %s requires a value\n", n); return nullptr; }
            return argv[++i];
        };
        if (x == "-h" || x == "--help") { a.help = true; return true; }
        else if (x == "-v" || x == "--verbose") a.verbose = true;
        else if (x == "--pairs")     { auto v = take("--pairs");     if (!v) return false; a.pairs_tsv = v; }
        else if (x == "--reference" || x == "-r") { auto v = take("--reference"); if (!v) return false; a.reference = v; }
        else if (x == "--reads")     { auto v = take("--reads");     if (!v) return false; a.reads = v; }
        else if (x == "-o" || x == "--output") { auto v = take("--output"); if (!v) return false; a.output = v; }
        else if (x == "--index-dir") { auto v = take("--index-dir"); if (!v) return false; a.index_dir = v; }
        else if (x == "--reference-index") { auto v = take("--reference-index"); if (!v) return false; a.reference_index = v; }
        else if (x == "--nahr-bed")  { auto v = take("--nahr-bed");  if (!v) return false; a.nahr_bed = v; }
        else if (x == "--routed-out") { auto v = take("--routed-out"); if (!v) return false; a.routed_out = v; }
        else if (x == "--max-pairs") { auto v = take("--max-pairs"); if (!v) return false; a.max_pairs = std::stoi(v); }
        else { std::fprintf(stderr, "error: unknown arg %s\n", x.c_str()); return false; }
    }
    if (a.pairs_tsv.empty() || a.reads.empty() || a.output.empty()) {
        std::fprintf(stderr, "error: --pairs, --reads, --output are all required\n");
        return false;
    }
    const bool have_routed = !a.reference_index.empty();
    if (have_routed) {
        if (a.index_dir.empty() || a.nahr_bed.empty()) {
            std::fprintf(stderr, "error: --reference-index requires --index-dir AND --nahr-bed\n");
            return false;
        }
    } else if (a.reference.empty() && a.index_dir.empty()) {
        std::fprintf(stderr,
            "error: provide --reference (legacy), --index-dir (cached), or "
            "--reference-index + --index-dir + --nahr-bed (routed)\n");
        return false;
    }
    return true;
}

// ----------------------- Legacy (build-from-scratch) -----------------------

int RunLegacy(const JhArgs& args, std::vector<junction_hunter::NahrPair>& pairs) {
    auto t0 = std::chrono::steady_clock::now();
    llmap::io::MmapFastaReader ref(args.reference);
    if (!ref.IsValid()) {
        std::fprintf(stderr, "error: cannot mmap reference %s: %s\n",
                      args.reference.c_str(), ref.LastError().c_str());
        return 1;
    }

    auto ccfg = junction_hunter::CascadeConfig::LongReadPreset();

    struct PairSeqs { std::string up, down, inn; };
    std::vector<PairSeqs> pair_seqs(pairs.size());
    std::vector<junction_hunter::CascadePairIndex> cidx(pairs.size());
    std::size_t tier1_unique = 0;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto& p = pairs[i];
        const auto up_len   = p.lcr_up_end   > p.lcr_up_start   ? p.lcr_up_end   - p.lcr_up_start   : 0;
        const auto down_len = p.lcr_down_end > p.lcr_down_start ? p.lcr_down_end - p.lcr_down_start : 0;
        const auto in_len   = p.interior_end > p.interior_start ? p.interior_end - p.interior_start : 0;
        pair_seqs[i].up   = ref.GetSubsequence(p.chrom, p.lcr_up_start,   up_len);
        pair_seqs[i].down = ref.GetSubsequence(p.chrom, p.lcr_down_start, down_len);
        pair_seqs[i].inn  = ref.GetSubsequence(p.chrom, p.interior_start, in_len);
        cidx[i].k_values = ccfg.k_values;
        tier1_unique += junction_hunter::BuildCascadeTier1(
            p, pair_seqs[i].up, pair_seqs[i].down, pair_seqs[i].inn, cidx[i]);
        if (args.verbose && (i + 1) % 500 == 0)
            std::fprintf(stderr, "[junction-hunt]   tier-1 built %zu/%zu pairs (%zu k-mers)\n",
                          i + 1, pairs.size(), tier1_unique);
    }
    if (args.verbose)
        std::fprintf(stderr, "[junction-hunt] tier-1 ready: %zu pairs, %zu k-mers (k=%u)\n",
                      pairs.size(), tier1_unique,
                      ccfg.k_values.empty() ? 0u : static_cast<unsigned>(ccfg.k_values[0]));

    std::ofstream out(args.output);
    if (!out) { std::fprintf(stderr, "error: cannot write %s\n", args.output.c_str()); return 1; }
    out << "read_id\tpair_id\tcall\tn_kmer_total\tn_up\tn_dn\tn_in\tn_amb\t"
           "up_mono\tdn_mono\tbreakpoint_read_pos\tbreakpoint_quality\n";

    llmap::io::FastaReader rd(args.reads);
    std::size_t n_reads = 0, n_junctions = 0, n_tier1_pass = 0;
    std::size_t n_promoted_pairs = 0;
    while (rd.HasMore()) {
        auto record = rd.Next();
        if (!record.IsValid()) continue;
        ++n_reads;
        for (std::size_t i = 0; i < cidx.size(); ++i) {
            auto first = junction_hunter::CascadeCall(
                record.name, record.sequence, cidx[i], ccfg);
            if (first.call == junction_hunter::JunctionCall::Unmapped) continue;
            bool any_built_now = false;
            for (std::size_t t = 1; t < cidx[i].k_values.size(); ++t) {
                if (!cidx[i].tier_built[t]) {
                    junction_hunter::BuildCascadeTierN(t,
                        pair_seqs[i].up, pair_seqs[i].down, pair_seqs[i].inn, cidx[i]);
                    any_built_now = true;
                }
            }
            if (any_built_now) ++n_promoted_pairs;
            auto rec = junction_hunter::CascadeCall(
                record.name, record.sequence, cidx[i], ccfg);
            if (rec.call == junction_hunter::JunctionCall::Unmapped) continue;
            ++n_tier1_pass;
            out << rec.read_id << '\t' << rec.pair_id << '\t'
                << junction_hunter::JunctionCallName(rec.call) << '\t'
                << rec.n_kmer_total << '\t'
                << rec.n_consensus_up << '\t' << rec.n_consensus_dn << '\t'
                << rec.n_consensus_in << '\t' << rec.n_ambiguous << '\t'
                << rec.up_monotonicity << '\t' << rec.dn_monotonicity << '\t'
                << rec.breakpoint_read_pos << '\t' << rec.breakpoint_quality << '\n';
            if (rec.call == junction_hunter::JunctionCall::JunctionReal) ++n_junctions;
        }
        if (args.verbose && n_reads % 10000 == 0)
            std::fprintf(stderr,
                "[junction-hunt]   %zu reads, %zu tier-1 passes, %zu promoted pairs, %zu junctions\n",
                n_reads, n_tier1_pass, n_promoted_pairs, n_junctions);
    }
    auto dt = std::chrono::steady_clock::now() - t0;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
    std::fprintf(stderr, "[junction-hunt] LEGACY done in %lds — %zu reads, %zu junction_real calls\n",
                  static_cast<long>(secs), n_reads, n_junctions);
    return 0;
}

// ------------------------- Cached (mmap inverted) --------------------------

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

/// Per-pair accumulator that mirrors the legacy data structures so we
/// can hand off to the existing `CallJunction` classifier unmodified.
struct PerPairAcc {
    /// One unordered_map<hash, KmerLoc> per tier we use (up to 5).
    std::array<junction_hunter::KmerClassMap, 5> per_k_class;
    /// Hit counts at tier-0 (used for promotion gate).
    std::uint32_t tier0_hits{0};
    bool         touched{false};
};

int RunCached(const JhArgs& args, std::vector<junction_hunter::NahrPair>& pairs) {
    auto t0 = std::chrono::steady_clock::now();
    auto ccfg = junction_hunter::CascadeConfig::LongReadPreset();

    // PairKmerIndex / ReadTiling are fixed-array<5> — use the first 5
    // tiers (mirrors what CascadeCall already does via AdaptToFlat).
    const std::size_t n_use_tiers = std::min<std::size_t>(5, ccfg.k_values.size());
    std::array<std::uint8_t, 5> emit_k{};
    for (std::size_t i = 0; i < n_use_tiers; ++i) emit_k[i] = ccfg.k_values[i];

    std::array<junction_hunter::PersistentKmerIndex, 5> tiers;
    for (std::size_t i = 0; i < n_use_tiers; ++i) {
        char p[1024];
        std::snprintf(p, sizeof(p), "%s/tier_k%u.bin",
                      args.index_dir.c_str(), static_cast<unsigned>(emit_k[i]));
        if (!tiers[i].Open(p)) {
            std::fprintf(stderr, "error: open tier %s failed: %s\n",
                          p, tiers[i].LastError().c_str());
            return 1;
        }
        if (args.verbose)
            std::fprintf(stderr, "[junction-hunt] tier k=%u: %llu entries\n",
                          static_cast<unsigned>(emit_k[i]),
                          static_cast<unsigned long long>(tiers[i].NumEntries()));
        tiers[i].HintRandomAccess();
    }

    std::ofstream out(args.output);
    if (!out) { std::fprintf(stderr, "error: cannot write %s\n", args.output.c_str()); return 1; }
    out << "read_id\tpair_id\tcall\tn_kmer_total\tn_up\tn_dn\tn_in\tn_amb\t"
           "up_mono\tdn_mono\tbreakpoint_read_pos\tbreakpoint_quality\n";

    llmap::io::FastaReader rd(args.reads);
    std::size_t n_reads = 0, n_records = 0, n_junctions = 0;

    while (rd.HasMore()) {
        auto record = rd.Next();
        if (!record.IsValid()) continue;
        ++n_reads;
        const std::string_view rseq = record.sequence;

        junction_hunter::ReadTiling tiling;
        tiling.k_values = emit_k;
        std::unordered_map<std::uint32_t, PerPairAcc> acc;
        acc.reserve(128);

        for (std::size_t ti = 0; ti < n_use_tiers; ++ti) {
            const std::uint8_t k = emit_k[ti];
            if (rseq.size() < k) continue;
            for (std::size_t pos = 0; pos + k <= rseq.size(); ++pos) {
                bool clean = true;
                for (std::size_t j = 0; j < k; ++j) {
                    if (!IsAcgt(rseq[pos + j])) { clean = false; break; }
                }
                if (!clean) continue;
                std::uint64_t h = HashKmer(rseq.substr(pos, k));
                tiling.per_k_hashes[ti].push_back(
                    {h, static_cast<std::uint32_t>(pos)});

                auto [first, last] = tiers[ti].Lookup(h);
                for (auto it = first; it != last; ++it) {
                    auto& a = acc[it->pair_id];
                    a.touched = true;
                    a.per_k_class[ti].try_emplace(
                        h, junction_hunter::KmerLoc{
                            static_cast<junction_hunter::LocusClass>(it->cls),
                            it->offset});
                    if (ti == 0) ++a.tier0_hits;
                }
            }
        }

        if (acc.empty()) continue;

        junction_hunter::MultiKConfig mk;
        mk.k_values = emit_k;
        mk.consensus_min = ccfg.consensus_min;
        mk.monotonicity_min = ccfg.monotonicity_min;
        mk.min_psv_switches = 3;

        for (auto& kv : acc) {
            if (kv.second.tier0_hits < ccfg.tier1_min_anchor_hits) continue;
            if (kv.first >= pairs.size()) continue;
            junction_hunter::PairKmerIndex pki;
            pki.pair_id  = pairs[kv.first].pair_id;
            pki.k_values = emit_k;
            for (std::size_t ti = 0; ti < 5; ++ti)
                pki.per_k_class[ti] = std::move(kv.second.per_k_class[ti]);

            auto rec = junction_hunter::CallJunction(
                record.name, tiling, pki, pairs[kv.first], mk);
            if (rec.call == junction_hunter::JunctionCall::Unmapped) continue;
            ++n_records;
            out << rec.read_id << '\t' << rec.pair_id << '\t'
                << junction_hunter::JunctionCallName(rec.call) << '\t'
                << rec.n_kmer_total << '\t'
                << rec.n_consensus_up << '\t' << rec.n_consensus_dn << '\t'
                << rec.n_consensus_in << '\t' << rec.n_ambiguous << '\t'
                << rec.up_monotonicity << '\t' << rec.dn_monotonicity << '\t'
                << rec.breakpoint_read_pos << '\t' << rec.breakpoint_quality << '\n';
            if (rec.call == junction_hunter::JunctionCall::JunctionReal) ++n_junctions;
        }

        if (args.verbose && n_reads % 5000 == 0)
            std::fprintf(stderr,
                "[junction-hunt]   %zu reads, %zu records, %zu junctions\n",
                n_reads, n_records, n_junctions);
    }

    auto dt = std::chrono::steady_clock::now() - t0;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
    std::fprintf(stderr,
        "[junction-hunt] CACHED done in %lds — %zu reads, %zu records, %zu junction_real\n",
        static_cast<long>(secs), n_reads, n_records, n_junctions);
    return 0;
}

// ------------------------------ Routed mode --------------------------------
// Gate: ReferenceGate.Classify(read) → Skip / RouteToNahr / Unmapped.
// On RouteToNahr, run the same inverted lookup as RunCached but ONLY for
// the read's nominated pair_id (if known). For Skip/Unmapped, optionally
// log to a routed-out audit file.

static int RunRouted(JhArgs& args, std::vector<junction_hunter::NahrPair>& pairs) {
    auto t0 = std::chrono::steady_clock::now();
    auto ccfg = junction_hunter::CascadeConfig::LongReadPreset();

    junction_hunter::ReferenceGate gate;
    if (!gate.Load(args.reference_index, args.nahr_bed)) {
        std::fprintf(stderr, "error: reference-gate load: %s\n",
                      gate.LastError().c_str());
        return 1;
    }

    // Open cache tiers (same as cached mode).
    const std::size_t n_use_tiers = std::min<std::size_t>(5, ccfg.k_values.size());
    std::array<std::uint8_t, 5> emit_k{};
    for (std::size_t i = 0; i < n_use_tiers; ++i) emit_k[i] = ccfg.k_values[i];

    std::array<junction_hunter::PersistentKmerIndex, 5> tiers;
    for (std::size_t i = 0; i < n_use_tiers; ++i) {
        char p[1024];
        std::snprintf(p, sizeof(p), "%s/tier_k%u.bin",
                      args.index_dir.c_str(), static_cast<unsigned>(emit_k[i]));
        if (!tiers[i].Open(p)) {
            std::fprintf(stderr, "error: open tier %s failed: %s\n",
                          p, tiers[i].LastError().c_str());
            return 1;
        }
        tiers[i].HintRandomAccess();
    }

    // Pair-id (string) → pairs[] index for fast routing lookup.
    std::unordered_map<std::string, std::uint32_t> pid_to_idx;
    pid_to_idx.reserve(pairs.size());
    for (std::uint32_t i = 0; i < pairs.size(); ++i) {
        pid_to_idx.emplace(pairs[i].pair_id, i);
    }

    std::ofstream out(args.output);
    if (!out) { std::fprintf(stderr, "error: cannot write %s\n", args.output.c_str()); return 1; }
    out << "read_id\tpair_id\tcall\tn_kmer_total\tn_up\tn_dn\tn_in\tn_amb\t"
           "up_mono\tdn_mono\tbreakpoint_read_pos\tbreakpoint_quality\n";

    std::ofstream routed_out;
    if (!args.routed_out.empty()) {
        routed_out.open(args.routed_out);
        if (routed_out)
            routed_out << "read_id\tverdict\tref_id\tref_start\tn_seeds\tpair_id\tarm\n";
    }

    llmap::io::FastaReader rd(args.reads);
    std::size_t n_reads = 0, n_skip = 0, n_route = 0, n_unmap = 0;
    std::size_t n_records = 0, n_junctions = 0;

    while (rd.HasMore()) {
        auto record = rd.Next();
        if (!record.IsValid()) continue;
        ++n_reads;
        const std::string_view rseq = record.sequence;

        auto g = gate.Classify(rseq);
        if (g.verdict == junction_hunter::GateResult::Skip) {
            ++n_skip;
            if (routed_out) routed_out << record.name << "\tSKIP\t" << g.ref_id
                                        << '\t' << g.ref_start << '\t'
                                        << g.n_seeds_in_bucket << "\t\t\n";
            continue;
        }
        if (g.verdict == junction_hunter::GateResult::Unmapped) {
            ++n_unmap;
            if (routed_out) routed_out << record.name << "\tUNMAPPED\t\t\t"
                                        << g.n_seeds_in_bucket << "\t\t\n";
            continue;
        }
        ++n_route;
        if (routed_out) routed_out << record.name << "\tROUTE\t" << g.ref_id
                                    << '\t' << g.ref_start << '\t'
                                    << g.n_seeds_in_bucket << '\t'
                                    << g.pair_id << '\t' << g.arm << '\n';

        // Read is routed to NAHR cache. Build per-read tiling + per-pair
        // accumulator across ALL tiers (mirrors RunCached's body). The
        // route hint (g.pair_id) is logged but we still scan all pairs
        // since a single read may overlap multiple NAHR loci.
        junction_hunter::ReadTiling tiling;
        tiling.k_values = emit_k;
        std::unordered_map<std::uint32_t, PerPairAcc> acc;
        acc.reserve(32);

        for (std::size_t ti = 0; ti < n_use_tiers; ++ti) {
            const std::uint8_t k = emit_k[ti];
            if (rseq.size() < k) continue;
            for (std::size_t pos = 0; pos + k <= rseq.size(); ++pos) {
                bool clean = true;
                for (std::size_t j = 0; j < k; ++j) {
                    if (!IsAcgt(rseq[pos + j])) { clean = false; break; }
                }
                if (!clean) continue;
                std::uint64_t h = HashKmer(rseq.substr(pos, k));
                tiling.per_k_hashes[ti].push_back(
                    {h, static_cast<std::uint32_t>(pos)});

                auto [first, last] = tiers[ti].Lookup(h);
                for (auto it = first; it != last; ++it) {
                    auto& a = acc[it->pair_id];
                    a.touched = true;
                    a.per_k_class[ti].try_emplace(
                        h, junction_hunter::KmerLoc{
                            static_cast<junction_hunter::LocusClass>(it->cls),
                            it->offset});
                    if (ti == 0) ++a.tier0_hits;
                }
            }
        }

        junction_hunter::MultiKConfig mk;
        mk.k_values = emit_k;
        mk.consensus_min = ccfg.consensus_min;
        mk.monotonicity_min = ccfg.monotonicity_min;
        mk.min_psv_switches = 3;

        for (auto& kv : acc) {
            if (kv.second.tier0_hits < ccfg.tier1_min_anchor_hits) continue;
            if (kv.first >= pairs.size()) continue;
            junction_hunter::PairKmerIndex pki;
            pki.pair_id  = pairs[kv.first].pair_id;
            pki.k_values = emit_k;
            for (std::size_t ti = 0; ti < 5; ++ti)
                pki.per_k_class[ti] = std::move(kv.second.per_k_class[ti]);
            auto rec = junction_hunter::CallJunction(
                record.name, tiling, pki, pairs[kv.first], mk);
            if (rec.call == junction_hunter::JunctionCall::Unmapped) continue;
            ++n_records;
            out << rec.read_id << '\t' << rec.pair_id << '\t'
                << junction_hunter::JunctionCallName(rec.call) << '\t'
                << rec.n_kmer_total << '\t'
                << rec.n_consensus_up << '\t' << rec.n_consensus_dn << '\t'
                << rec.n_consensus_in << '\t' << rec.n_ambiguous << '\t'
                << rec.up_monotonicity << '\t' << rec.dn_monotonicity << '\t'
                << rec.breakpoint_read_pos << '\t' << rec.breakpoint_quality << '\n';
            if (rec.call == junction_hunter::JunctionCall::JunctionReal) ++n_junctions;
        }

        if (args.verbose && n_reads % 1000 == 0)
            std::fprintf(stderr,
                "[junction-hunt]   %zu reads | skip=%zu route=%zu unmap=%zu | records=%zu jct=%zu\n",
                n_reads, n_skip, n_route, n_unmap, n_records, n_junctions);
    }

    auto dt = std::chrono::steady_clock::now() - t0;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(dt).count();
    std::fprintf(stderr,
        "[junction-hunt] ROUTED done in %lds — %zu reads | skip=%zu (%.1f%%) "
        "route=%zu (%.2f%%) unmap=%zu (%.1f%%) | records=%zu jct=%zu\n",
        static_cast<long>(secs), n_reads,
        n_skip,  n_reads ? 100.0 * n_skip  / n_reads : 0.0,
        n_route, n_reads ? 100.0 * n_route / n_reads : 0.0,
        n_unmap, n_reads ? 100.0 * n_unmap / n_reads : 0.0,
        n_records, n_junctions);
    return 0;
}

}  // namespace

int run_junction_hunt(int argc, char** argv) {
    JhArgs args;
    if (!ParseArgs(argc, argv, args)) { PrintUsage(); return 2; }
    if (args.help) { PrintUsage(); return 0; }

    std::vector<junction_hunter::NahrPair> pairs;
    auto pst = junction_hunter::LoadNahrPairsTsv(args.pairs_tsv, pairs);
    if (!pst.ok) {
        std::fprintf(stderr, "error: pair-tsv load failed: %s\n", pst.error.c_str());
        return 1;
    }
    if (args.max_pairs > 0 && static_cast<int>(pairs.size()) > args.max_pairs)
        pairs.resize(args.max_pairs);
    if (args.verbose)
        std::fprintf(stderr, "[junction-hunt] %zu pairs loaded\n", pairs.size());

    if (!args.reference_index.empty()) return RunRouted(args, pairs);
    if (!args.index_dir.empty())       return RunCached(args, pairs);
    return RunLegacy(args, pairs);
}

}  // namespace llmap::cli
