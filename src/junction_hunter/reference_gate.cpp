// LLmap — junction_hunter: reference-anchored routing gate impl.

#include "junction_hunter/reference_gate.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace llmap::junction_hunter {

namespace {

constexpr std::uint32_t kNoRefId = std::numeric_limits<std::uint32_t>::max();

// Parse a 5-column BED into AnnotationIntervals. Columns: chrom, start,
// end, pair_id, arm. Maps chrom→ref_id via the supplied map; rows with
// an unknown chrom are dropped (counted in `dropped`).
bool ParseNahrBed(const std::string& path,
                  const std::unordered_map<std::string, std::uint32_t>& chrom_to_id,
                  std::vector<annot::AnnotationInterval>& out,
                  std::size_t& dropped,
                  std::string& err) {
    std::ifstream fh(path);
    if (!fh) { err = "open BED failed: " + path; return false; }
    out.clear();
    dropped = 0;
    std::string line;
    while (std::getline(fh, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Split on '\t'
        std::array<std::string, 5> col;
        std::size_t pos = 0, idx = 0;
        for (std::size_t i = 0; i <= line.size() && idx < 5; ++i) {
            if (i == line.size() || line[i] == '\t') {
                col[idx++] = line.substr(pos, i - pos);
                pos = i + 1;
            }
        }
        if (idx < 5) continue;
        auto it = chrom_to_id.find(col[0]);
        if (it == chrom_to_id.end()) { ++dropped; continue; }
        annot::AnnotationInterval iv;
        iv.ref_id      = it->second;
        try {
            iv.start = static_cast<std::uint32_t>(std::stoul(col[1]));
            iv.end   = static_cast<std::uint32_t>(std::stoul(col[2]));
        } catch (...) { ++dropped; continue; }
        iv.region_name = col[3] + "|" + col[4];  // pair_id|arm
        iv.source      = "nahr_pair_arm";
        iv.layer       = annot::AnnotationLayer::SpecificLocus;
        out.push_back(std::move(iv));
    }
    return true;
}

}  // namespace

bool ReferenceGate::Load(const std::string& llmi_path,
                         const std::string& nahr_bed_path) {
    last_error_.clear();
    idx_ = classical::MinimizerIndex::Load(llmi_path);
    if (!idx_) {
        last_error_ = "MinimizerIndex::Load failed: " + llmi_path;
        return false;
    }

    chrom_to_ref_id_.clear();
    const auto& seqs = idx_->GetSequences();
    chrom_to_ref_id_.reserve(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        chrom_to_ref_id_.emplace(seqs[i].name,
                                  static_cast<std::uint32_t>(i));
    }

    std::vector<annot::AnnotationInterval> ivs;
    std::size_t dropped = 0;
    if (!ParseNahrBed(nahr_bed_path, chrom_to_ref_id_, ivs, dropped, last_error_)) {
        return false;
    }
    if (ivs.empty()) {
        last_error_ = "no NAHR intervals parsed from " + nahr_bed_path;
        return false;
    }
    nahr_iv_.Build(std::move(ivs));
    if (dropped > 0) {
        std::fprintf(stderr,
            "[reference-gate] %zu NAHR-bed rows dropped (chrom not in .llmi)\n",
            dropped);
    }
    std::fprintf(stderr,
        "[reference-gate] %zu NAHR intervals indexed across %zu contigs\n",
        nahr_iv_.Size(), chrom_to_ref_id_.size());
    return true;
}

std::size_t ReferenceGate::NumNahrIntervals() const noexcept {
    return nahr_iv_.Size();
}

std::uint32_t ReferenceGate::LookupRefId(std::string_view chrom) const noexcept {
    auto it = chrom_to_ref_id_.find(std::string(chrom));
    return it == chrom_to_ref_id_.end() ? kNoRefId : it->second;
}

GateResult ReferenceGate::Classify(std::string_view read_seq) const {
    GateResult r;
    if (!idx_ || read_seq.size() < 100) {
        r.verdict = GateResult::Unmapped;
        return r;
    }

    // Cap to a sane number of seed hits; assemblies have huge contigs.
    constexpr std::size_t kMaxSeedHits = 1u << 18;  // 262144
    auto hits = idx_->Query(read_seq, kMaxSeedHits);
    if (hits.empty()) {
        r.verdict = GateResult::Unmapped;
        return r;
    }

    // Bucket seeds by (ref_id, ref_pos / cluster_window). Find the
    // single bucket with the most hits — that's the best coarse-locate.
    // We use an open hash map keyed on (ref_id, bucket_idx) packed into
    // uint64_t to keep this branchless and small.
    std::unordered_map<std::uint64_t, std::uint32_t> bucket_counts;
    bucket_counts.reserve(hits.size() / 4 + 16);

    const std::uint64_t window = cfg_.cluster_window;
    for (const auto& h : hits) {
        std::uint64_t b = (static_cast<std::uint64_t>(h.ref_id) << 32) |
                          (h.ref_pos / window);
        ++bucket_counts[b];
    }

    std::uint64_t best_key = 0;
    std::uint32_t best_count = 0;
    for (const auto& kv : bucket_counts) {
        if (kv.second > best_count) {
            best_count = kv.second;
            best_key   = kv.first;
        }
    }

    if (best_count < cfg_.min_seeds) {
        r.verdict = GateResult::Unmapped;
        r.n_seeds_in_bucket = static_cast<std::uint16_t>(best_count);
        return r;
    }

    const std::uint32_t best_ref_id = static_cast<std::uint32_t>(best_key >> 32);
    const std::uint32_t bucket_idx  = static_cast<std::uint32_t>(best_key);
    const std::uint32_t bucket_start = bucket_idx * cfg_.cluster_window;
    const std::uint32_t bucket_end   = bucket_start + cfg_.cluster_window;

    r.ref_id   = best_ref_id;
    r.ref_start = bucket_start;
    r.ref_end   = bucket_end;
    r.n_seeds_in_bucket =
        best_count > std::numeric_limits<std::uint16_t>::max()
        ? std::numeric_limits<std::uint16_t>::max()
        : static_cast<std::uint16_t>(best_count);

    // Probe NAHR annotations at the bucket boundaries.
    auto hits_start = nahr_iv_.IntervalsAt(best_ref_id, bucket_start);
    auto hits_end   = nahr_iv_.IntervalsAt(best_ref_id, bucket_end > 0 ? bucket_end - 1 : 0);

    const annot::AnnotationInterval* hit = nullptr;
    if (!hits_start.empty()) hit = hits_start.front();
    else if (!hits_end.empty()) hit = hits_end.front();

    if (hit) {
        r.verdict = GateResult::RouteToNahr;
        const auto pipe = hit->region_name.find('|');
        if (pipe != std::string::npos) {
            r.pair_id = hit->region_name.substr(0, pipe);
            r.arm     = hit->region_name.substr(pipe + 1);
        } else {
            r.pair_id = hit->region_name;
        }
    } else {
        r.verdict = GateResult::Skip;
    }
    return r;
}

}  // namespace llmap::junction_hunter
