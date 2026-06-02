// LLmap — TranscriptKmerIndex core (storage, lookup, hashing).
//
// Build logic lives in transcript_kmer_index_build.cpp so this file stays
// under the modular LOC cap. Both translation units share the
// implementation-private members of TranscriptKmerIndex via the header.

#include "index/transcript_kmer_index.h"

#include <algorithm>
#include <cctype>

namespace llmap::index {

namespace {

constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime  = 1099511628211ULL;

inline std::uint64_t FnvHash(std::string_view s) noexcept {
    // FNV-1a — fast, decent distribution. Same primitive the existing
    // MinimizerIndex uses (kept here as a private copy so this module
    // can build standalone).
    std::uint64_t h = kFnvOffset;
    for (unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= kFnvPrime;
    }
    return h;
}

inline bool IsClean(char c) noexcept {
    const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return u == 'A' || u == 'C' || u == 'G' || u == 'T';
}

}  // namespace

// ---------------------------------------------------------------------------
// KmerOrigin name table.
// ---------------------------------------------------------------------------

const char* KmerOriginName(KmerOrigin o) noexcept {
    switch (o) {
        case KmerOrigin::IntraExon:          return "intra_exon";
        case KmerOrigin::JunctionSpanning:   return "junction_spanning";
        case KmerOrigin::BackSpliceSpanning: return "back_splice_spanning";
        case KmerOrigin::SterileIntronic:    return "sterile_intronic";
        case KmerOrigin::PreMrnaIntronic:    return "pre_mrna_intronic";
        case KmerOrigin::ShortRna:           return "short_rna";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// Hashing primitive (static).
// ---------------------------------------------------------------------------

std::uint64_t TranscriptKmerIndex::HashKmer(std::string_view k) noexcept {
    // Upper-case canonicalisation in-place via a tiny stack buffer for
    // small k-mers; falls back to a heap copy for very long k-mers.
    constexpr std::size_t kStackBuf = 128;
    if (k.size() <= kStackBuf) {
        char buf[kStackBuf];
        for (std::size_t i = 0; i < k.size(); ++i) {
            buf[i] = static_cast<char>(
                std::toupper(static_cast<unsigned char>(k[i])));
        }
        return FnvHash(std::string_view(buf, k.size()));
    } else {
        std::string copy(k);
        for (char& c : copy) {
            c = static_cast<char>(
                std::toupper(static_cast<unsigned char>(c)));
        }
        return FnvHash(copy);
    }
}

// ---------------------------------------------------------------------------
// AddHit — honour max_occ cap.
// ---------------------------------------------------------------------------

void TranscriptKmerIndex::AddHit(KmerOrigin origin,
                                  std::uint64_t hash,
                                  TranscriptKmerHit hit) {
    hit.origin = origin;
    auto& table = [this, origin]() -> std::unordered_map<
        std::uint64_t, std::vector<TranscriptKmerHit>>& {
        switch (origin) {
            case KmerOrigin::IntraExon:          return intra_exon_;
            case KmerOrigin::JunctionSpanning:   return junction_spanning_;
            case KmerOrigin::BackSpliceSpanning: return backsplice_;
            case KmerOrigin::SterileIntronic:    return sterile_intronic_;
            case KmerOrigin::PreMrnaIntronic:    return premrna_intronic_;
            case KmerOrigin::ShortRna:           return short_rna_;
        }
        // unreachable; default to intra_exon
        return intra_exon_;
    }();

    auto& bucket = table[hash];
    if (cfg_.max_occ == 0 || bucket.size() < cfg_.max_occ) {
        bucket.push_back(hit);
    }
    // else: silently drop — repeat-rich locus capped.
}

// ---------------------------------------------------------------------------
// Lookup APIs.
// ---------------------------------------------------------------------------

std::vector<TranscriptKmerHit>
TranscriptKmerIndex::QueryGenomeWindow(std::string_view window_seq) const {
    std::vector<TranscriptKmerHit> out;
    if (window_seq.empty()) return out;

    // Genome-direction queries are restricted to k-mers whose origin
    // is meaningful on contiguous genomic DNA: IntraExon (must match
    // by definition) + PreMrnaIntronic (intronic genomic positions).
    auto query_table = [&](std::uint8_t k,
                            const std::unordered_map<std::uint64_t,
                                std::vector<TranscriptKmerHit>>& table) {
        if (k == 0 || window_seq.size() < k) return;
        for (std::size_t i = 0; i + k <= window_seq.size(); ++i) {
            // skip windows containing non-ACGT
            bool clean = true;
            for (std::size_t j = 0; j < k; ++j) {
                if (!IsClean(window_seq[i + j])) { clean = false; break; }
            }
            if (!clean) continue;
            const auto h = HashKmer(window_seq.substr(i, k));
            auto it = table.find(h);
            if (it == table.end()) continue;
            for (const auto& hit : it->second) out.push_back(hit);
        }
    };

    query_table(cfg_.k_intra, intra_exon_);
    if (cfg_.include_premrna_intronic) {
        query_table(cfg_.k_intra, premrna_intronic_);
    }
    return out;
}

std::vector<TranscriptKmerHit>
TranscriptKmerIndex::QueryReadKmer(std::uint64_t hash) const {
    std::vector<TranscriptKmerHit> out;
    auto append = [&](const auto& table) {
        auto it = table.find(hash);
        if (it == table.end()) return;
        out.insert(out.end(), it->second.begin(), it->second.end());
    };
    append(intra_exon_);
    append(junction_spanning_);
    append(backsplice_);
    append(sterile_intronic_);
    append(premrna_intronic_);
    append(short_rna_);
    return out;
}

std::vector<TranscriptKmerHit>
TranscriptKmerIndex::QueryReadKmer(std::string_view kmer_seq) const {
    return QueryReadKmer(HashKmer(kmer_seq));
}

// ---------------------------------------------------------------------------
// Telemetry + maintenance.
// ---------------------------------------------------------------------------

std::size_t TranscriptKmerIndex::TableSize(KmerOrigin o) const noexcept {
    switch (o) {
        case KmerOrigin::IntraExon:          return intra_exon_.size();
        case KmerOrigin::JunctionSpanning:   return junction_spanning_.size();
        case KmerOrigin::BackSpliceSpanning: return backsplice_.size();
        case KmerOrigin::SterileIntronic:    return sterile_intronic_.size();
        case KmerOrigin::PreMrnaIntronic:    return premrna_intronic_.size();
        case KmerOrigin::ShortRna:           return short_rna_.size();
    }
    return 0;
}

std::size_t TranscriptKmerIndex::TotalKmers() const noexcept {
    return intra_exon_.size() + junction_spanning_.size()
        + backsplice_.size() + sterile_intronic_.size()
        + premrna_intronic_.size() + short_rna_.size();
}

void TranscriptKmerIndex::Clear() {
    intra_exon_.clear();
    junction_spanning_.clear();
    backsplice_.clear();
    sterile_intronic_.clear();
    premrna_intronic_.clear();
    short_rna_.clear();
    cfg_ = {};
}

// Save/Load are stubbed for now — they will be implemented in a
// follow-up once the index has a stable on-disk schema we want to pin.
// For the tests below we exercise build/query only.
bool TranscriptKmerIndex::Save(const std::filesystem::path&) const {
    return false;
}
bool TranscriptKmerIndex::Load(const std::filesystem::path&) {
    return false;
}

}  // namespace llmap::index
