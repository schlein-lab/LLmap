// LLmap — AnchorStore core implementation (storage + lookup, no I/O).
//
// I/O lives in companion translation units to keep this file under the
// 400-LOC modular cap and let each source loader evolve independently:
//   anchor_loader_gencode.cpp     — GENCODE GFF3
//   anchor_loader_mane.cpp        — MANE Select
//   anchor_loader_imgt.cpp        — IMGT/GENE-DB
//   anchor_loader_pangenome.cpp   — HPRC per-haplotype
//   anchor_loader_branch.cpp      — BRANCH bubble bridge
//
// This file owns the AnchorStore class itself (add, lookup, reindex,
// clear). Pure data movement; no parsing logic.

#include "anchor/anchor_store.h"

#include <algorithm>

namespace llmap::anchor {

// ---------------------------------------------------------------------------
// Mutation
// ---------------------------------------------------------------------------

std::uint32_t AnchorStore::AddAnchor(AnchorRecord rec) {
    const auto idx = static_cast<std::uint32_t>(anchors_.size());
    anchors_.push_back(std::move(rec));

    // Incremental index update — avoids the O(N) full Reindex() cost.
    const AnchorRecord& r = anchors_[idx];
    by_id_.emplace(r.anchor_id, idx);
    if (r.has_genomic_coords()) {
        by_chrom_[*r.ref_chrom].push_back(idx);
    }
    for (const auto& tag : r.tags) {
        by_tag_[tag].push_back(idx);
    }
    if (!r.transcript_id.empty()) {
        by_transcript_[r.transcript_id].push_back(idx);
    }
    return idx;
}

// ---------------------------------------------------------------------------
// Lookup
// ---------------------------------------------------------------------------

const AnchorRecord* AnchorStore::ById(std::string_view id) const {
    auto it = by_id_.find(std::string(id));
    if (it == by_id_.end()) return nullptr;
    return &anchors_[it->second];
}

std::vector<std::uint32_t>
AnchorStore::ByTag(std::string_view tag) const {
    auto it = by_tag_.find(std::string(tag));
    if (it == by_tag_.end()) return {};
    return it->second;
}

std::vector<std::uint32_t>
AnchorStore::ByRegion(std::string_view chrom,
                       std::int64_t start,
                       std::int64_t end) const {
    auto it = by_chrom_.find(std::string(chrom));
    if (it == by_chrom_.end()) return {};

    std::vector<std::uint32_t> hits;
    // by_chrom_ values are NOT sorted yet (Reindex sorts them).
    // We assume Reindex has been called; the linear scan is on a small
    // per-chrom subset (chr14 has ~30k transcripts in GENCODE, so ~200k
    // exon-level anchors — still trivial relative to alignment work).
    for (std::uint32_t idx : it->second) {
        const auto& a = anchors_[idx];
        if (!a.has_genomic_coords()) continue;
        // half-open overlap test
        if (*a.ref_end <= start) continue;
        if (*a.ref_start >= end) continue;
        hits.push_back(idx);
    }
    return hits;
}

std::vector<std::uint32_t>
AnchorStore::ByTranscriptId(std::string_view transcript_id) const {
    auto it = by_transcript_.find(std::string(transcript_id));
    if (it == by_transcript_.end()) return {};
    return it->second;
}

void AnchorStore::ForEach(
    const std::function<bool(const AnchorRecord&)>& pred,
    const std::function<void(std::uint32_t,
                              const AnchorRecord&)>& cb) const {
    for (std::uint32_t i = 0; i < anchors_.size(); ++i) {
        const auto& a = anchors_[i];
        if (pred(a)) cb(i, a);
    }
}

std::size_t AnchorStore::CountBySource(AnchorSource s) const {
    std::size_t n = 0;
    for (const auto& a : anchors_) {
        if (a.source == s) ++n;
    }
    return n;
}

// ---------------------------------------------------------------------------
// Maintenance
// ---------------------------------------------------------------------------

void AnchorStore::Reindex() {
    by_id_.clear();
    by_chrom_.clear();
    by_tag_.clear();
    by_transcript_.clear();

    for (std::uint32_t i = 0; i < anchors_.size(); ++i) {
        const auto& r = anchors_[i];
        by_id_.emplace(r.anchor_id, i);
        if (r.has_genomic_coords()) {
            by_chrom_[*r.ref_chrom].push_back(i);
        }
        for (const auto& tag : r.tags) {
            by_tag_[tag].push_back(i);
        }
        if (!r.transcript_id.empty()) {
            by_transcript_[r.transcript_id].push_back(i);
        }
    }

    // Sort per-chrom indices by ref_start so ByRegion can early-exit
    // once it walks past the query range. Not yet exploited by ByRegion
    // (current impl is O(n_per_chrom)), but sorting now makes a future
    // binary-search variant trivial.
    for (auto& [chrom, idx_vec] : by_chrom_) {
        (void)chrom;
        std::sort(idx_vec.begin(), idx_vec.end(),
                  [this](std::uint32_t a, std::uint32_t b) {
                      return *anchors_[a].ref_start
                          <  *anchors_[b].ref_start;
                  });
    }
}

void AnchorStore::Clear() {
    anchors_.clear();
    by_id_.clear();
    by_chrom_.clear();
    by_tag_.clear();
    by_transcript_.clear();
}

}  // namespace llmap::anchor
