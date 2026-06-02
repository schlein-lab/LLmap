// LLmap — SegDupCatalog core implementation (containers, indices, lookup).
//
// Per-file parser code lives in segdup_catalog_parse.cpp and bulk loader
// in bulk_loader.cpp to keep individual translation units under the
// repo's 400-LOC soft cap.

#include "catalog/segdup_catalog.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <system_error>

namespace llmap::catalog {

// ---------------------------------------------------------------------------
// SegDupCatalogEntry geometric predicates
// ---------------------------------------------------------------------------

bool SegDupCatalogEntry::contains(std::string_view assembly,
                                  std::string_view chrom,
                                  std::int64_t pos) const {
    auto it = coords_by_assembly.find(std::string(assembly));
    if (it == coords_by_assembly.end()) return false;
    const auto& c = it->second;
    if (c.chrom != chrom) return false;
    return pos >= c.start && pos < c.end;
}

bool SegDupCatalogEntry::overlaps(std::string_view assembly,
                                  std::string_view chrom,
                                  std::int64_t start,
                                  std::int64_t end) const {
    if (end <= start) return false;
    auto it = coords_by_assembly.find(std::string(assembly));
    if (it == coords_by_assembly.end()) return false;
    const auto& c = it->second;
    if (c.chrom != chrom) return false;
    // [start, end) overlaps [c.start, c.end) iff start < c.end && c.start < end.
    return start < c.end && c.start < end;
}

// ---------------------------------------------------------------------------
// SegDupCatalog — indexing
// ---------------------------------------------------------------------------

std::string SegDupCatalog::chrom_key(std::string_view assembly,
                                     std::string_view chrom) {
    std::string key;
    key.reserve(assembly.size() + 1 + chrom.size());
    key.append(assembly);
    key.push_back('\0');
    key.append(chrom);
    return key;
}

void SegDupCatalog::reindex_entry(std::size_t idx) {
    const auto& e = entries_[idx];
    by_locus_id_[e.locus_id].push_back(idx);
    for (const auto& [asm_name, coords] : e.coords_by_assembly) {
        if (coords.chrom.empty()) continue;
        const std::string key = chrom_key(asm_name, coords.chrom);
        auto& bucket = by_chrom_[key];
        // Insert maintaining sort by start.
        auto pos = std::lower_bound(
            bucket.begin(), bucket.end(), coords.start,
            [this, &asm_name](std::size_t a, std::int64_t s) {
                const auto& ec = entries_[a].coords_by_assembly.at(asm_name);
                return ec.start < s;
            });
        bucket.insert(pos, idx);
    }
}

void SegDupCatalog::add_entry(SegDupCatalogEntry entry) {
    entries_.push_back(std::move(entry));
    reindex_entry(entries_.size() - 1);
}

// ---------------------------------------------------------------------------
// SegDupCatalog — file loaders
// ---------------------------------------------------------------------------

namespace {

std::optional<std::string> read_file(const std::filesystem::path& path,
                                     std::string& err) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        err = "cannot open " + path.string();
        return std::nullopt;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

LoadStatus SegDupCatalog::load_curated_file(const std::filesystem::path& path) {
    LoadStatus st;
    std::string err;
    auto text = read_file(path, err);
    if (!text) {
        st.error = err;
        return st;
    }
    auto entry = parse_curated_json(*text, path, err);
    if (!entry) {
        st.error = err;
        return st;
    }
    add_entry(std::move(*entry));
    st.ok = true;
    st.records_loaded = 1;
    return st;
}

LoadStatus SegDupCatalog::load_curated_dir(const std::filesystem::path& dir) {
    LoadStatus st;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) {
        st.error = "not a directory: " + dir.string();
        return st;
    }
    std::vector<std::filesystem::path> json_files;
    for (const auto& de : std::filesystem::directory_iterator(dir, ec)) {
        if (!de.is_regular_file()) continue;
        if (de.path().extension() == ".json") json_files.push_back(de.path());
    }
    // Stable order for deterministic indices in tests.
    std::sort(json_files.begin(), json_files.end());
    for (const auto& p : json_files) {
        auto sub = load_curated_file(p);
        if (sub.ok) {
            st.records_loaded += sub.records_loaded;
        } else {
            st.records_skipped += 1;
            if (!st.error.empty()) st.error += "; ";
            st.error += p.filename().string() + ": " + sub.error;
        }
    }
    st.ok = (st.records_loaded > 0) || json_files.empty();
    return st;
}

// ---------------------------------------------------------------------------
// SegDupCatalog — lookups
// ---------------------------------------------------------------------------

std::optional<SegDupCatalogEntry>
SegDupCatalog::lookup_by_coords(std::string_view assembly,
                                std::string_view chrom,
                                std::int64_t pos) const {
    auto it = by_chrom_.find(chrom_key(assembly, chrom));
    if (it == by_chrom_.end()) return std::nullopt;
    // Multiple curated entries can legitimately overlap the same
    // coordinate — e.g. a small IGHG4 tandem-dup nested inside the
    // larger IGHG canondup NAHR block. Resolution rule: return the
    // most specific entry (smallest [start, end) span). Tie-broken by
    // earliest start, then by lexicographic locus_id for determinism.
    const SegDupCatalogEntry* best = nullptr;
    std::int64_t best_span = std::numeric_limits<std::int64_t>::max();
    for (std::size_t idx : it->second) {
        const auto& e = entries_[idx];
        if (!e.contains(assembly, chrom, pos)) continue;
        const auto& c = e.coords_by_assembly.at(std::string(assembly));
        const std::int64_t span = c.end - c.start;
        if (span < best_span ||
            (span == best_span && best &&
             (c.start < best->coords_by_assembly.at(std::string(assembly)).start ||
              (c.start == best->coords_by_assembly.at(std::string(assembly)).start &&
               e.locus_id < best->locus_id)))) {
            best = &e;
            best_span = span;
        }
    }
    if (!best) return std::nullopt;
    return *best;
}

std::vector<SegDupCatalogEntry>
SegDupCatalog::lookup_overlapping(std::string_view assembly,
                                  std::string_view chrom,
                                  std::int64_t start,
                                  std::int64_t end) const {
    std::vector<SegDupCatalogEntry> out;
    auto it = by_chrom_.find(chrom_key(assembly, chrom));
    if (it == by_chrom_.end()) return out;
    // Bucket is sorted by start. Once an entry's start >= query end no
    // *later* entry can begin earlier — but we cannot early-break on
    // containment because a long preceding entry could still contain a
    // later position. For overlap with [start, end) it is safe to break
    // when an entry starts at or past `end` because all further entries
    // start there or later.
    for (std::size_t idx : it->second) {
        const auto& e = entries_[idx];
        const auto& c = e.coords_by_assembly.at(std::string(assembly));
        if (c.start >= end) break;
        if (e.overlaps(assembly, chrom, start, end)) out.push_back(e);
    }
    return out;
}

std::vector<SegDupCatalogEntry>
SegDupCatalog::lookup_by_locus_id(std::string_view id) const {
    std::vector<SegDupCatalogEntry> out;
    auto it = by_locus_id_.find(std::string(id));
    if (it == by_locus_id_.end()) return out;
    out.reserve(it->second.size());
    for (std::size_t idx : it->second) out.push_back(entries_[idx]);
    return out;
}

std::size_t
SegDupCatalog::count_architecture(std::string_view architecture) const {
    std::size_t n = 0;
    for (const auto& e : entries_) {
        if (e.structural_architecture == architecture) ++n;
    }
    return n;
}

}  // namespace llmap::catalog
