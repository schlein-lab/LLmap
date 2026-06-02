// LLmap — ExpressionDb implementation.
//
// Best-effort TSV / GCT readers. Each loader returns true if at least
// the first non-empty data line parsed correctly; partial loads are
// preserved (any entries already in tpm_/cell_type_frac_ stay).

#include "fusion/expression_db.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace llmap::fusion {

namespace {

std::vector<std::string> SplitTab(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}

float ParseFloatTolerant(const std::string& s) {
    if (s.empty() || s == "NA" || s == "NaN") return -1.0f;
    try {
        return std::stof(s);
    } catch (...) {
        return -1.0f;
    }
}

}  // namespace

// ===========================================================================
// LoadGtexBulkGct
// ===========================================================================
// GCT format:
//   line 1: "#1.2" (version)
//   line 2: "<n_genes>\t<n_samples>"
//   line 3: header: Name  Description  Tissue1  Tissue2  ...
//   line 4+: gene_id  gene_name  tpm1  tpm2  ...
//
// We treat the "Name" column as transcript_id (for GTEx v8 medianTpm
// matrices it's actually gene_id; the user supplies an alternate file
// for transcript-level data via LoadGtexLongReadTsv).
bool ExpressionDb::LoadGtexBulkGct(const std::filesystem::path& gct) {
    std::ifstream in(gct);
    if (!in) return false;

    std::string line;
    // skip the two GCT preamble lines
    if (!std::getline(in, line)) return false;
    if (!std::getline(in, line)) return false;

    // header line
    if (!std::getline(in, line)) return false;
    auto cols = SplitTab(line);
    if (cols.size() < 3) return false;  // need at least 1 tissue column
    std::vector<std::string> tissues(cols.begin() + 2, cols.end());

    while (std::getline(in, line)) {
        auto row = SplitTab(line);
        if (row.size() < 3) continue;
        const std::string& tx_id = row[0];
        for (std::size_t i = 0; i + 2 < row.size() && i < tissues.size(); ++i) {
            const float v = ParseFloatTolerant(row[i + 2]);
            if (v < 0.0f) continue;
            tpm_[Key{tx_id, tissues[i]}] = v;
        }
    }
    return true;
}

// ===========================================================================
// LoadGtexLongReadTsv
//   columns: transcript_id  tissue  tpm
// ===========================================================================
bool ExpressionDb::LoadGtexLongReadTsv(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; if (line.find('\t') == std::string::npos) continue; }
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 3) continue;
        const float v = ParseFloatTolerant(cols[2]);
        if (v < 0.0f) continue;
        tpm_[Key{cols[0], cols[1]}] = v;
    }
    return true;
}

// ===========================================================================
// LoadTabulaSapiensTsv / LoadHcaTsv
//   columns: transcript_id  organ  cell_type  mean_expression
//
// We treat 'mean_expression' as a fraction-of-cells-expressing proxy
// and also fold it into tpm_ by treating organ as a tissue label.
// ===========================================================================
bool ExpressionDb::LoadTabulaSapiensTsv(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; if (line.find('\t') == std::string::npos) continue; }
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 4) continue;
        const float v = ParseFloatTolerant(cols[3]);
        if (v < 0.0f) continue;
        cell_type_frac_[CtKey{cols[0], cols[1], cols[2]}] = v;
        // Also fold into tpm_ as a tissue-level summary (max across
        // cell types within the organ) so callers without cell-type
        // info still get a useful expression prior.
        auto& m = tpm_[Key{cols[0], cols[1]}];
        m = std::max(m, v);
    }
    return true;
}

bool ExpressionDb::LoadHcaTsv(const std::filesystem::path& tsv) {
    // Same shape as Tabula Sapiens — re-use that parser.
    return LoadTabulaSapiensTsv(tsv);
}

// ===========================================================================
// LoadHpaTsv
//   Variable schema; we look for columns named "Gene" / "Tissue" / "TPM".
//   For simplicity we ingest the first four columns as
//   gene_id / gene_name / tissue / tpm.
// ===========================================================================
bool ExpressionDb::LoadHpaTsv(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; continue; }
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 4) continue;
        const float v = ParseFloatTolerant(cols[3]);
        if (v < 0.0f) continue;
        tpm_[Key{cols[0], cols[2]}] = v;
    }
    return true;
}

// ===========================================================================
// LoadRecount3Tsv
//   columns: gene_id  tissue  median_tpm
// ===========================================================================
bool ExpressionDb::LoadRecount3Tsv(const std::filesystem::path& tsv) {
    std::ifstream in(tsv);
    if (!in) return false;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; if (line.find('\t') == std::string::npos) continue; }
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 3) continue;
        const float v = ParseFloatTolerant(cols[2]);
        if (v < 0.0f) continue;
        tpm_[Key{cols[0], cols[1]}] = v;
    }
    return true;
}

// ===========================================================================
// Lookups
// ===========================================================================

float ExpressionDb::ExpectedTpm(std::string_view transcript_id,
                                  std::string_view tissue) const {
    auto it = tpm_.find(Key{std::string(transcript_id), std::string(tissue)});
    if (it == tpm_.end()) return -1.0f;
    return it->second;
}

float ExpressionDb::ExpectedCellTypeFraction(
    std::string_view transcript_id,
    std::string_view tissue,
    std::string_view cell_type) const {
    auto it = cell_type_frac_.find(CtKey{
        std::string(transcript_id),
        std::string(tissue),
        std::string(cell_type)});
    if (it == cell_type_frac_.end()) return -1.0f;
    return it->second;
}

float ExpressionDb::ExpressionEntropy(std::string_view transcript_id) const {
    // Shannon entropy over the per-tissue TPM values for the transcript.
    // Higher entropy = broader expression.
    std::vector<float> vals;
    float total = 0.0f;
    const std::string tx(transcript_id);
    for (const auto& [k, v] : tpm_) {
        if (k.transcript_id == tx && v > 0.0f) {
            vals.push_back(v);
            total += v;
        }
    }
    if (vals.empty() || total <= 0.0f) return 0.0f;
    float h = 0.0f;
    for (float v : vals) {
        const float p = v / total;
        if (p > 0.0f) h -= p * std::log2(p);
    }
    return h;
}

std::size_t ExpressionDb::TotalEntries() const noexcept {
    return tpm_.size() + cell_type_frac_.size();
}

std::size_t ExpressionDb::DistinctTranscripts() const noexcept {
    std::unordered_set<std::string> seen;
    for (const auto& [k, _] : tpm_) seen.insert(k.transcript_id);
    return seen.size();
}

std::size_t ExpressionDb::DistinctTissues() const noexcept {
    std::unordered_set<std::string> seen;
    for (const auto& [k, _] : tpm_) seen.insert(k.tissue);
    return seen.size();
}

void ExpressionDb::Clear() {
    tpm_.clear();
    cell_type_frac_.clear();
}

}  // namespace llmap::fusion
