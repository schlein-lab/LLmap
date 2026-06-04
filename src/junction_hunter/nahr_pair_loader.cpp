// LLmap — junction_hunter: NAHR-pair TSV loader implementation.

#include "junction_hunter/nahr_pair_loader.h"

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace llmap::junction_hunter {

namespace {

bool ParseUint(const std::string& s, std::uint64_t& out) {
    try { out = std::stoull(s); return true; } catch (...) { return false; }
}

bool ParseFloat(const std::string& s, float& out) {
    try { out = std::stof(s); return true; } catch (...) { return false; }
}

std::vector<std::string> SplitTabs(const std::string& line) {
    std::vector<std::string> cols;
    std::string cur;
    for (char c : line) {
        if (c == '\t') { cols.push_back(std::move(cur)); cur.clear(); }
        else { cur.push_back(c); }
    }
    cols.push_back(std::move(cur));
    return cols;
}

}  // namespace

PairLoaderStatus LoadNahrPairsTsv(const std::filesystem::path& tsv,
                                  std::vector<NahrPair>& out) {
    PairLoaderStatus st;
    std::ifstream in(tsv);
    if (!in) { st.error = "cannot open: " + tsv.string(); return st; }

    std::string line;
    if (!std::getline(in, line)) { st.error = "empty file"; return st; }
    auto header = SplitTabs(line);
    std::unordered_map<std::string, std::size_t> col;
    for (std::size_t i = 0; i < header.size(); ++i) col[header[i]] = i;

    const char* required[] = {
        "pair_id", "chrom", "lcr_up_start", "lcr_up_end",
        "lcr_down_start", "lcr_down_end",
        "interior_start", "interior_end", "identity"
    };
    for (auto* k : required) {
        if (!col.contains(k)) { st.error = std::string("missing column: ") + k; return st; }
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto cols = SplitTabs(line);
        if (cols.size() < header.size()) { ++st.records_skipped; continue; }

        NahrPair p;
        p.pair_id = cols[col["pair_id"]];
        p.chrom   = cols[col["chrom"]];
        bool good = ParseUint(cols[col["lcr_up_start"]],   p.lcr_up_start)
                 && ParseUint(cols[col["lcr_up_end"]],     p.lcr_up_end)
                 && ParseUint(cols[col["lcr_down_start"]], p.lcr_down_start)
                 && ParseUint(cols[col["lcr_down_end"]],   p.lcr_down_end)
                 && ParseUint(cols[col["interior_start"]], p.interior_start)
                 && ParseUint(cols[col["interior_end"]],   p.interior_end)
                 && ParseFloat(cols[col["identity"]],      p.lcr_identity);
        if (!good) { ++st.records_skipped; continue; }
        if (col.contains("interior_len_bp")) {
            std::uint64_t bp = 0;
            if (ParseUint(cols[col["interior_len_bp"]], bp)) p.interior_kb = static_cast<std::uint32_t>(bp / 1000);
        } else if (col.contains("interior_kb")) {
            std::uint64_t kb = 0;
            if (ParseUint(cols[col["interior_kb"]], kb)) p.interior_kb = static_cast<std::uint32_t>(kb);
        }
        out.push_back(std::move(p));
        ++st.records_loaded;
    }

    st.ok = true;
    return st;
}

}  // namespace llmap::junction_hunter
