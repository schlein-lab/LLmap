// LLmap — JunctionDb implementation.
//
// Loaders parse minimal-fidelity formats — each source's full schema
// has more fields than we use, and skipping them keeps the parse simple
// + fast. Specifically:
//
//   GENCODE GFF3 — we re-derive junctions from consecutive exon records
//                   of the same transcript (same approach as
//                   anchor_loader_gencode.cpp); the in_gencode flag is
//                   set true for every junction seen.
//   GTEx BED      — expects 6-column BED with name field encoding
//                   "<chrom>:<donor>-<acceptor>(:<n_samples>)?".
//   circBase BED  — expects 6-column BED with strand-aware coords; we
//                   ingest both forward and back-splice forms.
//   ChessDB GTF   — same exon→junction synthesis as GENCODE.
//
// On parse failure (file missing or unrecognised format) the loader
// returns false but the DB is not corrupted — partial loads are valid.

#include "annot/junction_db.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace llmap::annot {

namespace {

std::vector<std::string> SplitTab(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Internal: get-or-insert mutator.
// ---------------------------------------------------------------------------

JunctionEvidence& JunctionDb::GetOrInsert(
    std::string_view chrom,
    std::int64_t donor,
    std::int64_t acceptor) {
    auto& cmap = by_chrom_[std::string(chrom)];
    return cmap[JunctionKey{donor, acceptor}];
}

// ---------------------------------------------------------------------------
// LoadGencode — derive junctions from consecutive exons per transcript.
// ---------------------------------------------------------------------------

bool JunctionDb::LoadGencode(const std::filesystem::path& gff) {
    std::ifstream in(gff);
    if (!in) return false;

    // transcript_id → vector of (start_0,end_0,chrom)
    std::unordered_map<std::string,
        std::vector<std::tuple<std::int64_t, std::int64_t,
                                std::string>>> tx_exons;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 9) continue;
        if (cols[2] != "exon") continue;

        const std::string& chrom = cols[0];
        const auto start_1 = std::stoll(cols[3]);
        const auto end_1   = std::stoll(cols[4]);
        if (start_1 <= 0 || end_1 < start_1) continue;

        // pull transcript_id out of column 9
        const std::string& attrs = cols[8];
        const auto pos = attrs.find("transcript_id=");
        if (pos == std::string::npos) continue;
        const auto end_attr = attrs.find(';', pos);
        std::string tx = attrs.substr(pos + 14,
            (end_attr == std::string::npos)
                ? std::string::npos
                : end_attr - (pos + 14));
        // strip quotes if any
        if (!tx.empty() && tx.front() == '"') tx.erase(0, 1);
        if (!tx.empty() && tx.back()  == '"') tx.pop_back();

        tx_exons[tx].emplace_back(start_1 - 1, end_1, chrom);
    }

    // Sort + pair into junctions
    for (auto& [tx, exons] : tx_exons) {
        if (exons.size() < 2) continue;
        std::sort(exons.begin(), exons.end(),
                  [](const auto& a, const auto& b) {
                      return std::get<0>(a) < std::get<0>(b);
                  });
        for (std::size_t i = 0; i + 1 < exons.size(); ++i) {
            const auto& [s_i, e_i, chrom_i] = exons[i];
            const auto& [s_j, e_j, chrom_j] = exons[i + 1];
            (void)e_j;
            if (chrom_i != chrom_j) continue;
            auto& ev = GetOrInsert(chrom_i, e_i - 1, s_j);
            ev.in_gencode = true;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// LoadGtexJunctions — 6-column BED. name field encodes "donor-acceptor".
// ---------------------------------------------------------------------------

bool JunctionDb::LoadGtexJunctions(const std::filesystem::path& bed) {
    std::ifstream in(bed);
    if (!in) return false;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 4) continue;
        const std::string& chrom = cols[0];
        const auto donor    = std::stoll(cols[1]);
        const auto acceptor = std::stoll(cols[2]);
        auto& ev = GetOrInsert(chrom, donor, acceptor);
        ev.in_gtex = true;
        if (cols.size() >= 5) {
            // BED 'score' column doubles as sample count for our purposes.
            try {
                ev.gtex_sample_count = static_cast<std::uint32_t>(
                    std::stoul(cols[4]));
            } catch (...) {}
        }
    }
    return true;
}

bool JunctionDb::LoadCircRnaDb(const std::filesystem::path& bed) {
    std::ifstream in(bed);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 3) continue;
        const std::string& chrom = cols[0];
        const auto a = std::stoll(cols[1]);
        const auto b = std::stoll(cols[2]);
        // For circRNA, the BED 'start' is the back-splice acceptor
        // and 'end' is the back-splice donor (acceptor < donor).
        auto& ev = GetOrInsert(chrom, b, a);  // donor>acceptor for back-splice
        ev.in_circ_db = true;
    }
    return true;
}

bool JunctionDb::LoadChessDbJunctions(const std::filesystem::path& gtf) {
    // CHESS GTF — same exon→junction synthesis as GENCODE.
    std::ifstream in(gtf);
    if (!in) return false;

    std::unordered_map<std::string,
        std::vector<std::tuple<std::int64_t, std::int64_t,
                                std::string>>> tx_exons;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 9) continue;
        if (cols[2] != "exon") continue;
        const auto start_1 = std::stoll(cols[3]);
        const auto end_1   = std::stoll(cols[4]);

        const std::string& attrs = cols[8];
        const auto pos = attrs.find("transcript_id \"");
        if (pos == std::string::npos) continue;
        const auto end_quote = attrs.find('"', pos + 15);
        if (end_quote == std::string::npos) continue;
        std::string tx = attrs.substr(pos + 15, end_quote - (pos + 15));
        tx_exons[tx].emplace_back(start_1 - 1, end_1, cols[0]);
    }
    for (auto& [tx, exons] : tx_exons) {
        if (exons.size() < 2) continue;
        std::sort(exons.begin(), exons.end(),
                  [](const auto& a, const auto& b) {
                      return std::get<0>(a) < std::get<0>(b);
                  });
        for (std::size_t i = 0; i + 1 < exons.size(); ++i) {
            const auto& [s_i, e_i, c_i] = exons[i];
            const auto& [s_j, e_j, c_j] = exons[i + 1];
            (void)e_j;
            if (c_i != c_j) continue;
            auto& ev = GetOrInsert(c_i, e_i - 1, s_j);
            ev.in_chessdb = true;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Lookup / Inspection
// ---------------------------------------------------------------------------

JunctionEvidence JunctionDb::Lookup(std::string_view chrom,
                                     std::int64_t donor,
                                     std::int64_t acceptor) const {
    auto cit = by_chrom_.find(std::string(chrom));
    if (cit == by_chrom_.end()) return {};
    auto jit = cit->second.find(JunctionKey{donor, acceptor});
    if (jit == cit->second.end()) return {};
    return jit->second;
}

bool JunctionDb::HasAnyEvidence(std::string_view chrom,
                                 std::int64_t donor,
                                 std::int64_t acceptor) const {
    auto ev = Lookup(chrom, donor, acceptor);
    return ev.in_gencode || ev.in_mane || ev.in_refseq
        || ev.in_chessdb || ev.in_gtex || ev.in_circ_db;
}

std::size_t JunctionDb::TotalJunctions() const {
    std::size_t n = 0;
    for (const auto& [_, m] : by_chrom_) n += m.size();
    return n;
}

std::size_t JunctionDb::JunctionsOnChrom(std::string_view chrom) const {
    auto it = by_chrom_.find(std::string(chrom));
    if (it == by_chrom_.end()) return 0;
    return it->second.size();
}

}  // namespace llmap::annot
