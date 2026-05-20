// LLmap — IGH anchor catalog implementation.

#include "igh_locus/igh_anchor_catalog.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace llmap::igh_locus {

namespace {

const std::array<std::string_view, 9> kIsotypes = {
    "IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"};

char Complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

std::string ReverseComplement(std::string_view s) {
    std::string out(s.size(), 'N');
    for (std::size_t i = 0; i < s.size(); ++i) {
        out[s.size() - 1 - i] = Complement(s[i]);
    }
    return out;
}

std::string ToUpper(std::string_view s) {
    std::string out(s);
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    }
    return out;
}

// Split a header token like "IGHG4_G056511_hap1_CH1_295bp" into (gene, copy,
// exon). copy_label = everything up to (not including) the exon token. The exon
// token is the first token starting with "CH" (or "M"/"H" for membrane/hinge).
void ParseAnchorName(const std::string& name, std::string& gene,
                     std::string& copy_label, std::string& exon) {
    std::vector<std::string> tok;
    std::stringstream ss(name);
    std::string t;
    while (std::getline(ss, t, '_')) tok.push_back(t);
    gene = tok.empty() ? "" : tok.front();
    // find exon token index
    std::size_t ex_idx = tok.size();
    for (std::size_t i = 1; i < tok.size(); ++i) {
        if (tok[i].rfind("CH", 0) == 0 || tok[i] == "M1" || tok[i] == "M2" ||
            tok[i] == "H") {
            ex_idx = i;
            break;
        }
    }
    exon = (ex_idx < tok.size()) ? tok[ex_idx] : "CH?";
    std::string cl;
    std::size_t upto = (ex_idx < tok.size()) ? ex_idx : tok.size();
    for (std::size_t i = 0; i < upto; ++i) {
        if (i) cl += "_";
        cl += tok[i];
    }
    copy_label = cl.empty() ? gene : cl;
}

// Parse optional "loc=contig:start-end" attribute from a header line tail.
bool ParseLoc(const std::string& header_tail, std::string& contig,
              std::uint64_t& start, std::uint64_t& end) {
    auto pos = header_tail.find("loc=");
    if (pos == std::string::npos) return false;
    std::string v = header_tail.substr(pos + 4);
    auto sp = v.find_first_of(" \t");
    if (sp != std::string::npos) v = v.substr(0, sp);
    auto colon = v.find(':');
    auto dash = v.find('-', colon == std::string::npos ? 0 : colon);
    if (colon == std::string::npos || dash == std::string::npos) return false;
    contig = v.substr(0, colon);
    try {
        start = std::stoull(v.substr(colon + 1, dash - colon - 1));
        end = std::stoull(v.substr(dash + 1));
    } catch (...) {
        return false;
    }
    return true;
}

}  // namespace

bool IsIghGene(std::string_view token) {
    for (auto g : kIsotypes) {
        if (token == g) return true;
    }
    return false;
}

std::optional<IghAnchorCatalog> IghAnchorCatalog::LoadFasta(
    const std::string& path, int max_mismatch) {
    std::ifstream in(path);
    if (!in) return std::nullopt;

    IghAnchorCatalog cat;
    cat.max_mismatch_ = max_mismatch < 0 ? 0 : max_mismatch;

    std::string line, cur_name, cur_tail, cur_seq;
    auto flush = [&]() {
        if (cur_name.empty()) return;
        std::string gene, copy_label, exon;
        ParseAnchorName(cur_name, gene, copy_label, exon);
        if (!IsIghGene(gene) || cur_seq.empty()) {
            cur_name.clear(); cur_tail.clear(); cur_seq.clear();
            return;
        }
        IghAnchor a;
        a.name = cur_name;
        a.gene = gene;
        a.copy_label = copy_label;
        a.exon = exon;
        a.seq = ToUpper(cur_seq);
        a.has_locus = ParseLoc(cur_tail, a.contig, a.start, a.end);
        cat.copy_to_gene_[copy_label] = gene;
        cat.anchors_.push_back(std::move(a));
        cur_name.clear(); cur_tail.clear(); cur_seq.clear();
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            std::string body = line.substr(1);
            auto sp = body.find_first_of(" \t");
            cur_name = (sp == std::string::npos) ? body : body.substr(0, sp);
            cur_tail = (sp == std::string::npos) ? "" : body.substr(sp + 1);
        } else {
            cur_seq += line;
        }
    }
    flush();

    if (cat.anchors_.empty()) return std::nullopt;

    // Dedupe by sequence so paralog-shared exons are scanned once.
    std::unordered_map<std::string, std::size_t> seq_to_unique;
    for (std::size_t i = 0; i < cat.anchors_.size(); ++i) {
        const std::string& s = cat.anchors_[i].seq;
        auto it = seq_to_unique.find(s);
        if (it == seq_to_unique.end()) {
            seq_to_unique[s] = cat.unique_.size();
            UniqueSeq u;
            u.fwd = s;
            u.rev = ReverseComplement(s);
            u.anchor_idx.push_back(i);
            cat.unique_.push_back(std::move(u));
        } else {
            cat.unique_[it->second].anchor_idx.push_back(i);
        }
    }
    return cat;
}

bool IghAnchorCatalog::ContainsSeq(std::string_view read,
                                   const std::string& probe) const {
    if (probe.size() > read.size()) return false;
    if (max_mismatch_ <= 0) {
        return read.find(probe) != std::string_view::npos;
    }
    // Seeded Hamming search. By pigeonhole, a probe with <=k mismatches over L
    // bases, split into k+1 non-overlapping tiles, has at least one error-free
    // tile. Try every tile as an exact seed, then verify the full window.
    const std::size_t L = probe.size();
    const std::size_t n_tiles = static_cast<std::size_t>(max_mismatch_) + 1;
    const std::size_t tile = L / n_tiles;
    if (tile < 8) return read.find(probe) != std::string_view::npos;
    for (std::size_t t = 0; t < n_tiles; ++t) {
        const std::size_t toff = t * tile;
        std::string_view seed_sv(probe.data() + toff, tile);
        std::size_t from = 0;
        while (true) {
            std::size_t hit = read.find(seed_sv, from);
            if (hit == std::string_view::npos) break;
            if (hit >= toff) {
                const std::size_t wstart = hit - toff;  // where probe[0] aligns
                if (wstart + L <= read.size()) {
                    int mm = 0;
                    for (std::size_t j = 0; j < L && mm <= max_mismatch_; ++j) {
                        if (read[wstart + j] != probe[j]) ++mm;
                    }
                    if (mm <= max_mismatch_) return true;
                }
            }
            from = hit + 1;
        }
    }
    return false;
}

IghMatch IghAnchorCatalog::Match(std::string_view read_seq) const {
    IghMatch m;
    if (read_seq.empty()) return m;

    // copy_label -> (distinct exon set, total hits, reverse?, gene)
    struct CopyAgg {
        std::vector<std::string> exons;
        int hits{0};
        bool reverse{false};
        std::string gene;
    };
    std::unordered_map<std::string, CopyAgg> per_copy;

    for (const auto& u : unique_) {
        bool fwd = ContainsSeq(read_seq, u.fwd);
        bool rev = !fwd && ContainsSeq(read_seq, u.rev);
        if (!fwd && !rev) continue;
        for (std::size_t ai : u.anchor_idx) {
            const IghAnchor& a = anchors_[ai];
            auto& agg = per_copy[a.copy_label];
            agg.gene = a.gene;
            agg.hits += 1;
            agg.reverse = rev;
            if (std::find(agg.exons.begin(), agg.exons.end(), a.exon) ==
                agg.exons.end()) {
                agg.exons.push_back(a.exon);
            }
        }
    }
    if (per_copy.empty()) return m;

    // Winning copy: most distinct exons, tie-break on total hits.
    const std::string* best = nullptr;
    const CopyAgg* best_agg = nullptr;
    for (const auto& [label, agg] : per_copy) {
        if (!best_agg ||
            agg.exons.size() > best_agg->exons.size() ||
            (agg.exons.size() == best_agg->exons.size() &&
             agg.hits > best_agg->hits)) {
            best = &label;
            best_agg = &agg;
        }
    }

    // Gene-level ambiguity: do two different genes tie on max distinct exons?
    std::unordered_map<std::string, std::size_t> gene_best_exons;
    for (const auto& [label, agg] : per_copy) {
        auto& g = gene_best_exons[agg.gene];
        g = std::max(g, agg.exons.size());
    }
    std::size_t top = best_agg->exons.size();
    int n_genes_at_top = 0;
    for (const auto& [g, e] : gene_best_exons) {
        if (e == top) ++n_genes_at_top;
    }

    m.matched = true;
    m.gene = best_agg->gene;
    m.copy_label = *best;
    m.n_anchor_hits = best_agg->hits;
    m.n_distinct_exons = static_cast<int>(best_agg->exons.size());
    m.ambiguous_gene = n_genes_at_top > 1;
    m.is_reverse = best_agg->reverse;
    m.matched_exons = best_agg->exons;

    // Carry the winning copy's genomic locus if any anchor had one.
    for (const auto& a : anchors_) {
        if (a.copy_label == m.copy_label && a.has_locus) {
            m.contig = a.contig;
            m.start = a.start;
            m.end = a.end;
            break;
        }
    }
    return m;
}

std::optional<std::string> IghAnchorCatalog::GeneOfTarget(
    std::string_view target_id) const {
    // exact copy label
    auto it = copy_to_gene_.find(std::string(target_id));
    if (it != copy_to_gene_.end()) return it->second;
    // leading token is an isotype (covers transcript refs like "IGHG4_xxx")
    auto us = target_id.find('_');
    std::string_view lead =
        (us == std::string_view::npos) ? target_id : target_id.substr(0, us);
    if (IsIghGene(lead)) return std::string(lead);
    return std::nullopt;
}

}  // namespace llmap::igh_locus
