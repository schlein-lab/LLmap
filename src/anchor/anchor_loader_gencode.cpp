// LLmap — GENCODE GFF3 anchor loader.
//
// Streams a (possibly gzipped) GFF3 line-by-line and emits one
// AnchorRecord per (transcript_id, exon) tuple. Filters aggressively at
// the line-tokenisation stage so memory stays O(transcripts × exons)
// not O(GFF3 file size); GENCODE v46 GFF3 is ~1.5 GB uncompressed but
// only ~3 % of lines are exon features we care about here.
//
// We do NOT pull sequences from the reference FASTA in this pass — that
// happens on-demand from the k-mer index. Adding sequences here would
// blow memory for no benefit; the splice/junction-spanning k-mer
// indexer (Plan-Block 3) walks the FASTA once and fills sequences then.
//
// Implementation notes:
//
//   - Gzipped input via zlib's gzopen — same dependency we already pull
//     in for FASTA reading. No new dep.
//
//   - Attribute parsing is hand-rolled: GFF3 uses
//       key=value[;key=value]
//     with URL-encoded values. We decode the bare minimum: %3B → ';' and
//     %3D → '='. Anything fancier (CDATA, Unicode escapes) doesn't show
//     up in GENCODE so we deliberately don't slow down for it.
//
//   - Per-transcript intron derivation: GENCODE annotates exons but not
//     introns; we synthesise ExonBoundary records by sorting the exons
//     of each transcript by start and pairing consecutive ones. Strand
//     handling: for '-' strand transcripts we still record donor/acceptor
//     in genomic order (donor_genomic_pos < acceptor_genomic_pos for
//     canonical splice); downstream consumers know to flip if they
//     need transcript-orientation.

#include "anchor/anchor_store.h"

#include "core/transcript_kind.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <zlib.h>

namespace llmap::anchor {

namespace {

// ---------------------------------------------------------------------------
// Minimal URL-decode for GFF3 attribute values.
// ---------------------------------------------------------------------------
std::string UrlDecodeGff3(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '%' && i + 2 < s.size()) {
            const char hi = s[i + 1];
            const char lo = s[i + 2];
            auto hex = [](char c) -> int {
                if (c >= '0' && c <= '9') return c - '0';
                if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
                if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
                return -1;
            };
            const int a = hex(hi);
            const int b = hex(lo);
            if (a >= 0 && b >= 0) {
                out.push_back(static_cast<char>((a << 4) | b));
                i += 2;
                continue;
            }
        }
        out.push_back(s[i]);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Parse a GFF3 attribute column into key→value map.
// ---------------------------------------------------------------------------
std::unordered_map<std::string, std::string>
ParseAttributes(std::string_view attr_col) {
    std::unordered_map<std::string, std::string> kv;
    std::size_t i = 0;
    while (i < attr_col.size()) {
        // find ';' or end
        std::size_t j = attr_col.find(';', i);
        if (j == std::string_view::npos) j = attr_col.size();

        // within [i,j): key=value
        auto eq = attr_col.find('=', i);
        if (eq != std::string_view::npos && eq < j) {
            std::string k(attr_col.substr(i, eq - i));
            std::string v = UrlDecodeGff3(attr_col.substr(eq + 1, j - eq - 1));
            kv.emplace(std::move(k), std::move(v));
        }
        i = j + 1;
    }
    return kv;
}

// ---------------------------------------------------------------------------
// Map GENCODE biotype string → TranscriptKind.
//
// GENCODE uses many biotypes; we map them to our cleaner taxonomy.
// Unmatched biotypes fall through to TranscriptKind::Unknown — that's
// fine, the later TranscriptKindClassifier (Block 2.5) does sequence-
// based refinement.
// ---------------------------------------------------------------------------
core::TranscriptKind BiotypeToKind(std::string_view biotype) {
    using K = core::TranscriptKind;
    if (biotype == "protein_coding") return K::MatureMrna;
    if (biotype == "lncRNA" || biotype == "lincRNA"
        || biotype == "intronic_lncRNA"
        || biotype == "antisense_lncRNA") return K::Lncrna;
    if (biotype == "antisense") return K::Antisense;
    if (biotype == "miRNA") return K::Mirna;
    if (biotype == "piRNA") return K::Pirna;
    if (biotype == "siRNA") return K::Sirna;
    if (biotype == "snRNA")  return K::Snrna_Major;  // refined later
    if (biotype == "snoRNA") return K::Snorna_CDbox; // refined later
    if (biotype == "scaRNA") return K::Scarna;
    if (biotype == "rRNA")   return K::Rrna;
    if (biotype == "Mt_rRNA" || biotype == "Mt_tRNA") return K::Mitochondrial;
    if (biotype == "tRNA")   return K::Trna;
    if (biotype == "vaultRNA")   return K::Vaultrna;
    if (biotype == "Y_RNA")      return K::Yrna;
    if (biotype == "ribozyme")   return K::Rmrp_Rnasep;
    if (biotype == "IG_C_gene" || biotype == "IG_C_pseudogene"
        || biotype == "IG_J_gene" || biotype == "IG_J_pseudogene"
        || biotype == "IG_V_gene" || biotype == "IG_V_pseudogene"
        || biotype == "IG_D_gene"
        || biotype == "TR_C_gene" || biotype == "TR_J_gene"
        || biotype == "TR_V_gene" || biotype == "TR_D_gene"
        || biotype == "TR_V_pseudogene" || biotype == "TR_J_pseudogene") {
        return K::MatureMrna;  // VDJ kinds are mature mRNA at exon level;
                                // sterile-germline / class-switch happens
                                // at read level (in the classifier).
    }
    return K::Unknown;
}

}  // namespace

// ---------------------------------------------------------------------------
// LoadGencodeGff
// ---------------------------------------------------------------------------
//
// Two-pass over the GFF3:
//   pass 1 — build map<transcript_id, biotype>  from feature=transcript
//   pass 2 — emit AnchorRecord per feature=exon, pull biotype from map
//
// We could do it single-pass and back-fill biotypes, but two passes is
// simpler and the GFF3 is read-once cached in OS page cache between
// passes so the IO cost is one read.
// ---------------------------------------------------------------------------
LoadStatus AnchorStore::LoadGencodeGff(
    const std::filesystem::path& gff,
    const std::filesystem::path& /*ref_fa*/,
    bool /*with_sequence*/) {

    LoadStatus status;
    if (!std::filesystem::exists(gff)) {
        status.error = "GFF3 not found: " + gff.string();
        return status;
    }

    gzFile fh = gzopen(gff.c_str(), "rb");
    if (!fh) {
        status.error = "gzopen failed: " + gff.string();
        return status;
    }

    // pass 1: transcript_id → biotype (+ gene_id, transcript_name)
    struct TxMeta {
        std::string biotype;
        std::string gene_id;
        std::string gene_name;
        char strand{'.'};
    };
    std::unordered_map<std::string, TxMeta> tx_meta;

    constexpr std::size_t kBufSize = 1 << 16;
    char buf[kBufSize];

    auto split_columns = [](std::string_view line,
                            std::array<std::string_view, 9>& out) -> bool {
        std::size_t start = 0;
        for (std::size_t col = 0; col < 9; ++col) {
            std::size_t tab = (col == 8) ? line.size()
                                          : line.find('\t', start);
            if (tab == std::string_view::npos) return false;
            out[col] = line.substr(start, tab - start);
            start = tab + 1;
        }
        return true;
    };

    while (gzgets(fh, buf, kBufSize)) {
        if (buf[0] == '#') continue;
        std::string_view line(buf);
        // strip trailing newline
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.remove_suffix(1);
        }
        std::array<std::string_view, 9> col;
        if (!split_columns(line, col)) continue;
        if (col[2] != "transcript") continue;

        auto attrs = ParseAttributes(col[8]);
        auto tid_it = attrs.find("transcript_id");
        if (tid_it == attrs.end()) continue;

        TxMeta m;
        if (auto bt = attrs.find("transcript_type"); bt != attrs.end()) {
            m.biotype = bt->second;
        }
        if (auto gid = attrs.find("gene_id"); gid != attrs.end()) {
            m.gene_id = gid->second;
        }
        if (auto gn = attrs.find("gene_name"); gn != attrs.end()) {
            m.gene_name = gn->second;
        }
        if (!col[6].empty()) m.strand = col[6][0];
        tx_meta[std::move(tid_it->second)] = std::move(m);
    }
    gzrewind(fh);

    // pass 2: emit exon records grouped by transcript so we can build
    // ExonBoundary lists.
    struct TxExons {
        std::vector<std::pair<std::int64_t, std::int64_t>> exons;  // [start, end)
        std::string chrom;
        char strand{'.'};
    };
    std::unordered_map<std::string, TxExons> exons_by_tx;

    while (gzgets(fh, buf, kBufSize)) {
        if (buf[0] == '#') continue;
        std::string_view line(buf);
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.remove_suffix(1);
        }
        std::array<std::string_view, 9> col;
        if (!split_columns(line, col)) continue;
        if (col[2] != "exon") continue;

        auto attrs = ParseAttributes(col[8]);
        auto tid_it = attrs.find("transcript_id");
        if (tid_it == attrs.end()) continue;

        // GFF3 is 1-based inclusive; convert to 0-based half-open.
        std::int64_t start_1 = std::atoll(std::string(col[3]).c_str());
        std::int64_t end_1   = std::atoll(std::string(col[4]).c_str());
        if (start_1 <= 0 || end_1 < start_1) {
            ++status.records_skipped;
            continue;
        }
        std::int64_t start_0 = start_1 - 1;
        std::int64_t end_0   = end_1;        // GFF3 end is inclusive,
                                              // so half-open end is the
                                              // same numeric value.

        auto& bag = exons_by_tx[tid_it->second];
        bag.exons.emplace_back(start_0, end_0);
        // GFF3 columns: 0=seqid (chrom), 1=source, 2=type, 3=start, 4=end,
        // 5=score, 6=strand, 7=phase, 8=attributes.
        if (bag.chrom.empty()) bag.chrom = std::string(col[0]);
        if (bag.strand == '.') bag.strand = col[6].empty() ? '.' : col[6][0];
    }
    gzclose(fh);

    // Emit one AnchorRecord per exon plus synthesised ExonBoundary
    // records on the first exon-record per transcript.
    for (auto& [tid, bag] : exons_by_tx) {
        if (bag.exons.empty()) continue;

        // Sort exons by start so junction-pairing is correct
        std::sort(bag.exons.begin(), bag.exons.end());

        const auto meta_it = tx_meta.find(tid);
        const core::TranscriptKind kind = (meta_it != tx_meta.end())
            ? BiotypeToKind(meta_it->second.biotype)
            : core::TranscriptKind::Unknown;

        // Build ExonBoundary list (one per junction, n_exons - 1 entries).
        std::vector<ExonBoundary> boundaries;
        boundaries.reserve(bag.exons.size() > 1 ? bag.exons.size() - 1 : 0);
        std::uint32_t cumulative_tx_pos = 0;
        for (std::size_t i = 0; i + 1 < bag.exons.size(); ++i) {
            const auto [s_i, e_i] = bag.exons[i];
            cumulative_tx_pos += static_cast<std::uint32_t>(e_i - s_i);

            ExonBoundary b;
            b.pos_in_transcript    = cumulative_tx_pos;
            b.donor_genomic_pos    = static_cast<std::uint64_t>(e_i - 1);
            b.acceptor_genomic_pos = static_cast<std::uint64_t>(bag.exons[i + 1].first);
            // donor/acceptor motifs + scores filled later by SpliceSiteDb
            boundaries.push_back(std::move(b));
        }

        // Emit one AnchorRecord per exon. Transcript-level boundaries are
        // attached to ALL exons of that transcript so a single exon-anchor
        // hit can resolve the whole transcript architecture.
        for (std::size_t i = 0; i < bag.exons.size(); ++i) {
            AnchorRecord rec;
            rec.anchor_id = "GENCODE:" + tid + ":exon" + std::to_string(i + 1);
            rec.source = AnchorSource::Gencode;
            rec.kind = kind;
            rec.transcript_id = tid;
            if (meta_it != tx_meta.end()) {
                rec.host_gene_id = meta_it->second.gene_id;
            }
            rec.ref_chrom = bag.chrom;
            rec.ref_start = bag.exons[i].first;
            rec.ref_end   = bag.exons[i].second;
            rec.strand    = bag.strand;
            rec.exon_boundaries = boundaries;  // shared across exons of one tx

            // tag with gene name if known + biotype tag for downstream
            // filtering
            if (meta_it != tx_meta.end()) {
                if (!meta_it->second.gene_name.empty()) {
                    rec.tags.push_back("gene:" + meta_it->second.gene_name);
                }
                if (!meta_it->second.biotype.empty()) {
                    rec.tags.push_back("biotype:" + meta_it->second.biotype);
                }
            }

            AddAnchor(std::move(rec));
            ++status.records_loaded;
        }
    }

    // Final reindex so the per-chrom vectors are sorted by ref_start.
    Reindex();

    status.ok = true;
    return status;
}

// ---------------------------------------------------------------------------
// Stub implementations for sources we'll fill in subsequent blocks.
// ---------------------------------------------------------------------------

LoadStatus AnchorStore::LoadMane(const std::filesystem::path& /*tsv*/) {
    LoadStatus s;
    s.ok = true;  // empty no-op; MANE-loader lands in Block 9 surfaces.
    s.error = "MANE loader not yet implemented (Block 9 integration "
              "surfaces). Returning empty success so the pipeline keeps "
              "moving; tag GENCODE anchors manually if needed in interim.";
    return s;
}

LoadStatus AnchorStore::LoadImgtGeneDb(
    const std::filesystem::path& /*fasta*/) {
    LoadStatus s;
    s.ok = true;
    s.error = "IMGT/GENE-DB loader not yet implemented (Block 9).";
    return s;
}

LoadStatus AnchorStore::LoadPangenomeAnnotations(
    const std::filesystem::path& /*hprc_root*/,
    const std::vector<std::string>& /*sample_ids*/) {
    LoadStatus s;
    s.ok = true;
    s.error = "Pangenome per-Hap loader not yet implemented (Block 9).";
    return s;
}

LoadStatus AnchorStore::ImportBranchBubbles(
    const std::filesystem::path& /*bed*/) {
    LoadStatus s;
    s.ok = true;
    s.error = "BRANCH bubble import not yet implemented (Block 9).";
    return s;
}

}  // namespace llmap::anchor
