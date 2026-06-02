// LLmap — Misc-source anchor loaders.
//
// Replaces the stubs in anchor_loader_gencode.cpp for non-GENCODE
// sources: IMGT/GENE-DB FASTA, MANE TSV, BRANCH bubble BED, HPRC
// per-haplotype GFF. Lives in its own TU to keep the GENCODE loader
// under the 400-LOC modular cap and to let each source's parser
// evolve independently.

#include "anchor/anchor_store.h"

#include "core/transcript_kind.h"

#include <fstream>
#include <sstream>
#include <string>

namespace llmap::anchor {

namespace {

std::vector<std::string> SplitTab(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, '\t')) out.push_back(item);
    return out;
}

}  // namespace

// ===========================================================================
// LoadImgtGeneDb — IMGT/GENE-DB FASTA.
//
// Each record is a germline V/D/J/C segment for an Ig/TR gene. We use
// the FASTA header (everything after '>') as the anchor_id with an
// "IMGT:" prefix and tag the locus from the segment-name prefix:
//
//   IGHV1-2*01 → segment_type='V', locus='IGH'
//   IGHC*01    → segment_type='C', locus='IGH'
//
// Tagging strategy:
//   tags: ["IGH", "V_gene"]   for IGHV records
//   tags: ["IGH", "C_gene"]   for IGHC records
//   etc.
//
// The kind stays MatureMrna at the segment level — the classifier
// (Block 2.5) refines into SterileGermline when a sterile-transcript
// pattern is observed at read time.
// ===========================================================================

LoadStatus AnchorStore::LoadImgtGeneDb(const std::filesystem::path& fasta) {
    LoadStatus status;
    std::ifstream in(fasta);
    if (!in) {
        status.error = "IMGT FASTA not found: " + fasta.string();
        return status;
    }

    std::string line, current_id, current_seq;

    auto flush = [&]() {
        if (current_id.empty()) return;
        AnchorRecord r;
        r.anchor_id = "IMGT:" + current_id;
        r.source = AnchorSource::Imgt_GeneDb;
        r.kind = core::TranscriptKind::MatureMrna;
        r.sequence = current_seq;
        // tag from prefix: IGH / IGK / IGL / TRA / TRB / TRD / TRG
        const std::string locus = current_id.substr(0,
            current_id.size() >= 3 ? 3 : current_id.size());
        if (locus == "IGH" || locus == "IGK" || locus == "IGL"
            || locus == "TRA" || locus == "TRB"
            || locus == "TRD" || locus == "TRG") {
            r.tags.push_back(locus);
            // segment type from char 4 of the locus name when present
            if (current_id.size() >= 4) {
                const char seg = current_id[3];
                switch (seg) {
                    case 'V': r.tags.push_back("V_gene"); break;
                    case 'D': r.tags.push_back("D_gene"); break;
                    case 'J': r.tags.push_back("J_gene"); break;
                    case 'C': r.tags.push_back("C_gene"); break;
                    default: break;
                }
            }
        }
        AddAnchor(std::move(r));
        ++status.records_loaded;
    };

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            current_id.clear();
            current_seq.clear();
            // Take everything between '>' and the first whitespace.
            const auto end = line.find_first_of(" \t");
            current_id = line.substr(1,
                end == std::string::npos ? std::string::npos : end - 1);
        } else {
            current_seq += line;
        }
    }
    flush();

    Reindex();
    status.ok = true;
    return status;
}

// ===========================================================================
// LoadMane — MANE Select TSV.
//
// MANE Select consists of cross-references to existing GENCODE/RefSeq
// transcripts. We don't add new AnchorRecords here; instead we walk the
// existing GENCODE anchors and tag the MANE-selected ones with
// "mane_select" so downstream filters can pick them up.
//
// Expected TSV columns: GENCODE_id, RefSeq_id, gene_name (any layout
// containing a GENCODE-style ENST id will work).
// ===========================================================================

LoadStatus AnchorStore::LoadMane(const std::filesystem::path& tsv) {
    LoadStatus status;
    std::ifstream in(tsv);
    if (!in) {
        status.error = "MANE TSV not found: " + tsv.string();
        return status;
    }

    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; continue; }      // header
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        for (const auto& c : cols) {
            if (c.size() >= 4 && c.substr(0, 4) == "ENST") {
                // Walk existing GENCODE anchors with that transcript_id.
                auto idxs = ByTranscriptId(c);
                for (auto idx : idxs) {
                    auto& a = const_cast<AnchorRecord&>(anchors_[idx]);
                    a.tags.emplace_back("mane_select");
                    ++status.records_loaded;
                }
            }
        }
    }
    status.ok = true;
    return status;
}

// ===========================================================================
// ImportBranchBubbles — BRANCH bubble BED (see [[branch_tooling_state]]).
//
// BED columns (BRANCH 8-col emit):
//   1. chrom
//   2. start
//   3. end
//   4. bubble_id
//   5. confidence
//   6. per_alt_supports (comma-sep)
//   7. min_vaf
//   8. total_support
//
// We emit one AnchorRecord per bubble with source=Branch_Bubble and
// kind=NovelUnclassified — the classifier (Block 2.5) refines later if
// the bubble structure resembles a known biotype.
// ===========================================================================

LoadStatus AnchorStore::ImportBranchBubbles(
    const std::filesystem::path& bed) {

    LoadStatus status;
    std::ifstream in(bed);
    if (!in) {
        status.error = "BRANCH bubble BED not found: " + bed.string();
        return status;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto cols = SplitTab(line);
        if (cols.size() < 4) continue;

        try {
            const std::string& chrom = cols[0];
            const auto start = std::stoll(cols[1]);
            const auto end   = std::stoll(cols[2]);
            const std::string& bubble_id = cols[3];

            AnchorRecord r;
            r.anchor_id = "BRANCH:" + bubble_id;
            r.source = AnchorSource::Branch_Bubble;
            r.kind = core::TranscriptKind::NovelUnclassified;
            r.ref_chrom = chrom;
            r.ref_start = start;
            r.ref_end = end;
            r.strand = '+';
            r.tags.push_back("branch_bubble");
            if (cols.size() >= 7) {
                r.tags.push_back("vaf:" + cols[6]);
            }
            AddAnchor(std::move(r));
            ++status.records_loaded;
        } catch (...) {
            ++status.records_skipped;
        }
    }

    Reindex();
    status.ok = true;
    return status;
}

// ===========================================================================
// LoadPangenomeAnnotations — HPRC R2 per-Hap annotation GFFs.
//
// Walks <hprc_root>/<sample_id>/annotation/*.gff3{,.gz} for each
// sample in `sample_ids` and pulls exon records into anchors tagged
// with the sample/hap id. We reuse the GENCODE GFF3 parse logic via
// LoadGencodeGff with an empty reference — the parser tolerates that
// (sequences left empty until k-mer-index pass).
// ===========================================================================

LoadStatus AnchorStore::LoadPangenomeAnnotations(
    const std::filesystem::path& hprc_root,
    const std::vector<std::string>& sample_ids) {

    LoadStatus total;
    total.ok = true;
    for (const auto& sample : sample_ids) {
        const auto annot_dir = hprc_root / sample / "annotation";
        if (!std::filesystem::exists(annot_dir)) {
            ++total.records_skipped;
            continue;
        }
        for (const auto& entry :
             std::filesystem::directory_iterator(annot_dir)) {
            if (!entry.is_regular_file()) continue;
            const auto path = entry.path();
            const auto ext = path.extension().string();
            if (ext != ".gff3" && ext != ".gz") continue;

            // Capture the pre-load size so we can tag the new records.
            const std::size_t pre = anchors_.size();
            auto s = LoadGencodeGff(path, /*ref_fa=*/"",
                                     /*with_sequence=*/false);
            if (!s.ok) {
                ++total.records_skipped;
                continue;
            }
            total.records_loaded += s.records_loaded;
            // Tag the newly inserted anchors with sample+source provenance.
            for (std::size_t i = pre; i < anchors_.size(); ++i) {
                auto& a = anchors_[i];
                a.source = AnchorSource::Pangenome_PerHap;
                a.tags.push_back("PANGEN_" + sample);
            }
        }
    }
    Reindex();
    return total;
}

}  // namespace llmap::anchor
