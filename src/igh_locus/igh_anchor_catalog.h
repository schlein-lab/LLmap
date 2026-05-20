// LLmap — IGH paralog-specific exon anchor catalog.
//
// Each anchor is a full constant-region exon sequence (CH1/CH2/CH3/CH4, ~250-360
// bp) extracted from a specific genomic copy. Because a single SNP between two
// IGH paralogs is enough to break an exact match, exact-substring matching of a
// full exon is highly specific: a read that contains a copy's exact CH exon can
// only have come from that copy (or a CDS-identical sibling). This is the same
// alignment-free principle validated on HPRC iso-seq: minimap2 mis-assigns
// paralog reads, but exact exon anchors recover the true copy of origin.
//
// FASTA header convention (tokens split on '_'):
//   <GENE>_<COPY>_<HAP>_<EXON>[_<LEN>bp]   e.g. IGHG4_G056511_hap1_CH1_295bp
// The leading token is the isotype/gene; everything up to the exon token is the
// copy label (the re-sort target). An optional " loc=<contig>:<start>-<end>"
// suffix on the header line carries the copy's genomic interval for genomic-mode
// re-placement.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace llmap::igh_locus {

// One exon anchor from one genomic copy.
struct IghAnchor {
    std::string name;        // full header token (unique)
    std::string gene;        // isotype, e.g. "IGHG4"
    std::string copy_label;  // header minus exon suffix, e.g. "IGHG4_G056511_hap1"
    std::string exon;        // e.g. "CH1"
    std::string seq;         // forward exon sequence (uppercase)
    bool has_locus{false};
    std::string contig;
    std::uint64_t start{0};
    std::uint64_t end{0};
};

// Result of matching one read against the catalog.
struct IghMatch {
    bool matched{false};
    std::string gene;           // best-supported isotype
    std::string copy_label;     // best-supported copy
    int n_anchor_hits{0};       // total anchor sequences hit
    int n_distinct_exons{0};    // distinct exons of the winning copy
    bool ambiguous_gene{false}; // >1 gene tied on distinct-exon support
    bool is_reverse{false};     // winning copy matched on the reverse strand
    std::vector<std::string> matched_exons;  // exons of the winning copy
    std::optional<std::string> contig;       // winning copy genomic contig
    std::optional<std::uint64_t> start;      // winning copy genomic start
    std::optional<std::uint64_t> end;
};

// Returns true if a token is a recognised IGH constant-region isotype.
[[nodiscard]] bool IsIghGene(std::string_view token);

class IghAnchorCatalog {
public:
    // Load anchors from a FASTA. max_mismatch>0 enables a seeded Hamming search
    // (tolerates sequencing error in e.g. minimap2 reads); 0 = exact only.
    // Returns nullopt on read/parse failure or if no IGH anchors were found.
    [[nodiscard]] static std::optional<IghAnchorCatalog> LoadFasta(
        const std::string& path, int max_mismatch = 0);

    // Best paralog/copy supported by exact (or <=max_mismatch) exon anchors.
    [[nodiscard]] IghMatch Match(std::string_view read_seq) const;

    // Gene of a reference target name, if that target is a known IGH copy.
    // Recognises exact copy labels and any target whose leading token is an
    // IGH isotype (covers transcript references named e.g. "IGHG4_..."/"IGHG4").
    [[nodiscard]] std::optional<std::string> GeneOfTarget(
        std::string_view target_id) const;

    [[nodiscard]] std::size_t size() const { return anchors_.size(); }
    [[nodiscard]] const std::vector<IghAnchor>& anchors() const {
        return anchors_;
    }
    [[nodiscard]] int max_mismatch() const { return max_mismatch_; }

private:
    struct UniqueSeq {
        std::string fwd;
        std::string rev;
        std::vector<std::size_t> anchor_idx;  // anchors sharing this sequence
    };

    [[nodiscard]] bool ContainsSeq(std::string_view read,
                                   const std::string& probe) const;

    std::vector<IghAnchor> anchors_;
    std::vector<UniqueSeq> unique_;
    std::unordered_map<std::string, std::string> copy_to_gene_;
    int max_mismatch_{0};
};

}  // namespace llmap::igh_locus
