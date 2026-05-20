// LLmap — post-hoc IGH locus re-sort implementation.

#include "igh_locus/igh_resort.h"

#include <cstdio>
#include <optional>
#include <string>

namespace llmap::igh_locus {

namespace {

ParalogCall MakeParalogCall(const IghMatch& m) {
    ParalogCall pc;
    pc.inter_paralog.push_back({m.gene, 1.0f});
    pc.n_discriminating_psvs = static_cast<std::uint32_t>(m.n_distinct_exons);
    return pc;
}

}  // namespace

ResortStats ApplyResort(const IghAnchorCatalog& catalog,
                        std::vector<AlignmentRecord>& records,
                        const std::vector<std::string>& read_seqs,
                        const ResortOptions& opts) {
    ResortStats st;
    for (std::size_t i = 0; i < records.size(); ++i) {
        ++st.n_examined;
        AlignmentRecord& rec = records[i];
        if (!rec.primary.has_value()) continue;  // unmapped: not rescued in v1

        AlignmentHit& hit = *rec.primary;

        // Is this an IGH read? transcriptomic (target is an IGH copy/gene) or
        // genomic (alignment falls inside the IGHC region).
        std::optional<std::string> cur_gene =
            catalog.GeneOfTarget(hit.target_id);
        const bool genomic =
            !cur_gene.has_value() && opts.region.Contains(hit.target_id, hit.start);
        if (!cur_gene.has_value() && !genomic) continue;
        ++st.n_in_locus;

        if (i >= read_seqs.size() || read_seqs[i].empty()) continue;

        IghMatch m = catalog.Match(read_seqs[i]);
        if (!m.matched || m.ambiguous_gene ||
            m.n_distinct_exons < opts.min_exons) {
            continue;
        }
        ++st.n_matched;

        if (opts.set_paralog) {
            rec.paralog_assignment = MakeParalogCall(m);
            ++st.n_paralog_set;
        }

        if (!opts.relabel) continue;

        bool changed = false;
        if (cur_gene.has_value()) {
            // Transcriptomic: relabel target to the true copy when it differs.
            if (hit.target_id != m.copy_label) {
                hit.target_id = m.copy_label;
                hit.is_reverse = m.is_reverse;
                changed = true;
            }
        } else if (genomic && m.contig.has_value() && m.start.has_value()) {
            // Genomic: re-place onto the true copy's interval when the current
            // position does not already fall inside it.
            const std::uint64_t cs = *m.start;
            const std::uint64_t ce = m.end.value_or(cs);
            const bool same_contig = (hit.target_id == *m.contig);
            const bool inside = same_contig && hit.start >= cs && hit.start < ce;
            if (!inside) {
                hit.target_id = *m.contig;
                hit.start = cs;
                hit.end = ce;
                hit.is_reverse = m.is_reverse;
                changed = true;
            }
        }

        if (changed) {
            ++st.n_resorted;
            if (opts.verbose) {
                std::fprintf(stderr,
                    "[igh] resort %s -> %s (%s, %d exons)\n",
                    rec.read_id.c_str(), m.copy_label.c_str(),
                    m.gene.c_str(), m.n_distinct_exons);
            }
        }
    }

    if (opts.verbose) {
        std::fprintf(stderr,
            "[igh] examined=%zu in_locus=%zu matched=%zu resorted=%zu paralog_set=%zu\n",
            st.n_examined, st.n_in_locus, st.n_matched, st.n_resorted,
            st.n_paralog_set);
    }
    return st;
}

}  // namespace llmap::igh_locus
