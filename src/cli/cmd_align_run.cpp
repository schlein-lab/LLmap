// LLmap — cmd_align batch processing logic.
// Contains the main alignment loop and batch streaming.

#include "cli/cmd_align_internal.h"

#include <chrono>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "classical/classical_pipeline.h"
#include "core/alignment_record.h"
#include "core/thread_pool.h"
#include "igh_locus/igh_anchor_catalog.h"
#include "igh_locus/igh_resort.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"
#include "psv/psv_catalog.h"

namespace llmap::cli::align_internal {

namespace {

std::string CigarToString(const std::vector<classical::CigarElement>& cigar) {
    std::ostringstream oss;
    for (const auto& elem : cigar) {
        oss << elem.ToString();
    }
    return oss.str();
}

}  // namespace

AlignmentRecord ConvertClassicalAlignment(
    const classical::ClassicalAlignment& aln,
    std::uint32_t read_len) {
    AlignmentHit hit;
    hit.target_id = aln.ref_name;
    hit.start = static_cast<std::uint64_t>(aln.ref_start);
    hit.end = static_cast<std::uint64_t>(aln.ref_end);
    hit.cigar = CigarString{CigarToString(aln.cigar)};
    hit.score = aln.score;
    hit.nm = static_cast<std::uint32_t>((1.0f - aln.identity) * aln.AlignedBases());
    hit.is_reverse = !aln.is_forward;

    return make_mapped(aln.query_name, read_len, std::move(hit));
}

BatchAlignResult RunAlignmentBatches(
    classical::ClassicalPipeline& pipeline,
    io::FastqReader& read_reader,
    output::BamWriter& bam_writer,
    output::ParquetWriter* parquet_writer,
    core::ThreadPool* thread_pool,
    const std::optional<psv::PsvCatalog>& psv_catalog,
    const AlignArgs& args) {

    constexpr std::size_t kBatchSize = 50000;

    BatchAlignResult result;

    // Load the IGH paralog anchor catalog once (post-hoc re-sort stage). The
    // stage is ON by default but only acts when an anchor FASTA is supplied.
    std::optional<igh_locus::IghAnchorCatalog> igh_catalog;
    if (args.enable_igh_locus && !args.igh_anchors.empty()) {
        igh_catalog = igh_locus::IghAnchorCatalog::LoadFasta(
            args.igh_anchors, args.igh_max_mismatch);
        if (igh_catalog) {
            std::fprintf(stderr,
                "[igh] post-hoc re-sort enabled: %zu anchors from %s (max_mismatch=%d)\n",
                igh_catalog->size(), args.igh_anchors.c_str(),
                args.igh_max_mismatch);
        } else {
            std::fprintf(stderr,
                "[igh] warning: could not load anchors from %s; skipping re-sort\n",
                args.igh_anchors.c_str());
        }
    } else if (args.enable_igh_locus && args.verbose) {
        std::fprintf(stderr,
            "[igh] re-sort enabled but no --igh-anchors given; stage is a no-op\n");
    }

    while (read_reader.HasMore()) {
        std::vector<std::string> batch_names;
        std::vector<std::string> batch_seqs;
        std::vector<std::uint32_t> batch_lens;
        batch_names.reserve(kBatchSize);
        batch_seqs.reserve(kBatchSize);
        batch_lens.reserve(kBatchSize);

        while (batch_names.size() < kBatchSize && read_reader.HasMore()) {
            auto record = read_reader.Next();
            if (record && record->IsValid()) {
                batch_names.push_back(record->id);
                batch_seqs.push_back(record->sequence);
                batch_lens.push_back(
                    static_cast<std::uint32_t>(record->sequence.size()));
            }
        }
        if (batch_names.empty()) break;

        std::vector<classical::ReadAlignmentResult> results;
        if (thread_pool) {
            results = pipeline.AlignReadsParallel(
                batch_names, batch_seqs, *thread_pool);
        } else {
            results = pipeline.AlignReads(batch_names, batch_seqs);
        }

        const auto& bs = pipeline.Stats();
        result.agg_stats.total_hits += bs.total_hits;
        result.agg_stats.total_chains += bs.total_chains;
        result.agg_stats.total_extensions += bs.total_extensions;
        result.agg_stats.alignments_filtered_by_identity +=
            bs.alignments_filtered_by_identity;
        result.agg_stats.alignments_filtered_by_length +=
            bs.alignments_filtered_by_length;
        result.agg_stats.seeding_time_ms += bs.seeding_time_ms;
        result.agg_stats.chaining_time_ms += bs.chaining_time_ms;
        result.agg_stats.extension_time_ms += bs.extension_time_ms;
        result.agg_stats.reads_aligned += bs.reads_aligned;
        result.agg_stats.reads_unmapped += bs.reads_unmapped;
        result.identity_sum_weighted += static_cast<double>(bs.avg_identity) *
                                        static_cast<double>(bs.reads_aligned);

        std::vector<AlignmentRecord> records;
        records.reserve(results.size());
        for (std::size_t i = 0; i < results.size(); ++i) {
            const auto& res = results[i];
            if (res.HasAlignment()) {
                const auto* primary = res.PrimaryAlignment();
                if (primary) {
                    records.push_back(ConvertClassicalAlignment(*primary, batch_lens[i]));
                    ++result.n_mapped;
                }
            } else {
                records.push_back(make_unmapped(
                    res.query_name, batch_lens[i],
                    RejectionReason::NoSeeds));
                ++result.n_unmapped;
            }
        }

        if (psv_catalog) {
            ApplyPsvAssignments(
                *psv_catalog, args, records, batch_seqs, args.verbose);
        }

        if (igh_catalog) {
            igh_locus::ResortOptions ropts;
            ropts.verbose = args.verbose;
            igh_locus::ResortStats rs =
                igh_locus::ApplyResort(*igh_catalog, records, batch_seqs, ropts);
            result.igh_resorted += rs.n_resorted;
            result.igh_paralog_set += rs.n_paralog_set;
        }

        if (!bam_writer.WriteBatch(records)) {
            std::fprintf(stderr, "Error: failed to write alignments: %s\n",
                         bam_writer.LastError().c_str());
            result.error = true;
            return result;
        }
        if (parquet_writer) {
            parquet_writer->WriteBatch(records);
        }

        result.total_reads += batch_names.size();
        if (args.verbose) {
            std::fprintf(stderr,
                "  Processed %zu reads (mapped %zu / unmapped %zu)\n",
                result.total_reads, result.n_mapped, result.n_unmapped);
        }
    }

    return result;
}

void FinalizeAlignStats(BatchAlignResult& result) {
    if (result.agg_stats.reads_aligned > 0) {
        result.agg_stats.avg_identity =
            static_cast<float>(result.identity_sum_weighted /
                               static_cast<double>(result.agg_stats.reads_aligned));
    }
}

}  // namespace llmap::cli::align_internal
