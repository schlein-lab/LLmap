// LLmap — `llmap align` CLI command.

#include "cli/commands.h"
#include "cli/cmd_align_internal.h"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "classical/classical_pipeline.h"
#include "core/alignment_record.h"
#include "core/thread_pool.h"
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"
#include "psv/psv_catalog.h"

namespace llmap::cli {

namespace {

std::string CigarToString(const std::vector<classical::CigarElement>& cigar) {
    std::ostringstream oss;
    for (const auto& elem : cigar) {
        oss << elem.ToString();
    }
    return oss.str();
}

AlignmentRecord ConvertAlignment(const classical::ClassicalAlignment& aln,
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

}  // namespace

int run_align(int argc, char** argv) {
    using namespace align_internal;

    AlignArgs args;
    if (!ParseAlignArgs(argc, argv, args)) {
        PrintAlignUsage();
        return 1;
    }

    if (args.help) {
        PrintAlignUsage();
        return 0;
    }

    if (args.reads.empty()) {
        std::fprintf(stderr, "Error: --reads is required\n");
        PrintAlignUsage();
        return 1;
    }

    if (args.reference.empty()) {
        std::fprintf(stderr, "Error: --reference is required\n");
        PrintAlignUsage();
        return 1;
    }

    if (args.output.empty()) {
        std::fprintf(stderr, "Error: --output is required\n");
        PrintAlignUsage();
        return 1;
    }

    if (!std::filesystem::exists(args.reads)) {
        std::fprintf(stderr, "Error: reads file not found: %s\n", args.reads.c_str());
        return 1;
    }

    if (!std::filesystem::exists(args.reference)) {
        std::fprintf(stderr, "Error: reference file not found: %s\n", args.reference.c_str());
        return 1;
    }

    if (!args.index.empty() && !std::filesystem::exists(args.index)) {
        std::fprintf(stderr, "Error: index file not found: %s\n", args.index.c_str());
        return 1;
    }

    auto start_time = std::chrono::steady_clock::now();

    // Load PSV catalog if specified
    std::optional<psv::PsvCatalog> psv_catalog;
    if (!args.psv_catalog.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Loading PSV catalog: %s\n", args.psv_catalog.c_str());
        }
        psv_catalog = LoadPsvCatalog(args);
        if (!psv_catalog) {
            return 1;
        }
        if (args.verbose) {
            std::fprintf(stderr, "Loaded %zu PSV sites\n", psv_catalog->Size());
        }
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loading reference: %s\n", args.reference.c_str());
    }

    io::FastaReader ref_reader(args.reference);
    std::vector<std::string> ref_names;
    std::vector<std::string> ref_seqs;

    while (ref_reader.HasMore()) {
        auto record = ref_reader.Next();
        if (record.IsValid()) {
            ref_names.push_back(record.name);
            ref_seqs.push_back(record.sequence);
        }
    }

    if (ref_names.empty()) {
        std::fprintf(stderr, "Error: no sequences in reference\n");
        return 1;
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loaded %zu reference sequences\n", ref_names.size());
    }

    std::unique_ptr<classical::MinimizerIndex> index;
    classical::MinimizerConfig mini_cfg;

    if (!args.index.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Loading pre-built index: %s\n", args.index.c_str());
        }
        index = classical::MinimizerIndex::Load(args.index);
        if (!index) {
            std::fprintf(stderr, "Error: failed to load index: %s\n", args.index.c_str());
            return 1;
        }
        mini_cfg = index->GetConfig();
        if (args.verbose) {
            std::fprintf(stderr, "Index loaded: %zu minimizers (k=%d, w=%d)\n",
                         index->Size(), mini_cfg.k, mini_cfg.w);
        }
    } else {
        if (args.verbose) {
            std::fprintf(stderr, "Building minimizer index (k=%d, w=%d)\n",
                         args.kmer_size, args.window_size);
        }
        mini_cfg.k = args.kmer_size;
        mini_cfg.w = args.window_size;

        classical::MinimizerIndex::Builder idx_builder(mini_cfg);
        for (std::size_t i = 0; i < ref_names.size(); ++i) {
            idx_builder.AddSequence(ref_names[i], ref_seqs[i]);
        }
        index = idx_builder.Build();

        if (args.verbose) {
            std::fprintf(stderr, "Index built: %zu minimizers\n", index->Size());
        }
    }

    classical::ClassicalPipelineConfig pipe_cfg;
    pipe_cfg.minimizer_config = mini_cfg;
    pipe_cfg.chain_config.min_chain_score = args.min_chain;
    pipe_cfg.extension_config.match_score = 2;
    pipe_cfg.extension_config.mismatch_penalty = 4;
    pipe_cfg.extension_config.gap_open = 4;
    pipe_cfg.extension_config.gap_extend = 2;
    pipe_cfg.min_identity = args.min_identity;
    pipe_cfg.max_chains_to_extend = args.max_chains;

    classical::ClassicalPipeline pipeline(pipe_cfg);
    pipeline.SetIndex(std::move(index));
    pipeline.SetReferenceSequences(ref_seqs);

    if (args.verbose) {
        std::fprintf(stderr, "Loading reads: %s\n", args.reads.c_str());
    }

    auto read_reader = io::FastqReader::Open(args.reads);
    if (!read_reader) {
        std::fprintf(stderr, "Error: failed to open reads file: %s\n", args.reads.c_str());
        return 1;
    }

    // Open output writers up-front so we can stream batches through.
    std::vector<output::ReferenceSequence> ref_info;
    ref_info.reserve(ref_names.size());
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        ref_info.push_back({ref_names[i], ref_seqs[i].size()});
    }

    output::BamWriterConfig bam_cfg;
    bam_cfg.format = args.use_bam ? output::BamOutputFormat::BAM
                                   : output::BamOutputFormat::SAM;
    bam_cfg.include_wavecollapse_tags = false;
    bam_cfg.include_paralog_tags = psv_catalog.has_value();

    auto bam_writer = output::BamWriter::Create(args.output, ref_info, bam_cfg);
    if (!bam_writer) {
        std::fprintf(stderr, "Error: failed to create output file: %s\n",
                     args.output.c_str());
        return 1;
    }

    std::unique_ptr<output::ParquetWriter> parquet_writer;
    if (!args.parquet_output.empty()) {
        output::ParquetWriterConfig parquet_cfg;
        parquet_cfg.min_probability = 0.0f;
        parquet_writer = output::ParquetWriter::Create(
            args.parquet_output, parquet_cfg);
        if (!parquet_writer) {
            std::fprintf(stderr, "Warning: failed to create Parquet output\n");
        }
    }

    // Single thread pool reused across batches.
    std::unique_ptr<llmap::core::ThreadPool> thread_pool;
    if (args.threads > 1) {
        thread_pool = std::make_unique<llmap::core::ThreadPool>(
            static_cast<size_t>(args.threads));
    }

    // Stream reads in batches to bound RAM independent of input size.
    // 50k * ~15 kb HiFi = ~750 MB peak; safe even for tens-of-millions-read WGS.
    constexpr std::size_t kBatchSize = 50000;

    classical::ClassicalPipelineStats agg_stats;
    std::size_t total_reads = 0;
    std::size_t n_mapped = 0;
    std::size_t n_unmapped = 0;
    double identity_sum_weighted = 0.0;  // \Sigma (batch_avg * batch_aligned)

    auto align_start = std::chrono::steady_clock::now();

    while (read_reader->HasMore()) {
        std::vector<std::string> batch_names;
        std::vector<std::string> batch_seqs;
        std::vector<std::uint32_t> batch_lens;
        batch_names.reserve(kBatchSize);
        batch_seqs.reserve(kBatchSize);
        batch_lens.reserve(kBatchSize);

        while (batch_names.size() < kBatchSize && read_reader->HasMore()) {
            auto record = read_reader->Next();
            if (record && record->IsValid()) {
                batch_names.push_back(record->id);
                batch_seqs.push_back(record->sequence);
                batch_lens.push_back(
                    static_cast<std::uint32_t>(record->sequence.size()));
            }
        }
        if (batch_names.empty()) break;

        std::vector<llmap::classical::ReadAlignmentResult> results;
        if (thread_pool) {
            results = pipeline.AlignReadsParallel(
                batch_names, batch_seqs, *thread_pool);
        } else {
            results = pipeline.AlignReads(batch_names, batch_seqs);
        }

        const auto& bs = pipeline.Stats();
        agg_stats.total_hits += bs.total_hits;
        agg_stats.total_chains += bs.total_chains;
        agg_stats.total_extensions += bs.total_extensions;
        agg_stats.alignments_filtered_by_identity +=
            bs.alignments_filtered_by_identity;
        agg_stats.alignments_filtered_by_length +=
            bs.alignments_filtered_by_length;
        agg_stats.seeding_time_ms += bs.seeding_time_ms;
        agg_stats.chaining_time_ms += bs.chaining_time_ms;
        agg_stats.extension_time_ms += bs.extension_time_ms;
        agg_stats.reads_aligned += bs.reads_aligned;
        agg_stats.reads_unmapped += bs.reads_unmapped;
        identity_sum_weighted += static_cast<double>(bs.avg_identity) *
                                 static_cast<double>(bs.reads_aligned);

        std::vector<AlignmentRecord> records;
        records.reserve(results.size());
        for (std::size_t i = 0; i < results.size(); ++i) {
            const auto& result = results[i];
            if (result.HasAlignment()) {
                const auto* primary = result.PrimaryAlignment();
                if (primary) {
                    records.push_back(ConvertAlignment(*primary, batch_lens[i]));
                    ++n_mapped;
                }
            } else {
                records.push_back(make_unmapped(
                    result.query_name, batch_lens[i],
                    RejectionReason::NoSeeds));
                ++n_unmapped;
            }
        }

        if (psv_catalog) {
            ApplyPsvAssignments(
                *psv_catalog, args, records, batch_seqs, args.verbose);
        }

        if (!bam_writer->WriteBatch(records)) {
            std::fprintf(stderr, "Error: failed to write alignments: %s\n",
                         bam_writer->LastError().c_str());
            return 1;
        }
        if (parquet_writer) {
            parquet_writer->WriteBatch(records);
        }

        total_reads += batch_names.size();
        if (args.verbose) {
            std::fprintf(stderr,
                "  Processed %zu reads (mapped %zu / unmapped %zu)\n",
                total_reads, n_mapped, n_unmapped);
        }
    }

    auto align_end = std::chrono::steady_clock::now();
    float align_time_ms = std::chrono::duration<float, std::milli>(
        align_end - align_start).count();

    bam_writer->Close();
    if (parquet_writer) {
        parquet_writer->Close();
    }

    if (agg_stats.reads_aligned > 0) {
        agg_stats.avg_identity =
            static_cast<float>(identity_sum_weighted /
                               static_cast<double>(agg_stats.reads_aligned));
    }

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    float mapping_rate = total_reads > 0
        ? static_cast<float>(n_mapped) / static_cast<float>(total_reads)
        : 0.0f;

    std::printf("Alignment complete:\n");
    std::printf("  Input reads:    %zu\n", total_reads);
    std::printf("  Mapped:         %zu (%.1f%%)\n",
                n_mapped, 100.0f * mapping_rate);
    std::printf("  Unmapped:       %zu\n", n_unmapped);
    std::printf("  Pipeline drop breakdown:\n");
    std::printf("    minimizer hits found: %zu (avg %.1f/read)\n",
                agg_stats.total_hits,
                static_cast<float>(agg_stats.total_hits) /
                    std::max<size_t>(1, total_reads));
    std::printf("    chains formed:        %zu (avg %.2f/read)\n",
                agg_stats.total_chains,
                static_cast<float>(agg_stats.total_chains) /
                    std::max<size_t>(1, total_reads));
    std::printf("    extensions accepted:  %zu\n", agg_stats.total_extensions);
    std::printf("    rejected: identity:   %zu\n",
                agg_stats.alignments_filtered_by_identity);
    std::printf("    rejected: length:     %zu\n",
                agg_stats.alignments_filtered_by_length);
    std::printf("    avg identity:         %.3f\n", agg_stats.avg_identity);
    std::printf("  Phase breakdown (sum across threads):\n");
    std::printf("    seeding:     %.2f s\n", agg_stats.seeding_time_ms / 1000.0f);
    std::printf("    chaining:    %.2f s\n", agg_stats.chaining_time_ms / 1000.0f);
    std::printf("    extension:   %.2f s\n", agg_stats.extension_time_ms / 1000.0f);
    std::printf("  Align time:     %.2f s\n", align_time_ms / 1000.0f);
    std::printf("  Total time:     %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Throughput:     %.1f reads/s\n",
                total_reads / (align_time_ms / 1000.0f));
    std::printf("  Output:         %s\n", args.output.c_str());
    if (!args.parquet_output.empty()) {
        std::printf("  Parquet:        %s\n", args.parquet_output.c_str());
    }

    // --classical-only overrides --llm (skip probabilistic framework)
    bool llm_enabled = args.enable_llm && !args.classical_only;

    if (args.enable_llm && args.classical_only && args.verbose) {
        std::fprintf(stderr, "Note: --llm ignored due to --classical-only mode\n");
    }

    if (llm_enabled && mapping_rate < args.llm_threshold) {
        if (args.verbose) {
            std::fprintf(stderr, "\nMapping rate %.1f%% below threshold %.1f%%, "
                         "running LLM diagnostics...\n",
                         100.0f * mapping_rate, 100.0f * args.llm_threshold);
        }

        auto agent = CreatePipelineAgent(args, args.verbose);
        if (agent) {
            float avg_identity = agg_stats.avg_identity;
            // Streaming mode discards per-batch ReadAlignmentResults; LLM
            // diagnostics only see aggregate counters. Re-running with a
            // smaller subsample would be needed to inspect individual reads.
            std::vector<llmap::classical::ReadAlignmentResult> empty_results;
            RunLlmDiagnostics(*agent, args, empty_results,
                              n_mapped, n_unmapped, avg_identity);
        }
    } else if (llm_enabled && args.verbose) {
        std::fprintf(stderr, "Mapping rate %.1f%% meets threshold, "
                     "skipping LLM diagnostics\n", 100.0f * mapping_rate);
    }

    return 0;
}

}  // namespace llmap::cli
