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
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"

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

    auto start_time = std::chrono::steady_clock::now();

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

    if (args.verbose) {
        std::fprintf(stderr, "Building minimizer index (k=%d, w=%d)\n",
                     args.kmer_size, args.window_size);
    }

    classical::MinimizerConfig mini_cfg;
    mini_cfg.k = args.kmer_size;
    mini_cfg.w = args.window_size;

    classical::MinimizerIndex::Builder idx_builder(mini_cfg);
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        idx_builder.AddSequence(ref_names[i], ref_seqs[i]);
    }
    auto index = idx_builder.Build();

    if (args.verbose) {
        std::fprintf(stderr, "Index built: %zu minimizers\n", index->Size());
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

    if (args.verbose) {
        std::fprintf(stderr, "Loading reads: %s\n", args.reads.c_str());
    }

    auto read_reader = io::FastqReader::Open(args.reads);
    if (!read_reader) {
        std::fprintf(stderr, "Error: failed to open reads file: %s\n", args.reads.c_str());
        return 1;
    }

    std::vector<std::string> read_names;
    std::vector<std::string> read_seqs;
    std::vector<std::uint32_t> read_lens;

    while (read_reader->HasMore()) {
        auto record = read_reader->Next();
        if (record && record->IsValid()) {
            read_names.push_back(record->id);
            read_seqs.push_back(record->sequence);
            read_lens.push_back(static_cast<std::uint32_t>(record->sequence.size()));
        }
    }

    if (read_names.empty()) {
        std::fprintf(stderr, "Error: no reads in input\n");
        return 1;
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loaded %zu reads\n", read_names.size());
    }

    if (args.verbose) {
        std::fprintf(stderr, "Aligning...\n");
    }

    auto align_start = std::chrono::steady_clock::now();
    auto results = pipeline.AlignReads(read_names, read_seqs);
    auto align_end = std::chrono::steady_clock::now();

    float align_time_ms = std::chrono::duration<float, std::milli>(
        align_end - align_start).count();

    std::vector<AlignmentRecord> records;
    records.reserve(results.size());

    std::size_t n_mapped = 0;
    std::size_t n_unmapped = 0;

    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];

        if (result.HasAlignment()) {
            const auto* primary = result.PrimaryAlignment();
            if (primary) {
                records.push_back(ConvertAlignment(*primary, read_lens[i]));
                n_mapped++;
            }
        } else {
            records.push_back(make_unmapped(
                result.query_name, read_lens[i], RejectionReason::NoSeeds));
            n_unmapped++;
        }
    }

    if (args.verbose) {
        std::fprintf(stderr, "Writing output: %s\n", args.output.c_str());
    }

    std::vector<output::ReferenceSequence> ref_info;
    ref_info.reserve(ref_names.size());
    for (std::size_t i = 0; i < ref_names.size(); ++i) {
        ref_info.push_back({ref_names[i], ref_seqs[i].size()});
    }

    output::BamWriterConfig bam_cfg;
    bam_cfg.format = args.use_bam ? output::BamOutputFormat::BAM
                                   : output::BamOutputFormat::SAM;
    bam_cfg.include_wavecollapse_tags = false;
    bam_cfg.include_paralog_tags = false;

    auto bam_writer = output::BamWriter::Create(args.output, ref_info, bam_cfg);
    if (!bam_writer) {
        std::fprintf(stderr, "Error: failed to create output file: %s\n",
                     args.output.c_str());
        return 1;
    }

    if (!bam_writer->WriteBatch(records)) {
        std::fprintf(stderr, "Error: failed to write alignments: %s\n",
                     bam_writer->LastError().c_str());
        return 1;
    }

    bam_writer->Close();

    if (!args.parquet_output.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Writing Parquet: %s\n", args.parquet_output.c_str());
        }

        output::ParquetWriterConfig parquet_cfg;
        parquet_cfg.min_probability = 0.0f;

        auto parquet_writer = output::ParquetWriter::Create(
            args.parquet_output, parquet_cfg);
        if (!parquet_writer) {
            std::fprintf(stderr, "Warning: failed to create Parquet output\n");
        } else {
            parquet_writer->WriteBatch(records);
            parquet_writer->Close();
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    float mapping_rate = static_cast<float>(n_mapped) /
        static_cast<float>(read_names.size());

    std::printf("Alignment complete:\n");
    std::printf("  Input reads:    %zu\n", read_names.size());
    std::printf("  Mapped:         %zu (%.1f%%)\n",
                n_mapped, 100.0f * mapping_rate);
    std::printf("  Unmapped:       %zu\n", n_unmapped);
    std::printf("  Align time:     %.2f s\n", align_time_ms / 1000.0f);
    std::printf("  Total time:     %.2f s\n", total_time_ms / 1000.0f);
    std::printf("  Throughput:     %.1f reads/s\n",
                read_names.size() / (align_time_ms / 1000.0f));
    std::printf("  Output:         %s\n", args.output.c_str());
    if (!args.parquet_output.empty()) {
        std::printf("  Parquet:        %s\n", args.parquet_output.c_str());
    }

    if (args.enable_llm && mapping_rate < args.llm_threshold) {
        if (args.verbose) {
            std::fprintf(stderr, "\nMapping rate %.1f%% below threshold %.1f%%, "
                         "running LLM diagnostics...\n",
                         100.0f * mapping_rate, 100.0f * args.llm_threshold);
        }

        auto agent = CreatePipelineAgent(args, args.verbose);
        if (agent) {
            float avg_identity = pipeline.Stats().avg_identity;
            RunLlmDiagnostics(*agent, args, results, n_mapped, n_unmapped, avg_identity);
        }
    } else if (args.enable_llm && args.verbose) {
        std::fprintf(stderr, "Mapping rate %.1f%% meets threshold, "
                     "skipping LLM diagnostics\n", 100.0f * mapping_rate);
    }

    return 0;
}

}  // namespace llmap::cli
