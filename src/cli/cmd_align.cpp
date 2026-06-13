// LLmap — `llmap align` CLI command.

#include "cli/commands.h"
#include "cli/cmd_align_internal.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "annot/annotation_store.h"
#include "checkpoint/checkpoint_cache.h"
#include "checkpoint/checkpoint_dispatcher.h"
#include "claude_agent/anthropic_client.h"
#include "checkpoint/checkpoint_types.h"
#include "classical/classical_pipeline.h"
#include "core/logging.h"
#include "core/thread_pool.h"
#include "core/transcript_mode.h"
#include "io/fasta_reader.h"
#include "io/fastq_reader.h"
#include "io/input_sniffer.h"
#include "output/bam_writer.h"
#include "output/parquet_writer.h"
#include "psv/psv_catalog.h"

namespace llmap::cli {

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

    // Resolve the input mode (Transcript-Mode architecture). `--mode` (if not
    // Auto) wins; otherwise the input sniffer decides from format / basename /
    // length stats. Emitted as a [mode-detect] log line.
    {
        const bool has_reads = !args.reads.empty();
        const bool has_assembly = !args.assembly.empty();
        const io::SniffResult sniff =
            io::ResolveMode(args.reads, args.mode, has_reads, has_assembly);
        args.resolved_mode = sniff.mode;
        char msg[512];
        std::snprintf(msg, sizeof(msg),
            "[mode-detect] detected=%s format=%s reason=\"%s\"",
            core::TranscriptModeName(sniff.mode),
            io::FileFormatName(sniff.format), sniff.reason.c_str());
        LLMAP_LOG_INFO(msg);
        if (args.verbose) {
            std::fprintf(stderr, "%s\n", msg);
        }
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

    // Optional region annotation
    std::unique_ptr<annot::AnnotationStore> region_store;
    if (!args.region_annot.empty()) {
        if (args.verbose) {
            std::fprintf(stderr, "Loading region annotation: %s\n",
                         args.region_annot.c_str());
        }
        region_store = annot::AnnotationStore::Load(
            args.region_annot, ref_names);
        if (!region_store) {
            std::fprintf(stderr,
                "Warning: failed to load region annotation %s, "
                "continuing without\n",
                args.region_annot.c_str());
        } else if (args.verbose) {
            std::fprintf(stderr, "  %zu annotation intervals loaded\n",
                         region_store->Size());
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
    if (region_store) {
        pipe_cfg.chain_config.annot = &region_store->Index();
    }
    // Transcript-Mode: a spliced read is reported by the chainer as several
    // per-exon chains. Report them all so the spliced stage can join them
    // (the joiner re-collapses them into one spliced alignment).
    if (args.resolved_mode == core::TranscriptMode::Transcript) {
        pipe_cfg.report_secondary = true;
        // Lossless / no-drop: short terminal or novel exons (alternative 3'
        // ends, NMD-escape last exons, short UTR exons) must survive the
        // length filter so the spliced stage can attach them — otherwise the
        // informative exon block is dropped before the joiner ever sees it and
        // the read collapses onto the wrong (annotated) isoform. Micro-exons
        // below the seed footprint are a separate (seeding) concern.
        pipe_cfg.min_aligned_bases =
            std::min(pipe_cfg.min_aligned_bases, 15);
        // R-B: break the chain at every intron-signature gap (large
        // reference-only gap), not just the ones exceeding max_gap_diff. This
        // turns each exon into a clean single-exon sub-chain — the precondition
        // the spliced stage needs to re-join them with N-ops. Without it, small
        // introds (<= max_gap_diff) link into one low-identity composite chain
        // that the identity filter then drops, leaving only one exon aligned
        // and the rest soft-clipped. Matches JoinerConfig::min_intron_bp (50).
        pipe_cfg.chain_config.intron_break_min = 50;
        // R-B (cont.): cap chain-end WFA extension so a sub-chain ending at an
        // intron break is NOT extended across the intron into the neighbouring
        // exon (that produced the `7=364M225I` identity-0.62 artefact the filter
        // then dropped). Beyond this flank the read is soft-clipped and the next
        // exon is handled as its own sub-chain. Kept above the seed window (k+w)
        // so the genuine exon-boundary flank is still recovered.
        pipe_cfg.extension_max_span = 50;
        // R-B (cont.): the per-exon sub-chains are COMPLEMENTARY (each covers a
        // different exon of the same read), not competing alignments of one
        // region. The default min_score_fraction (keep chains >= 0.5x best)
        // therefore wrongly drops the SHORTER exons — e.g. a 120 bp exon next
        // to a 240 bp one (0.5x240=120) is borderline-culled, so it never
        // reaches the joiner and the transcript loses that exon + its junction.
        // Drop the relative cull in transcript mode; the absolute floor
        // (min_chain_score) still removes spurious tiny chains.
        // Keep short exons (min_score_fraction culls relative to the best/longest
        // exon's chain) but not at 0.0 — on a real chromosome a 4 kb read throws
        // hundreds of scattered low-score minimizer-hit chains; 0.0 keeps them
        // all and the per-read WFA extension load explodes (300 reads vs chr14
        // did not finish in 120 s). A small floor drops the scattered junk while
        // the genuine exon chains (the top-scoring, colinear ones) survive.
        pipe_cfg.chain_config.min_score_fraction = 0.05f;
        // A transcript has many exons (sometimes dozens), so the default cap of
        // 5 is too low — but 256 over-extends on real references. 64 bounds the
        // per-read WFA work while still covering every exon of all but the most
        // extreme transcripts (the chains are extended best-score-first).
        pipe_cfg.max_chains_to_extend = 64;
        pipe_cfg.max_alignments = 64;
        // On a real chromosome a 4 kb read throws hundreds of scattered spurious
        // chains (repeat/paralog hits); the real exons cluster in one gene. Pull
        // that dominant colinear locus ahead of the noise so the 64-chain budget
        // hits the exons (else only 1/11 reads spliced vs minimap2's 18/24). The
        // span bounds the intra-locus intron gap — generous enough for the
        // largest human introns (~1 Mb) while separating far-scattered noise.
        pipe_cfg.transcript_locus_span = 1'000'000;
    }

    classical::ClassicalPipeline pipeline(pipe_cfg);
    pipeline.SetIndex(std::move(index));
    pipeline.SetReferenceSequences(ref_seqs);

    // Layer 3 active LLM consult dispatcher. Owned here; passed by pointer
    // into the chain-DP / extension layer (TODO: actual chain-DP-side
    // checkpoint firing -- see src/classical/chain_dp.cpp).
    checkpoint::LlmMode llm_mode = checkpoint::LlmMode::Off;
    if (args.llm_mode == "auto")     llm_mode = checkpoint::LlmMode::Auto;
    else if (args.llm_mode == "required") llm_mode = checkpoint::LlmMode::Required;

    std::unique_ptr<checkpoint::CheckpointCache> checkpoint_cache;
    std::unique_ptr<checkpoint::CheckpointDispatcher> checkpoint_dispatcher;
    std::unique_ptr<claude_agent::AnthropicClient> checkpoint_client;
    if (llm_mode != checkpoint::LlmMode::Off) {
        checkpoint_cache = std::make_unique<checkpoint::CheckpointCache>(
            args.llm_cache_dir.empty()
                ? std::filesystem::path{}
                : std::filesystem::path(args.llm_cache_dir));

        // Direct AnthropicClient for synchronous per-checkpoint calls.
        // Picks up the key from --llm-api-key, --psv-llm-api-key flow, or
        // env ANTHROPIC_API_KEY. If no key, the dispatcher falls back per
        // the LlmMode policy.
        std::string api_key = GetApiKey(args);
        if (!api_key.empty()) {
            claude_agent::AgentConfig acfg;
            acfg.api_key = api_key;
            acfg.model = "claude-sonnet-4-6";
            checkpoint_client =
                std::make_unique<claude_agent::AnthropicClient>(acfg);
            checkpoint_dispatcher =
                std::make_unique<checkpoint::CheckpointDispatcher>(
                    llm_mode, checkpoint_cache.get(),
                    checkpoint_client.get(),
                    checkpoint::CheckpointDispatcher::DirectClientTag{});
        } else {
            checkpoint_dispatcher =
                std::make_unique<checkpoint::CheckpointDispatcher>(
                    llm_mode, checkpoint_cache.get(),
                    static_cast<claude_agent::PipelineAgent*>(nullptr));
            if (llm_mode == checkpoint::LlmMode::Required) {
                std::fprintf(stderr,
                    "ERROR: --llm-mode=required but no API key (set "
                    "ANTHROPIC_API_KEY or use --llm-api-key)\n");
                return 1;
            }
        }
        if (args.verbose) {
            std::fprintf(stderr,
                "Layer 3 checkpoint dispatcher: mode=%s, cache=%s, "
                "api_key=%s\n",
                checkpoint::LlmModeName(llm_mode),
                checkpoint_cache->RootDir().c_str(),
                api_key.empty() ? "absent" : "present");
        }
    }

    if (args.verbose) {
        std::fprintf(stderr, "Loading reads: %s\n", args.reads.c_str());
    }

    auto read_reader = io::FastqReader::Open(args.reads);
    if (!read_reader) {
        std::fprintf(stderr, "Error: failed to open reads file: %s\n", args.reads.c_str());
        return 1;
    }

    // Open output writers
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
    // R-D: in Transcript-Mode unjoined exon fragments are emitted as
    // SUPPLEMENTARY (0x800) lines + SA links rather than an XA tag, so a
    // split spliced read keeps every exon block as a first-class alignment.
    bam_cfg.emit_split_as_supplementary =
        args.resolved_mode == core::TranscriptMode::Transcript;

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

    std::unique_ptr<core::ThreadPool> thread_pool;
    if (args.threads > 1) {
        thread_pool = std::make_unique<core::ThreadPool>(
            static_cast<size_t>(args.threads));
    }

    auto align_start = std::chrono::steady_clock::now();

    BatchAlignResult result = RunAlignmentBatches(
        pipeline, *read_reader, *bam_writer,
        parquet_writer.get(), thread_pool.get(),
        psv_catalog, args, ref_names, ref_seqs);

    if (result.error) {
        return 1;
    }

    auto align_end = std::chrono::steady_clock::now();
    float align_time_ms = std::chrono::duration<float, std::milli>(
        align_end - align_start).count();

    bam_writer->Close();
    if (parquet_writer) {
        parquet_writer->Close();
    }

    FinalizeAlignStats(result);

    auto end_time = std::chrono::steady_clock::now();
    float total_time_ms = std::chrono::duration<float, std::milli>(
        end_time - start_time).count();

    float mapping_rate = result.total_reads > 0
        ? static_cast<float>(result.n_mapped) / static_cast<float>(result.total_reads)
        : 0.0f;

    PrintAlignmentSummary(args, result.agg_stats, result.total_reads,
                          result.n_mapped, result.n_unmapped,
                          align_time_ms, total_time_ms);

    if (args.enable_llm && args.classical_only && args.verbose) {
        std::fprintf(stderr, "Note: --llm ignored due to --classical-only mode\n");
    }

    if (ShouldRunLlmDiagnostics(args, mapping_rate)) {
        if (args.verbose) {
            std::fprintf(stderr, "\nMapping rate %.1f%% below threshold %.1f%%, "
                         "running LLM diagnostics...\n",
                         100.0f * mapping_rate, 100.0f * args.llm_threshold);
        }

        auto agent = CreatePipelineAgent(args, args.verbose);
        if (agent) {
            std::vector<classical::ReadAlignmentResult> empty_results;
            RunLlmDiagnostics(*agent, args, empty_results,
                              result.n_mapped, result.n_unmapped,
                              result.agg_stats.avg_identity);
        }
    } else if (IsLlmEnabledButSkipped(args, mapping_rate) && args.verbose) {
        std::fprintf(stderr, "Mapping rate %.1f%% meets threshold, "
                     "skipping LLM diagnostics\n", 100.0f * mapping_rate);
    }

    return 0;
}

}  // namespace llmap::cli
