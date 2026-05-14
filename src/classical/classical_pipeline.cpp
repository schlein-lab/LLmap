// LLmap — ClassicalPipeline: seed-chain-extend orchestration.

#include "classical/classical_pipeline.h"
#include "classical/classical_pipeline_internal.h"
#include "core/thread_pool.h"

#include <algorithm>
#include <atomic>
#include <sstream>

namespace llmap::classical {

std::string ClassicalAlignment::CigarString() const {
    std::ostringstream oss;
    for (const auto& elem : cigar) {
        oss << elem.ToString();
    }
    return oss.str();
}

uint32_t ClassicalAlignment::AlignedBases() const {
    uint32_t bases = 0;
    for (const auto& elem : cigar) {
        if (elem.op == CigarOp::Match || elem.op == CigarOp::Equal ||
            elem.op == CigarOp::Diff) {
            bases += elem.length;
        }
    }
    return bases;
}

const ClassicalAlignment* ReadAlignmentResult::PrimaryAlignment() const {
    for (const auto& aln : alignments) {
        if (aln.is_primary) {
            return &aln;
        }
    }
    return nullptr;
}

ClassicalPipeline::ClassicalPipeline() : ClassicalPipeline(ClassicalPipelineConfig{}) {}

ClassicalPipeline::ClassicalPipeline(const ClassicalPipelineConfig& config)
    : config_(config), aligner_(config.extension_config) {}

ClassicalPipeline::~ClassicalPipeline() = default;

ClassicalPipeline::ClassicalPipeline(ClassicalPipeline&&) noexcept = default;
ClassicalPipeline& ClassicalPipeline::operator=(ClassicalPipeline&&) noexcept = default;

void ClassicalPipeline::SetIndex(std::unique_ptr<MinimizerIndex> index) {
    owned_index_ = std::move(index);
    index_ = owned_index_.get();
}

void ClassicalPipeline::SetIndex(const MinimizerIndex* index) {
    owned_index_.reset();
    index_ = index;
}

void ClassicalPipeline::SetReferenceSequences(std::span<const std::string> ref_seqs) {
    owned_ref_seqs_.clear();
    ref_seqs_ = ref_seqs;
}

void ClassicalPipeline::SetReferenceSequences(std::vector<std::string> ref_seqs) {
    owned_ref_seqs_ = std::move(ref_seqs);
    ref_seqs_ = owned_ref_seqs_;
}

bool ClassicalPipeline::HasIndex() const {
    return index_ != nullptr && !index_->Empty();
}

ReadAlignmentResult ClassicalPipeline::AlignRead(
    std::string_view query_name,
    std::string_view query_seq) const {

    ReadAlignmentResult result;
    result.query_name = std::string(query_name);

    if (!HasIndex() || query_seq.empty()) {
        return result;
    }

    // Phase 1: Seeding
    Timer seed_timer;
    auto hits = index_->Query(query_seq);
    result.seeding_time_us = seed_timer.ElapsedUs();
    result.num_hits = hits.size();

    if (hits.empty()) {
        return result;
    }

    // Phase 2: Chaining
    Timer chain_timer;
    auto chain_result = ExtractChains(
        std::span<const MinimizerHit>(hits),
        static_cast<uint32_t>(query_seq.size()),
        config_.chain_config);
    result.chaining_time_us = chain_timer.ElapsedUs();
    result.num_chains = chain_result.chains.size();

    if (chain_result.chains.empty()) {
        return result;
    }

    // Convert hits to anchors for extension
    std::vector<Anchor> anchors;
    anchors.reserve(hits.size());
    for (const auto& hit : hits) {
        anchors.push_back({
            .ref_id = hit.ref_id,
            .ref_pos = hit.ref_pos,
            .query_pos = hit.query_pos,
            .same_strand = hit.same_strand
        });
    }

    // Phase 3: Extension (for top chains)
    Timer ext_timer;
    const auto& seqs = index_->GetSequences();
    uint32_t chains_to_try = std::min(
        config_.max_chains_to_extend,
        static_cast<uint32_t>(chain_result.chains.size()));

    bool found_primary = false;
    for (uint32_t i = 0; i < chains_to_try; ++i) {
        const auto& chain = chain_result.chains[i];

        auto alignment = ExtendChain(query_seq, chain, anchors);
        if (!alignment) continue;

        // Fill in metadata
        alignment->query_name = std::string(query_name);
        if (chain.ref_id < seqs.size()) {
            alignment->ref_name = seqs[chain.ref_id].name;
        }
        alignment->ref_id = chain.ref_id;
        alignment->chain_anchors = chain.NumAnchors();
        alignment->chain_score = chain.score;
        alignment->is_forward = chain.is_forward;
        alignment->is_primary = !found_primary;

        // Compute MAPQ (simplified: based on score gap to secondary)
        if (!found_primary) {
            found_primary = true;
            if (chain_result.chains.size() > 1 &&
                chain_result.chains[1].score > 0) {
                float score_ratio = static_cast<float>(chain_result.chains[1].score) /
                                    static_cast<float>(chain.score);
                alignment->mapq = static_cast<uint32_t>((1.0f - score_ratio) * 60);
            } else {
                alignment->mapq = 60;  // Unique mapping
            }
        } else {
            alignment->mapq = 0;  // Secondary
        }

        // Filter by identity
        if (alignment->identity >= config_.min_identity &&
            alignment->AlignedBases() >= static_cast<uint32_t>(config_.min_aligned_bases)) {
            result.alignments.push_back(std::move(*alignment));
            ++result.chains_extended;

            if (!config_.report_secondary && result.alignments.size() >= 1) {
                break;
            }
            if (result.alignments.size() >= config_.max_alignments) {
                break;
            }
        }
    }

    result.extension_time_us = ext_timer.ElapsedUs();
    return result;
}

std::vector<ReadAlignmentResult> ClassicalPipeline::AlignReads(
    std::span<const std::string> query_names,
    std::span<const std::string> query_seqs) const {

    Timer total_timer;
    stats_ = {};
    stats_.total_reads = query_names.size();

    std::vector<ReadAlignmentResult> results;
    results.reserve(query_names.size());

    for (size_t i = 0; i < query_names.size(); ++i) {
        auto result = AlignRead(query_names[i], query_seqs[i]);

        stats_.total_hits += result.num_hits;
        stats_.total_chains += result.num_chains;
        stats_.total_extensions += result.chains_extended;
        stats_.seeding_time_ms += result.seeding_time_us / 1000.0f;
        stats_.chaining_time_ms += result.chaining_time_us / 1000.0f;
        stats_.extension_time_ms += result.extension_time_us / 1000.0f;

        if (result.HasAlignment()) {
            ++stats_.reads_aligned;
            stats_.avg_identity += result.alignments[0].identity;
        } else {
            ++stats_.reads_unmapped;
        }

        results.push_back(std::move(result));
    }

    if (stats_.reads_aligned > 0) {
        stats_.avg_identity /= static_cast<float>(stats_.reads_aligned);
    }
    stats_.total_time_ms = total_timer.ElapsedMs();

    return results;
}

std::vector<ReadAlignmentResult> ClassicalPipeline::AlignReadsParallel(
    std::span<const std::string> query_names,
    std::span<const std::string> query_seqs,
    core::ThreadPool& pool) const {

    Timer total_timer;
    stats_ = {};
    stats_.total_reads = query_names.size();

    std::vector<ReadAlignmentResult> results(query_names.size());

    // Thread-safe counters for stats aggregation
    std::atomic<size_t> total_hits{0};
    std::atomic<size_t> total_chains{0};
    std::atomic<size_t> total_extensions{0};
    std::atomic<size_t> reads_aligned{0};
    std::atomic<size_t> reads_unmapped{0};
    std::atomic<uint64_t> seeding_time_us{0};
    std::atomic<uint64_t> chaining_time_us{0};
    std::atomic<uint64_t> extension_time_us{0};
    std::atomic<uint64_t> identity_sum_scaled{0};  // Scaled by 1e6 for precision

    core::ParallelFor(pool, query_names.size(), [&](size_t i) {
        auto result = AlignRead(query_names[i], query_seqs[i]);

        total_hits.fetch_add(result.num_hits, std::memory_order_relaxed);
        total_chains.fetch_add(result.num_chains, std::memory_order_relaxed);
        total_extensions.fetch_add(result.chains_extended, std::memory_order_relaxed);
        seeding_time_us.fetch_add(
            static_cast<uint64_t>(result.seeding_time_us), std::memory_order_relaxed);
        chaining_time_us.fetch_add(
            static_cast<uint64_t>(result.chaining_time_us), std::memory_order_relaxed);
        extension_time_us.fetch_add(
            static_cast<uint64_t>(result.extension_time_us), std::memory_order_relaxed);

        if (result.HasAlignment()) {
            reads_aligned.fetch_add(1, std::memory_order_relaxed);
            identity_sum_scaled.fetch_add(
                static_cast<uint64_t>(result.alignments[0].identity * 1e6),
                std::memory_order_relaxed);
        } else {
            reads_unmapped.fetch_add(1, std::memory_order_relaxed);
        }

        results[i] = std::move(result);
    });

    // Aggregate stats
    stats_.total_hits = total_hits.load();
    stats_.total_chains = total_chains.load();
    stats_.total_extensions = total_extensions.load();
    stats_.reads_aligned = reads_aligned.load();
    stats_.reads_unmapped = reads_unmapped.load();
    stats_.seeding_time_ms = static_cast<float>(seeding_time_us.load()) / 1000.0f;
    stats_.chaining_time_ms = static_cast<float>(chaining_time_us.load()) / 1000.0f;
    stats_.extension_time_ms = static_cast<float>(extension_time_us.load()) / 1000.0f;

    if (stats_.reads_aligned > 0) {
        stats_.avg_identity = static_cast<float>(identity_sum_scaled.load()) /
                             (static_cast<float>(stats_.reads_aligned) * 1e6f);
    }
    stats_.total_time_ms = total_timer.ElapsedMs();

    return results;
}

std::vector<ReadAlignmentResult> AlignWithClassicalPath(
    std::span<const std::string> ref_names,
    std::span<const std::string> ref_seqs,
    std::span<const std::string> query_names,
    std::span<const std::string> query_seqs,
    const ClassicalPipelineConfig& config) {

    // Build index
    MinimizerIndex::Builder builder(config.minimizer_config);
    for (size_t i = 0; i < ref_names.size(); ++i) {
        builder.AddSequence(ref_names[i], ref_seqs[i]);
    }
    auto index = builder.Build();

    // Create pipeline and align
    ClassicalPipeline pipeline(config);
    pipeline.SetIndex(std::move(index));
    pipeline.SetReferenceSequences(ref_seqs);  // Enable WFA2 extension
    return pipeline.AlignReads(query_names, query_seqs);
}

}  // namespace llmap::classical
