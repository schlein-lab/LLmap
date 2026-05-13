// LLmap — ClassicalPipeline: seed-chain-extend orchestration.

#include "classical/classical_pipeline.h"

#include <algorithm>
#include <chrono>
#include <sstream>

namespace llmap::classical {

namespace {

// Helper to measure time
class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}

    float ElapsedUs() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<float, std::micro>(now - start_).count();
    }

    float ElapsedMs() const {
        return ElapsedUs() / 1000.0f;
    }

private:
    std::chrono::steady_clock::time_point start_;
};

}  // namespace

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

std::optional<ClassicalAlignment> ClassicalPipeline::ExtendChain(
    std::string_view query_seq,
    const Chain& chain,
    const std::vector<Anchor>& anchors) const {

    if (chain.anchors.empty()) {
        return std::nullopt;
    }

    ClassicalAlignment aln;
    aln.ref_start = chain.ref_start;
    aln.ref_end = chain.ref_end;
    aln.query_start = chain.query_start;
    aln.query_end = chain.query_end;
    aln.score = chain.score;

    // Build a simple CIGAR from chain anchors
    // For a proper implementation, we would do full WFA2 extension
    // between anchors. For now, we approximate with M operations
    // connecting the anchors.

    uint32_t prev_query = chain.query_start;
    uint32_t prev_ref = chain.ref_start;
    size_t matches = 0;
    size_t insertions = 0;
    size_t deletions = 0;

    for (uint32_t anchor_idx : chain.anchors) {
        if (anchor_idx >= anchors.size()) continue;
        const auto& anchor = anchors[anchor_idx];

        // Gap between previous position and this anchor
        int32_t query_gap = static_cast<int32_t>(anchor.query_pos) - static_cast<int32_t>(prev_query);
        int32_t ref_gap = static_cast<int32_t>(anchor.ref_pos) - static_cast<int32_t>(prev_ref);

        if (query_gap > 0 && ref_gap > 0) {
            // Aligned segment
            uint32_t aligned = static_cast<uint32_t>(std::min(query_gap, ref_gap));
            if (aligned > 0) {
                aln.cigar.push_back({CigarOp::Match, aligned});
                matches += aligned;
            }

            // Handle gaps
            if (query_gap > ref_gap) {
                uint32_t ins = static_cast<uint32_t>(query_gap - ref_gap);
                aln.cigar.push_back({CigarOp::Insertion, ins});
                insertions += ins;
            } else if (ref_gap > query_gap) {
                uint32_t del = static_cast<uint32_t>(ref_gap - query_gap);
                aln.cigar.push_back({CigarOp::Deletion, del});
                deletions += del;
            }
        } else if (query_gap > 0) {
            aln.cigar.push_back({CigarOp::Insertion, static_cast<uint32_t>(query_gap)});
            insertions += query_gap;
        } else if (ref_gap > 0) {
            aln.cigar.push_back({CigarOp::Deletion, static_cast<uint32_t>(ref_gap)});
            deletions += ref_gap;
        }

        // k-mer match
        uint8_t k = config_.minimizer_config.k;
        aln.cigar.push_back({CigarOp::Match, k});
        matches += k;

        prev_query = anchor.query_pos + k;
        prev_ref = anchor.ref_pos + k;
    }

    // Handle trailing gap
    if (prev_query < chain.query_end) {
        uint32_t trail = chain.query_end - prev_query;
        aln.cigar.push_back({CigarOp::Match, trail});
        matches += trail;
    }

    // Merge adjacent CIGAR operations
    std::vector<CigarElement> merged;
    for (const auto& elem : aln.cigar) {
        if (elem.length == 0) continue;
        if (!merged.empty() && merged.back().op == elem.op) {
            merged.back().length += elem.length;
        } else {
            merged.push_back(elem);
        }
    }
    aln.cigar = std::move(merged);

    // Statistics
    size_t total_aligned = matches + insertions + deletions;
    aln.identity = total_aligned > 0
        ? static_cast<float>(matches) / static_cast<float>(total_aligned)
        : 0.0f;

    return aln;
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
    return pipeline.AlignReads(query_names, query_seqs);
}

}  // namespace llmap::classical
