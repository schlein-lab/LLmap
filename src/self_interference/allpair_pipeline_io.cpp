// LLmap — AllpairPipeline I/O stages.
//
// Contains: LoadReads, ComputeEmbeddings, WriteOutput

#include "self_interference/allpair_pipeline.h"
#include "self_interference/allpair_pipeline_internal.h"

#include <cmath>
#include <cstdio>
#include <span>

#include "ai/foundation_embedder.h"
#include "io/fastq_reader.h"
#include "output/cluster_writer.h"

namespace llmap::self_interference {

bool AllpairPipeline::LoadReads() {
    io::FastqReaderConfig reader_config;
    reader_config.min_length = config_.min_read_length;
    reader_config.max_length = config_.max_read_length;
    reader_config.max_records = config_.max_reads;

    auto reader = io::FastqReader::Open(config_.reads_a, reader_config);
    if (!reader) {
        state_->last_error = "Failed to open FASTQ file: " + config_.reads_a.string();
        return false;
    }

    state_->reads = reader->ReadAll();

    if (state_->reads.empty()) {
        state_->last_error = "No reads loaded from FASTQ file";
        return false;
    }

    auto stats = reader->GetStats();
    if (config_.verbose) {
        std::fprintf(stderr, "Loaded %zu reads, avg length %.1f bp\n",
                     stats.total_records, stats.avg_length);
    }

    return true;
}

bool AllpairPipeline::ComputeEmbeddings() {
    // If no model path provided, use placeholder embeddings for testing
    if (config_.model_path.empty() || !std::filesystem::exists(config_.model_path)) {
        // Generate deterministic pseudo-embeddings based on sequence content
        const size_t num_reads = state_->reads.size();
        const size_t dim = config_.embedding_dim;
        state_->embeddings.resize(num_reads * dim);

        for (size_t i = 0; i < num_reads; ++i) {
            const auto& seq = state_->reads[i].sequence;

            // Simple hash-based embedding for testing
            std::hash<std::string> hasher;
            size_t h = hasher(seq);

            for (size_t d = 0; d < dim; ++d) {
                h = h * 2654435761u + d;
                state_->embeddings[i * dim + d] =
                    static_cast<float>(h % 1000) / 500.0f - 1.0f;
            }

            // Normalize
            float norm = 0.0f;
            for (size_t d = 0; d < dim; ++d) {
                norm += state_->embeddings[i * dim + d] * state_->embeddings[i * dim + d];
            }
            norm = std::sqrt(norm);
            if (norm > 0) {
                for (size_t d = 0; d < dim; ++d) {
                    state_->embeddings[i * dim + d] /= norm;
                }
            }
        }

        return true;
    }

    // Use real embedder
    ai::EmbedderConfig embedder_config;
    embedder_config.model_path = config_.model_path;
    embedder_config.embedding_dim = config_.embedding_dim;
    embedder_config.batch_size = config_.embedding_batch_size;

    if (config_.use_gpu_embedder) {
        embedder_config.provider = ai::ExecutionProvider::CUDA;
        embedder_config.device_id = config_.gpu_device_id;
    }

    auto embedder = ai::FoundationEmbedder::Create(embedder_config);
    if (!embedder) {
        state_->last_error = "Failed to create embedder";
        return false;
    }

    const size_t num_reads = state_->reads.size();
    const size_t dim = embedder->EmbeddingDim();
    state_->embeddings.resize(num_reads * dim);

    // Embed in batches
    std::vector<std::string_view> batch_seqs;
    batch_seqs.reserve(config_.embedding_batch_size);

    for (size_t i = 0; i < num_reads; i += config_.embedding_batch_size) {
        batch_seqs.clear();
        const size_t batch_end = std::min(i + config_.embedding_batch_size, num_reads);

        for (size_t j = i; j < batch_end; ++j) {
            batch_seqs.push_back(state_->reads[j].sequence);
        }

        std::span<float> output_span(
            state_->embeddings.data() + i * dim,
            (batch_end - i) * dim);

        if (!embedder->EmbedBatchInto(batch_seqs, output_span)) {
            state_->last_error = "Embedding failed";
            return false;
        }
    }

    return true;
}

bool AllpairPipeline::WriteOutput() {
    if (!state_->swc_result) {
        state_->last_error = "No results to write";
        return false;
    }

    output::ClusterWriterConfig writer_config;
    writer_config.format = output::ClusterOutputFormat::TSV;
    writer_config.include_header = true;

    auto writer = output::ClusterWriter::Create(config_.output, writer_config);
    if (!writer) {
        state_->last_error = "Failed to create output writer";
        return false;
    }

    const size_t num_reads = state_->reads.size();
    std::vector<size_t> cluster_sizes;
    if (state_->swc_result) {
        cluster_sizes = state_->swc_result->GetClusterSizes();
    }

    for (size_t i = 0; i < num_reads; ++i) {
        output::ClusterAssignment assignment;
        assignment.read_id = state_->reads[i].id;
        assignment.read_length = state_->reads[i].sequence.size();

        // Find this read's assignment
        bool found = false;
        for (const auto& a : state_->swc_result->assignments) {
            if (a.read_idx == i) {
                assignment.cluster_id = a.cluster_id;
                assignment.confidence = a.confidence;
                assignment.cluster_size = (a.cluster_id < cluster_sizes.size())
                    ? cluster_sizes[a.cluster_id] : 0;
                found = true;
                break;
            }
        }

        if (!found) {
            // Fallback: use Leiden assignment
            if (state_->clustering && i < state_->clustering->labels.size()) {
                assignment.cluster_id = state_->clustering->labels[i];
                assignment.confidence = 0.5f;
            } else {
                assignment.cluster_id = 0;
                assignment.confidence = 0.0f;
            }
            assignment.cluster_size = 0;
        }

        // Check if representative
        assignment.is_representative = false;
        if (state_->rep_result) {
            assignment.is_representative = state_->rep_result->IsRepresentative(
                static_cast<uint32_t>(i));
        }

        if (!writer->Write(assignment)) {
            state_->last_error = "Failed to write assignment: " + writer->LastError();
            return false;
        }
    }

    writer->Close();
    return true;
}

}  // namespace llmap::self_interference
