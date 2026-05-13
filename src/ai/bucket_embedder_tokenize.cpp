#include "ai/bucket_embedder.h"

#include <algorithm>
#include <cmath>

namespace llmap::ai {

namespace {

// DNA nucleotide to token ID mapping
constexpr int64_t kPadToken = 0;
constexpr int64_t kClsToken = 2;
constexpr int64_t kSepToken = 3;
constexpr int64_t kAToken = 4;
constexpr int64_t kCToken = 5;
constexpr int64_t kGToken = 6;
constexpr int64_t kTToken = 7;
constexpr int64_t kNToken = 8;
constexpr int64_t kUnkToken = 1;

int64_t NucleotideToToken(char c) {
    switch (c) {
        case 'A': case 'a': return kAToken;
        case 'C': case 'c': return kCToken;
        case 'G': case 'g': return kGToken;
        case 'T': case 't': return kTToken;
        case 'N': case 'n': return kNToken;
        default: return kUnkToken;
    }
}

}  // namespace

std::vector<int64_t> BucketEmbedder::Tokenize(std::string_view sequence) const {
    std::vector<int64_t> tokens;
    size_t max_len = config_.max_sequence_length;
    tokens.reserve(max_len + 2);

    tokens.push_back(kClsToken);

    size_t seq_len = std::min(sequence.length(), max_len);
    for (size_t i = 0; i < seq_len; ++i) {
        tokens.push_back(NucleotideToToken(sequence[i]));
    }

    tokens.push_back(kSepToken);

    while (tokens.size() < max_len + 2) {
        tokens.push_back(kPadToken);
    }

    return tokens;
}

void BucketEmbedder::TokenizeBatch(
    std::span<const std::string_view> sequences,
    std::vector<int64_t>& token_ids,
    std::vector<int64_t>& attention_mask) const {

    size_t batch_size = sequences.size();
    size_t seq_len = config_.max_sequence_length + 2;

    token_ids.resize(batch_size * seq_len);
    attention_mask.resize(batch_size * seq_len);

    for (size_t b = 0; b < batch_size; ++b) {
        size_t offset = b * seq_len;
        const auto& seq = sequences[b];

        token_ids[offset] = kClsToken;
        attention_mask[offset] = 1;

        size_t actual_len = std::min(seq.length(), config_.max_sequence_length);
        for (size_t i = 0; i < actual_len; ++i) {
            token_ids[offset + 1 + i] = NucleotideToToken(seq[i]);
            attention_mask[offset + 1 + i] = 1;
        }

        token_ids[offset + 1 + actual_len] = kSepToken;
        attention_mask[offset + 1 + actual_len] = 1;

        for (size_t i = actual_len + 2; i < seq_len; ++i) {
            token_ids[offset + i] = kPadToken;
            attention_mask[offset + i] = 0;
        }
    }
}

void BucketEmbedder::NormalizeEmbedding(std::span<float> embedding) const {
    float norm_sq = 0.0f;
    for (float v : embedding) {
        norm_sq += v * v;
    }
    if (norm_sq > 1e-12f) {
        float inv_norm = 1.0f / std::sqrt(norm_sq);
        for (float& v : embedding) {
            v *= inv_norm;
        }
    }
}

std::vector<float> BucketEmbedder::PoolChunkEmbeddings(
    const std::vector<std::vector<float>>& chunk_embeddings) const {

    if (chunk_embeddings.empty()) {
        return std::vector<float>(config_.embedding_dim, 0.0f);
    }

    if (chunk_embeddings.size() == 1) {
        return chunk_embeddings[0];
    }

    size_t dim = chunk_embeddings[0].size();
    std::vector<float> pooled(dim, 0.0f);

    for (const auto& chunk : chunk_embeddings) {
        for (size_t i = 0; i < dim; ++i) {
            pooled[i] += chunk[i];
        }
    }

    float inv_n = 1.0f / static_cast<float>(chunk_embeddings.size());
    for (float& v : pooled) {
        v *= inv_n;
    }

    if (config_.normalize_output) {
        NormalizeEmbedding(pooled);
    }

    return pooled;
}

}  // namespace llmap::ai
