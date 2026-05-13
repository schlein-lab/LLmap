#include "ai/foundation_embedder.h"
#include "ai/foundation_embedder_impl.h"

#include <algorithm>

namespace llmap::ai {

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

bool IsValidNucleotide(char c) {
    switch (c) {
        case 'A': case 'a':
        case 'C': case 'c':
        case 'G': case 'g':
        case 'T': case 't':
        case 'N': case 'n':
            return true;
        default:
            return false;
    }
}

std::vector<int64_t> FoundationEmbedder::Tokenize(std::string_view sequence) const {
    std::vector<int64_t> tokens;
    tokens.reserve(config_.max_sequence_length + 2);  // +2 for CLS and SEP

    // Add CLS token
    tokens.push_back(kClsToken);

    // Tokenize sequence
    size_t seq_len = std::min(sequence.length(), config_.max_sequence_length);
    for (size_t i = 0; i < seq_len; ++i) {
        tokens.push_back(NucleotideToToken(sequence[i]));
    }

    // Add SEP token
    tokens.push_back(kSepToken);

    // Pad to max length
    while (tokens.size() < config_.max_sequence_length + 2) {
        tokens.push_back(kPadToken);
    }

    return tokens;
}

void FoundationEmbedder::TokenizeBatch(
    std::span<const std::string_view> sequences,
    std::vector<int64_t>& token_ids,
    std::vector<int64_t>& attention_mask) const {

    size_t batch_size = sequences.size();
    size_t seq_len = config_.max_sequence_length + 2;  // +2 for CLS and SEP

    token_ids.resize(batch_size * seq_len);
    attention_mask.resize(batch_size * seq_len);

    for (size_t b = 0; b < batch_size; ++b) {
        size_t offset = b * seq_len;
        const auto& seq = sequences[b];

        // CLS token
        token_ids[offset] = kClsToken;
        attention_mask[offset] = 1;

        // Sequence tokens
        size_t actual_len = std::min(seq.length(), config_.max_sequence_length);
        for (size_t i = 0; i < actual_len; ++i) {
            token_ids[offset + 1 + i] = NucleotideToToken(seq[i]);
            attention_mask[offset + 1 + i] = 1;
        }

        // SEP token
        token_ids[offset + 1 + actual_len] = kSepToken;
        attention_mask[offset + 1 + actual_len] = 1;

        // Padding
        for (size_t i = actual_len + 2; i < seq_len; ++i) {
            token_ids[offset + i] = kPadToken;
            attention_mask[offset + i] = 0;
        }
    }
}

}  // namespace llmap::ai
