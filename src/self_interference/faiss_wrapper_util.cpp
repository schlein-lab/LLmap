#include "self_interference/faiss_wrapper_impl.h"

#include <cmath>

#ifdef LLMAP_HAS_FAISS
#include <faiss/index_io.h>
#endif

namespace llmap::self_interference {

bool FaissIndex::IsReady() const {
#ifdef LLMAP_HAS_FAISS
    return impl_->index != nullptr && impl_->stats.is_trained;
#else
    return false;
#endif
}

bool FaissIndex::IsTrained() const {
    return impl_->stats.is_trained;
}

size_t FaissIndex::NumVectors() const {
    return impl_->stats.num_vectors;
}

FaissIndexStats FaissIndex::GetStats() const {
    return impl_->stats;
}

bool FaissIndex::IsGpuEnabled() const {
    return active_provider_ == FaissProvider::GPU;
}

std::string FaissIndex::LastError() const {
    return impl_->last_error;
}

bool FaissIndex::Save(const std::string& path) const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) {
        return false;
    }
    try {
        faiss::write_index(impl_->index.get(), path.c_str());
        return true;
    } catch (const std::exception&) {
        return false;
    }
#else
    return false;
#endif
}

std::unique_ptr<FaissIndex> FaissIndex::Load(
    const std::string& path,
    const FaissIndexConfig& config)
{
#ifdef LLMAP_HAS_FAISS
    try {
        auto wrapper = std::unique_ptr<FaissIndex>(new FaissIndex(config));
        wrapper->impl_->index.reset(faiss::read_index(path.c_str()));

        if (!wrapper->impl_->index) {
            return nullptr;
        }

        wrapper->impl_->stats.num_vectors = wrapper->impl_->index->ntotal;
        wrapper->impl_->stats.embedding_dim = config.embedding_dim;
        wrapper->impl_->stats.is_trained = wrapper->impl_->index->is_trained;

        return wrapper;
    } catch (const std::exception&) {
        return nullptr;
    }
#else
    return nullptr;
#endif
}

void FaissIndex::NormalizeVectors(std::span<float> vectors, size_t num_vectors) const {
    const size_t dim = config_.embedding_dim;
    for (size_t i = 0; i < num_vectors; ++i) {
        float* vec = vectors.data() + i * dim;
        float norm = 0.0f;
        for (size_t j = 0; j < dim; ++j) {
            norm += vec[j] * vec[j];
        }
        norm = std::sqrt(norm);
        if (norm > 1e-8f) {
            const float inv_norm = 1.0f / norm;
            for (size_t j = 0; j < dim; ++j) {
                vec[j] *= inv_norm;
            }
        }
    }
}

std::vector<float> FaissIndex::NormalizeVectorsCopy(
    std::span<const float> vectors,
    size_t num_vectors) const
{
    std::vector<float> normalized(vectors.begin(),
        vectors.begin() + num_vectors * config_.embedding_dim);
    NormalizeVectors(normalized, num_vectors);
    return normalized;
}

bool IsFaissAvailable() {
#ifdef LLMAP_HAS_FAISS
    return true;
#else
    return false;
#endif
}

bool IsFaissGpuAvailable() {
#ifdef LLMAP_HAS_FAISS_GPU
    return true;
#else
    return false;
#endif
}

std::string GetFaissVersion() {
#ifdef LLMAP_HAS_FAISS
    return "FAISS (compiled)";
#else
    return "FAISS (not available)";
#endif
}

float ComputeRecallAtK(
    std::span<const int64_t> approx_indices,
    std::span<const int64_t> exact_indices,
    size_t num_queries,
    size_t k)
{
    if (num_queries == 0 || k == 0) {
        return 0.0f;
    }

    if (approx_indices.size() < num_queries * k ||
        exact_indices.size() < num_queries * k) {
        return 0.0f;
    }

    size_t total_found = 0;
    for (size_t q = 0; q < num_queries; ++q) {
        const auto* approx = approx_indices.data() + q * k;
        const auto* exact = exact_indices.data() + q * k;

        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                if (approx[i] == exact[j] && approx[i] >= 0) {
                    ++total_found;
                    break;
                }
            }
        }
    }

    return static_cast<float>(total_found) / static_cast<float>(num_queries * k);
}

}  // namespace llmap::self_interference
