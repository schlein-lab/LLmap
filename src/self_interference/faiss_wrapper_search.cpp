#include "self_interference/faiss_wrapper_impl.h"

#include <chrono>

#ifdef LLMAP_HAS_FAISS
#include <faiss/IndexIVF.h>
#ifdef LLMAP_HAS_FAISS_GPU
#include <faiss/gpu/GpuIndexIVF.h>
#endif
#endif

namespace llmap::self_interference {

std::optional<KnnResult> FaissIndex::Search(std::span<const float> query, size_t k) const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        return std::nullopt;
    }

    if (query.size() < config_.embedding_dim) {
        return std::nullopt;
    }

    try {
        KnnResult result;
        result.indices.resize(k);
        result.distances.resize(k);

        const float* data_ptr = query.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(query, 1);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            1,
            data_ptr,
            static_cast<faiss::idx_t>(k),
            result.distances.data(),
            result.indices.data());

        return result;

    } catch (const std::exception&) {
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

std::optional<BatchKnnResult> FaissIndex::SearchBatch(
    std::span<const float> queries,
    size_t num_queries,
    size_t k) const
{
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        impl_->last_error = "Index empty or not initialized";
        return std::nullopt;
    }

    if (queries.size() < num_queries * config_.embedding_dim) {
        impl_->last_error = "Query buffer too small";
        return std::nullopt;
    }

    try {
        auto start = std::chrono::high_resolution_clock::now();

        BatchKnnResult result;
        result.num_queries = num_queries;
        result.k = k;
        result.indices.resize(num_queries * k);
        result.distances.resize(num_queries * k);

        const float* data_ptr = queries.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(queries, num_queries);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            static_cast<faiss::idx_t>(num_queries),
            data_ptr,
            static_cast<faiss::idx_t>(k),
            result.distances.data(),
            result.indices.data());

        auto end = std::chrono::high_resolution_clock::now();
        result.search_time_ms = std::chrono::duration<float, std::milli>(end - start).count();

        return result;

    } catch (const std::exception& e) {
        impl_->last_error = std::string("SearchBatch failed: ") + e.what();
        return std::nullopt;
    }
#else
    return std::nullopt;
#endif
}

bool FaissIndex::SearchBatchInto(
    std::span<const float> queries,
    size_t num_queries,
    size_t k,
    std::span<int64_t> indices,
    std::span<float> distances) const
{
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index || impl_->stats.num_vectors == 0) {
        return false;
    }

    if (indices.size() < num_queries * k || distances.size() < num_queries * k) {
        return false;
    }

    try {
        const float* data_ptr = queries.data();
        std::vector<float> normalized;
        if (config_.normalize_vectors) {
            normalized = NormalizeVectorsCopy(queries, num_queries);
            data_ptr = normalized.data();
        }

        impl_->index->search(
            static_cast<faiss::idx_t>(num_queries),
            data_ptr,
            static_cast<faiss::idx_t>(k),
            distances.data(),
            indices.data());

        return true;

    } catch (const std::exception&) {
        return false;
    }
#else
    return false;
#endif
}

void FaissIndex::SetNprobe(size_t nprobe) {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) return;

    if (auto* ivf = dynamic_cast<faiss::IndexIVF*>(impl_->index.get())) {
        ivf->nprobe = nprobe;
    }
#ifdef LLMAP_HAS_FAISS_GPU
    if (auto* gpu_ivf = dynamic_cast<faiss::gpu::GpuIndexIVF*>(impl_->index.get())) {
        gpu_ivf->setNumProbes(static_cast<int>(nprobe));
    }
#endif
#endif
}

size_t FaissIndex::GetNprobe() const {
#ifdef LLMAP_HAS_FAISS
    if (!impl_->index) return 0;

    if (auto* ivf = dynamic_cast<faiss::IndexIVF*>(impl_->index.get())) {
        return ivf->nprobe;
    }
#ifdef LLMAP_HAS_FAISS_GPU
    if (auto* gpu_ivf = dynamic_cast<faiss::gpu::GpuIndexIVF*>(impl_->index.get())) {
        return gpu_ivf->getNumProbes();
    }
#endif
#endif
    return config_.nprobe;
}

}  // namespace llmap::self_interference
