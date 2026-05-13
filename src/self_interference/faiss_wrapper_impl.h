#pragma once

#include "self_interference/faiss_wrapper.h"

#ifdef LLMAP_HAS_FAISS
#include <faiss/Index.h>
#ifdef LLMAP_HAS_FAISS_GPU
#include <faiss/gpu/StandardGpuResources.h>
#endif
#endif

namespace llmap::self_interference {

class FaissIndexImpl {
public:
#ifdef LLMAP_HAS_FAISS
    std::unique_ptr<faiss::Index> index;
#ifdef LLMAP_HAS_FAISS_GPU
    std::unique_ptr<faiss::gpu::StandardGpuResources> gpu_resources;
#endif
#endif
    FaissIndexStats stats;
    std::string last_error;
};

}  // namespace llmap::self_interference
