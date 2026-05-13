add_library(llmap_core STATIC
    core/alignment_record.cpp
    core/bucket_pyramid.cpp
    core/wave_state.cpp
    synthetic/igh_locus_generator.cpp
)

target_include_directories(llmap_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_core PUBLIC cxx_std_23)

# IO module — file readers (FASTQ, BAM stubs)
add_library(llmap_io STATIC
    io/fastq_reader.cpp
)

target_include_directories(llmap_io PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_io PUBLIC cxx_std_23)
target_link_libraries(llmap_io PUBLIC llmap_core)

# Output module — writers (cluster assignments, BAM/Parquet stubs)
add_library(llmap_output STATIC
    output/cluster_writer.cpp
)

target_include_directories(llmap_output PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_output PUBLIC cxx_std_23)
target_link_libraries(llmap_output PUBLIC llmap_core)

# AI module — foundation model embedders
add_library(llmap_ai STATIC
    ai/foundation_embedder.cpp
    ai/bucket_embedder.cpp
)

target_include_directories(llmap_ai PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_ai PUBLIC cxx_std_23)
target_link_libraries(llmap_ai PUBLIC llmap_core)

# Link ONNX Runtime if available
if(LLMAP_ONNXRUNTIME_AVAILABLE)
    target_compile_definitions(llmap_ai PRIVATE LLMAP_HAS_ONNXRUNTIME)
    if(TARGET onnxruntime::onnxruntime)
        target_link_libraries(llmap_ai PUBLIC onnxruntime::onnxruntime)
    elseif(ONNXRUNTIME_LIBRARY)
        target_link_libraries(llmap_ai PUBLIC ${ONNXRUNTIME_LIBRARY})
        target_include_directories(llmap_ai PUBLIC ${ONNXRUNTIME_INCLUDE_DIR})
    endif()
endif()

# Self-interference module — Stage 1 (read-vs-read ANN clustering)
add_library(llmap_self_interference STATIC
    self_interference/faiss_wrapper_index.cpp
    self_interference/faiss_wrapper_search.cpp
    self_interference/faiss_wrapper_util.cpp
    self_interference/similarity_graph.cpp
    self_interference/leiden_clustering_core.cpp
    self_interference/leiden_clustering_phases.cpp
    self_interference/leiden_clustering_result.cpp
    self_interference/self_wavecollapse.cpp
    self_interference/cluster_rep_selector.cpp
    self_interference/cluster_rep_methods.cpp
    self_interference/cluster_rep_result.cpp
    self_interference/allpair_pipeline.cpp
    self_interference/allpair_pipeline_io.cpp
    self_interference/allpair_pipeline_stages.cpp
)

target_include_directories(llmap_self_interference PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_self_interference PUBLIC cxx_std_23)
target_link_libraries(llmap_self_interference PUBLIC llmap_core llmap_ai llmap_io llmap_output)

# Link FAISS if available
if(LLMAP_FAISS_AVAILABLE)
    target_compile_definitions(llmap_self_interference PRIVATE LLMAP_HAS_FAISS)
    if(TARGET faiss::faiss)
        target_link_libraries(llmap_self_interference PUBLIC faiss::faiss)
    elseif(FAISS_LIBRARY)
        target_link_libraries(llmap_self_interference PUBLIC ${FAISS_LIBRARY})
        target_include_directories(llmap_self_interference PUBLIC ${FAISS_INCLUDE_DIR})
    endif()
    if(LLMAP_FAISS_GPU_AVAILABLE)
        target_compile_definitions(llmap_self_interference PRIVATE LLMAP_HAS_FAISS_GPU)
        target_link_libraries(llmap_self_interference PUBLIC ${FAISS_GPU_LIBRARY})
    endif()
endif()

# Reference collapse module — Stage 2 (reference-based WaveCollapse)
add_library(llmap_reference_collapse STATIC
    reference_collapse/reference_index.cpp
    reference_collapse/reference_index_builder.cpp
    reference_collapse/reference_index_io.cpp
    reference_collapse/em_iterator.cpp
    reference_collapse/collapse_check.cpp
    reference_collapse/refinement.cpp
    reference_collapse/member_propagation.cpp
    reference_collapse/stage2_pipeline.cpp
)

target_include_directories(llmap_reference_collapse PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_reference_collapse PUBLIC cxx_std_23)
target_link_libraries(llmap_reference_collapse PUBLIC llmap_core llmap_ai llmap_self_interference)

# CLI — full CLI surface
add_executable(llmap cli/llmap_main.cpp)
target_link_libraries(llmap PRIVATE llmap_core llmap_ai llmap_io llmap_output llmap_self_interference llmap_reference_collapse)
