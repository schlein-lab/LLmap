add_library(llmap_core STATIC
    core/alignment_record.cpp
    core/bucket_pyramid.cpp
    core/wave_state.cpp
    synthetic/igh_locus_generator_core.cpp
    synthetic/igh_locus_generator_presets.cpp
)

target_include_directories(llmap_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_core PUBLIC cxx_std_23)

# IO module — file readers (FASTA, FASTQ, BAM stubs)
add_library(llmap_io STATIC
    io/fasta_reader.cpp
    io/fastq_reader.cpp
)

target_include_directories(llmap_io PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_io PUBLIC cxx_std_23)
target_link_libraries(llmap_io PUBLIC llmap_core)

# Output module — writers (cluster assignments, BAM/Parquet)
add_library(llmap_output STATIC
    output/cluster_writer.cpp
    output/bam_writer.cpp
    output/bam_writer_cigar.cpp
    output/parquet_writer.cpp
    output/parquet_writer_convert.cpp
    output/parquet_writer_write.cpp
    output/parquet_reader.cpp
    output/parquet_reader_read.cpp
    output/parquet_reader_util.cpp
)

target_include_directories(llmap_output PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_output PUBLIC cxx_std_23)
target_link_libraries(llmap_output PUBLIC llmap_core)

# AI module — foundation model embedders
add_library(llmap_ai STATIC
    ai/foundation_embedder_core.cpp
    ai/foundation_embedder_tokenize.cpp
    ai/foundation_embedder_embed.cpp
    ai/bucket_embedder_core.cpp
    ai/bucket_embedder_tokenize.cpp
    ai/bucket_embedder_embed.cpp
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
    self_interference/similarity_graph_build.cpp
    self_interference/similarity_graph_query.cpp
    self_interference/similarity_graph_io.cpp
    self_interference/leiden_clustering_core.cpp
    self_interference/leiden_clustering_phases.cpp
    self_interference/leiden_clustering_result.cpp
    self_interference/self_wavecollapse_core.cpp
    self_interference/self_wavecollapse_em.cpp
    self_interference/self_wavecollapse_result.cpp
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
    reference_collapse/stage2_pipeline_core.cpp
    reference_collapse/stage2_pipeline_init.cpp
    reference_collapse/stage2_pipeline_em.cpp
)

target_include_directories(llmap_reference_collapse PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_reference_collapse PUBLIC cxx_std_23)
target_link_libraries(llmap_reference_collapse PUBLIC llmap_core llmap_ai llmap_self_interference)

# Classical module — minimizers, chains, WFA2 extension
add_library(llmap_classical STATIC
    classical/minimizer_index_core.cpp
    classical/minimizer_index_io.cpp
    classical/minimizer_index_extract.cpp
    classical/chain_score.cpp
    classical/chain_dp.cpp
    classical/wfa2_aligner_core.cpp
    classical/wfa2_aligner_align.cpp
    classical/classical_pipeline.cpp
)

target_include_directories(llmap_classical PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_classical PUBLIC cxx_std_23)
target_link_libraries(llmap_classical PUBLIC llmap_core)

# Link WFA2-lib if available
if(LLMAP_WFA2_AVAILABLE)
    target_compile_definitions(llmap_classical PRIVATE LLMAP_HAS_WFA2)
    target_link_libraries(llmap_classical PUBLIC ${WFA2_LIBRARY})
    target_include_directories(llmap_classical PUBLIC ${WFA2_INCLUDE_DIR})
endif()

# Validation module — kill-switch validation framework
add_library(llmap_validation STATIC
    validation/killswitch.cpp
    validation/end_to_end.cpp
    validation/real_reference_parse.cpp
    validation/real_reference_validate.cpp
    validation/real_reference_slurm.cpp
)

target_include_directories(llmap_validation PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_validation PUBLIC cxx_std_23)
target_link_libraries(llmap_validation PUBLIC llmap_core llmap_classical llmap_io)

# Claude agent module — LLM integration for adaptive alignment
add_library(llmap_claude_agent STATIC
    claude_agent/agent_types.cpp
    claude_agent/anthropic_client.cpp
    claude_agent/agent_session.cpp
    claude_agent/biology_prior.cpp
    claude_agent/cuda_sandbox_analyzer.cpp
    claude_agent/cuda_sandbox_parse.cpp
    claude_agent/cuda_sandbox_patterns.cpp
    claude_agent/cuda_sandbox_compile.cpp
    claude_agent/cuda_sandbox_load.cpp
    claude_agent/pipeline_agent.cpp
)

target_include_directories(llmap_claude_agent PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_claude_agent PUBLIC cxx_std_23)
target_link_libraries(llmap_claude_agent PUBLIC llmap_core)

# CLI — full CLI surface
add_executable(llmap
    cli/llmap_main.cpp
    cli/cmd_allpair.cpp
    cli/cmd_generate_synth.cpp
    cli/cmd_validate_real.cpp
    cli/cmd_align.cpp
    cli/cmd_align_args.cpp
    cli/cmd_align_llm.cpp
)
target_link_libraries(llmap PRIVATE llmap_core llmap_ai llmap_io llmap_output llmap_self_interference llmap_reference_collapse llmap_classical llmap_validation llmap_claude_agent)
