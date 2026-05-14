include(GoogleTest)

add_executable(test_alignment_record
    unit/test_alignment_record.cpp
)
target_link_libraries(test_alignment_record
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_alignment_record)

add_executable(test_bucket_pyramid
    unit/test_bucket_pyramid.cpp
)
target_link_libraries(test_bucket_pyramid
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_bucket_pyramid)

add_executable(test_wave_state
    unit/test_wave_state.cpp
)
target_link_libraries(test_wave_state
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_wave_state)

add_executable(test_igh_locus_generator
    unit/test_igh_locus_generator.cpp
)
target_link_libraries(test_igh_locus_generator
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_igh_locus_generator)

add_executable(test_foundation_embedder
    unit/test_foundation_embedder.cpp
)
target_link_libraries(test_foundation_embedder
    PRIVATE
    llmap_ai
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_foundation_embedder)

add_executable(test_bucket_embedder
    unit/test_bucket_embedder.cpp
)
target_link_libraries(test_bucket_embedder
    PRIVATE
    llmap_ai
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_bucket_embedder)

add_executable(test_faiss_wrapper
    unit/test_faiss_wrapper.cpp
)
target_link_libraries(test_faiss_wrapper
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_faiss_wrapper)

add_executable(test_similarity_graph
    unit/test_similarity_graph.cpp
)
target_link_libraries(test_similarity_graph
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_similarity_graph)

add_executable(test_leiden_clustering
    unit/test_leiden_clustering.cpp
)
target_link_libraries(test_leiden_clustering
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_leiden_clustering)

add_executable(test_self_wavecollapse
    unit/test_self_wavecollapse.cpp
)
target_link_libraries(test_self_wavecollapse
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_self_wavecollapse)

add_executable(test_cluster_rep
    unit/test_cluster_rep.cpp
)
target_link_libraries(test_cluster_rep
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cluster_rep)

# IO module tests
add_executable(test_fastq_reader
    unit/test_fastq_reader.cpp
)
target_link_libraries(test_fastq_reader
    PRIVATE
    llmap_io
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_fastq_reader)

add_executable(test_mmap_fasta
    unit/test_mmap_fasta.cpp
)
target_link_libraries(test_mmap_fasta
    PRIVATE
    llmap_io
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_mmap_fasta)

# Output module tests
add_executable(test_cluster_writer
    unit/test_cluster_writer.cpp
)
target_link_libraries(test_cluster_writer
    PRIVATE
    llmap_output
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cluster_writer)

add_executable(test_bam_writer
    unit/test_bam_writer.cpp
)
target_link_libraries(test_bam_writer
    PRIVATE
    llmap_output
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_bam_writer)

add_executable(test_parquet_writer
    unit/test_parquet_writer.cpp
)
target_link_libraries(test_parquet_writer
    PRIVATE
    llmap_output
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_parquet_writer)

add_executable(test_parquet_reader
    unit/test_parquet_reader.cpp
)
target_link_libraries(test_parquet_reader
    PRIVATE
    llmap_output
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_parquet_reader)

# Allpair pipeline tests
add_executable(test_allpair_pipeline
    unit/test_allpair_pipeline.cpp
)
target_link_libraries(test_allpair_pipeline
    PRIVATE
    llmap_self_interference
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_allpair_pipeline)

# CLI tests (integration-style, uses subprocess)
add_executable(test_llmap_cli
    unit/test_llmap_cli.cpp
)
target_link_libraries(test_llmap_cli
    PRIVATE
    GTest::gtest
    GTest::gtest_main
)
# CLI tests need the llmap binary to exist
add_dependencies(test_llmap_cli llmap)
gtest_discover_tests(test_llmap_cli)

# Reference collapse module tests
add_executable(test_reference_index
    unit/test_reference_index.cpp
)
target_link_libraries(test_reference_index
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_reference_index)

add_executable(test_em_iterator
    unit/test_em_iterator.cpp
)
target_link_libraries(test_em_iterator
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_em_iterator)

add_executable(test_collapse_check
    unit/test_collapse_check.cpp
)
target_link_libraries(test_collapse_check
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_collapse_check)

add_executable(test_refinement
    unit/test_refinement.cpp
)
target_link_libraries(test_refinement
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_refinement)

add_executable(test_member_propagation
    unit/test_member_propagation.cpp
)
target_link_libraries(test_member_propagation
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_member_propagation)

add_executable(test_stage2_pipeline
    unit/test_stage2_pipeline.cpp
)
target_link_libraries(test_stage2_pipeline
    PRIVATE
    llmap_reference_collapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_stage2_pipeline)

# Classical module tests
add_executable(test_minimizer_index
    unit/test_minimizer_index.cpp
)
target_link_libraries(test_minimizer_index
    PRIVATE
    llmap_classical
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_minimizer_index)

add_executable(test_chain
    unit/test_chain.cpp
)
target_link_libraries(test_chain
    PRIVATE
    llmap_classical
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_chain)

add_executable(test_wfa2_aligner
    unit/test_wfa2_aligner.cpp
)
target_link_libraries(test_wfa2_aligner
    PRIVATE
    llmap_classical
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_wfa2_aligner)

add_executable(test_classical_pipeline
    unit/test_classical_pipeline.cpp
)
target_link_libraries(test_classical_pipeline
    PRIVATE
    llmap_classical
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_classical_pipeline)

# Validation module tests
add_executable(test_killswitch
    unit/test_killswitch.cpp
)
target_link_libraries(test_killswitch
    PRIVATE
    llmap_validation
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_killswitch)

add_executable(test_end_to_end
    unit/test_end_to_end.cpp
)
target_link_libraries(test_end_to_end
    PRIVATE
    llmap_validation
    llmap_classical
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_end_to_end)

add_executable(test_real_reference
    unit/test_real_reference.cpp
)
target_link_libraries(test_real_reference
    PRIVATE
    llmap_validation
    llmap_classical
    llmap_io
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_real_reference)

# Profiler tests
add_executable(test_profiler
    unit/test_profiler.cpp
)
target_link_libraries(test_profiler
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_profiler)

# Arena allocator tests
add_executable(test_arena
    unit/test_arena.cpp
)
target_link_libraries(test_arena
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_arena)

# SIMD utilities tests
add_executable(test_simd
    unit/test_simd.cpp
)
target_link_libraries(test_simd
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_simd)

# Thread pool tests
add_executable(test_thread_pool
    unit/test_thread_pool.cpp
)
target_link_libraries(test_thread_pool
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_thread_pool)

# Cache layout tests
add_executable(test_cache_layout
    unit/test_cache_layout.cpp
)
target_link_libraries(test_cache_layout
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cache_layout)

# Logging tests
add_executable(test_logging
    unit/test_logging.cpp
)
target_link_libraries(test_logging
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_logging)

# Error handling tests
add_executable(test_error
    unit/test_error.cpp
)
target_link_libraries(test_error
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_error)

# Configuration tests
add_executable(test_config
    unit/test_config.cpp
)
target_link_libraries(test_config
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_config)

# Version tests
add_executable(test_version
    unit/test_version.cpp
)
target_link_libraries(test_version
    PRIVATE
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
# Version CLI tests need the llmap binary to exist
add_dependencies(test_version llmap)
gtest_discover_tests(test_version)

# Embedder throughput benchmark (standalone, not using Google Benchmark)
if(LLMAP_ENABLE_BENCH)
    add_executable(bench_embedder_throughput
        bench/bench_embedder_throughput.cpp
    )
    target_link_libraries(bench_embedder_throughput
        PRIVATE
        llmap_ai
    )

    add_executable(bench_classical_pipeline
        bench/bench_classical_pipeline.cpp
    )
    target_link_libraries(bench_classical_pipeline
        PRIVATE
        llmap_classical
        llmap_core
    )
    target_compile_definitions(bench_classical_pipeline PRIVATE LLMAP_ENABLE_PROFILING)
endif()

# ========== Claude Agent Module Tests ==========

add_executable(test_agent_types
    unit/test_agent_types.cpp
)
target_link_libraries(test_agent_types
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_agent_types)

add_executable(test_anthropic_client
    unit/test_anthropic_client.cpp
)
target_link_libraries(test_anthropic_client
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_anthropic_client)

add_executable(test_agent_session
    unit/test_agent_session.cpp
)
target_link_libraries(test_agent_session
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_agent_session)

add_executable(test_biology_prior
    unit/test_biology_prior.cpp
)
target_link_libraries(test_biology_prior
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_biology_prior)

add_executable(test_cuda_sandbox
    unit/test_cuda_sandbox.cpp
)
target_link_libraries(test_cuda_sandbox
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cuda_sandbox)

add_executable(test_pipeline_agent
    unit/test_pipeline_agent.cpp
)
target_link_libraries(test_pipeline_agent
    PRIVATE
    llmap_claude_agent
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_pipeline_agent)

# ========== PSV Module Tests ==========

add_executable(test_psv
    unit/test_psv.cpp
)
target_link_libraries(test_psv
    PRIVATE
    llmap_psv
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_psv)

# ========== Single-cell Module Tests ==========

add_executable(test_cb_preservation
    unit/test_cb_preservation.cpp
)
target_link_libraries(test_cb_preservation
    PRIVATE
    llmap_singlecell
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cb_preservation)

add_executable(test_per_cell_paralog
    unit/test_per_cell_paralog.cpp
)
target_link_libraries(test_per_cell_paralog
    PRIVATE
    llmap_singlecell
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_per_cell_paralog)

add_executable(test_sc_paralog_qc
    unit/test_sc_paralog_qc.cpp
)
target_link_libraries(test_sc_paralog_qc
    PRIVATE
    llmap_singlecell
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_sc_paralog_qc)

# ========== Integration Tests ==========
# These test end-to-end workflows across multiple modules

add_executable(test_align_integration
    integration/test_align_integration.cpp
)
target_link_libraries(test_align_integration
    PRIVATE
    llmap_classical
    llmap_io
    llmap_output
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_align_integration)
