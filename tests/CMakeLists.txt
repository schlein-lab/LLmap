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

# Embedder throughput benchmark (standalone, not using Google Benchmark)
if(LLMAP_ENABLE_BENCH)
    add_executable(bench_embedder_throughput
        bench/bench_embedder_throughput.cpp
    )
    target_link_libraries(bench_embedder_throughput
        PRIVATE
        llmap_ai
    )
endif()
