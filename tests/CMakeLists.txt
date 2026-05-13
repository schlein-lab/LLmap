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
