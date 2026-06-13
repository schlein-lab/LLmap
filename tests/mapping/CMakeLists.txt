# LLmap — Tests for src/mapping (Spliced-chain joiner).

include(GoogleTest)

add_executable(test_chain_spliced
    test_chain_spliced.cpp
)
target_link_libraries(test_chain_spliced
    PRIVATE
    llmap_mapping
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_chain_spliced)

add_executable(test_transcript_stage
    test_transcript_stage.cpp
)
target_link_libraries(test_transcript_stage
    PRIVATE
    llmap_mapping
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_transcript_stage)

add_executable(test_splice_snap
    test_splice_snap.cpp
)
target_link_libraries(test_splice_snap
    PRIVATE
    llmap_mapping
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_splice_snap)

add_executable(test_splice_determinism
    test_splice_determinism.cpp
)
target_link_libraries(test_splice_determinism
    PRIVATE
    llmap_mapping
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_splice_determinism)
