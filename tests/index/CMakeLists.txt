# LLmap — Tests for src/index (TranscriptKmerIndex).

include(GoogleTest)

add_executable(test_transcript_kmer_index
    test_transcript_kmer_index.cpp
)
target_link_libraries(test_transcript_kmer_index
    PRIVATE
    llmap_index
    llmap_anchor
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_transcript_kmer_index)
