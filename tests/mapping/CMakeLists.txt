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
