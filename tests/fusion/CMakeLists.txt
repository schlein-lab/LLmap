# LLmap — Tests for src/fusion (Multi-Signal Fusion Engine).

include(GoogleTest)

add_executable(test_likelihood_factors
    test_likelihood_factors.cpp
)
target_link_libraries(test_likelihood_factors
    PRIVATE
    llmap_fusion
    llmap_anchor
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_likelihood_factors)
