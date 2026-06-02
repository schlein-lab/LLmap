# LLmap — Tests for src/wavecollapse (Transcript-Mode spliced likelihood).

include(GoogleTest)

add_executable(test_spliced_likelihood
    test_spliced_likelihood.cpp
)
target_link_libraries(test_spliced_likelihood
    PRIVATE
    llmap_wavecollapse
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_spliced_likelihood)
