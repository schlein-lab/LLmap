# LLmap — Tests for src/chimera (Chimera detector + VDJ-mask).

include(GoogleTest)

add_executable(test_chimera_detector
    test_chimera_detector.cpp
)
target_link_libraries(test_chimera_detector
    PRIVATE
    llmap_chimera
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_chimera_detector)
