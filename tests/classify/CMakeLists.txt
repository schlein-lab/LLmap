# Classify module — unit tests.

add_executable(test_upstream_anchor
    test_upstream_anchor.cpp
)
target_link_libraries(test_upstream_anchor
    PRIVATE
    llmap_classify
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_upstream_anchor)
