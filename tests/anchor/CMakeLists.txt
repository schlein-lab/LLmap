# LLmap — Tests for src/anchor (Multi-Source Anchor Store).

include(GoogleTest)

add_executable(test_anchor_store
    test_anchor_store.cpp
)
target_link_libraries(test_anchor_store
    PRIVATE
    llmap_anchor
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_anchor_store)

# Block 4 — cluster-anchor promotion (Stage-1 → AnchorStore).
add_executable(test_cluster_anchor_promotion
    test_cluster_anchor_promotion.cpp
)
target_link_libraries(test_cluster_anchor_promotion
    PRIVATE
    llmap_anchor
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_cluster_anchor_promotion)
