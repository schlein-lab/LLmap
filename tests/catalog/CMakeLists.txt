# LLmap — Tests for src/catalog (SegDup catalog loader).

include(GoogleTest)

add_executable(test_curated_loader
    test_curated_loader.cpp
)
target_link_libraries(test_curated_loader
    PRIVATE
    llmap_catalog
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_curated_loader)

add_executable(test_coord_lookup
    test_coord_lookup.cpp
)
target_link_libraries(test_coord_lookup
    PRIVATE
    llmap_catalog
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_coord_lookup)

add_executable(test_bulk_loader
    test_bulk_loader.cpp
)
target_link_libraries(test_bulk_loader
    PRIVATE
    llmap_catalog
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_bulk_loader)

# for_each_curated() + adapter no-IGHG-hardcode proof.
# Links llmap_classify so it can exercise SegDupCatalogAdapter.
add_executable(test_for_each_curated
    test_for_each_curated.cpp
)
target_link_libraries(test_for_each_curated
    PRIVATE
    llmap_catalog
    llmap_classify
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_for_each_curated)
