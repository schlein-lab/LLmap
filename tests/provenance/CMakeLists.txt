# LLmap — Tests for src/provenance.

include(GoogleTest)

add_executable(test_provenance
    test_provenance.cpp
)
target_link_libraries(test_provenance
    PRIVATE
    llmap_provenance
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_provenance)

add_executable(test_mapping_confusion
    test_mapping_confusion.cpp
)
target_link_libraries(test_mapping_confusion
    PRIVATE
    llmap_provenance
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_mapping_confusion)

add_executable(test_provenance_resolver
    test_provenance_resolver.cpp
)
target_link_libraries(test_provenance_resolver
    PRIVATE
    llmap_provenance
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_provenance_resolver)

add_executable(test_pseudogene_catalog
    test_pseudogene_catalog.cpp
)
target_link_libraries(test_pseudogene_catalog
    PRIVATE
    llmap_provenance
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_pseudogene_catalog)

add_executable(test_chimera_provenance
    test_chimera_provenance.cpp
)
target_link_libraries(test_chimera_provenance
    PRIVATE
    llmap_provenance
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_chimera_provenance)
