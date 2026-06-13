# LLmap — Mapper module tests (fallback hierarchy).

add_executable(test_fallback_stages
    test_fallback_stages.cpp
)
target_link_libraries(test_fallback_stages
    PRIVATE
    llmap_mapper
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_fallback_stages)

add_executable(test_sufficiency
    test_sufficiency.cpp
)
target_link_libraries(test_sufficiency
    PRIVATE
    llmap_mapper
    llmap_core
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_sufficiency)
