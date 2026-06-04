# LLmap — Tests for src/junction_hunter (Mode-5 NAHR junction caller).

include(GoogleTest)

add_executable(test_junction_hunter
    test_junction_hunter.cpp
)
target_link_libraries(test_junction_hunter
    PRIVATE
    llmap_junction_hunter
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_junction_hunter)
