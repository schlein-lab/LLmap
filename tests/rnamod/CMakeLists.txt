# LLmap — Tests for src/rnamod (Modifications + Splicing-State + AID).

include(GoogleTest)

add_executable(test_modification_db
    test_modification_db.cpp
)
target_link_libraries(test_modification_db
    PRIVATE
    llmap_rnamod
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_modification_db)

add_executable(test_splicing_state
    test_splicing_state.cpp
)
target_link_libraries(test_splicing_state
    PRIVATE
    llmap_rnamod
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_splicing_state)

add_executable(test_aid_footprint
    test_aid_footprint.cpp
)
target_link_libraries(test_aid_footprint
    PRIVATE
    llmap_rnamod
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_aid_footprint)
