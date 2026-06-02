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

# Block 4.5b — ExpressionDb (GTEx + Tabula + HCA + HPA + Recount3).
add_executable(test_expression_db
    test_expression_db.cpp
)
target_link_libraries(test_expression_db
    PRIVATE
    llmap_fusion
    llmap_anchor
    GTest::gtest
    GTest::gtest_main
)
gtest_discover_tests(test_expression_db)
