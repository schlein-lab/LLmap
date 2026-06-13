# LLmap — mapper module.
#
# Currently only contains the SegDup fallback hierarchy. Additional
# mapper concerns (decision logic, locus dispatch) land here later.

add_library(llmap_mapper STATIC
    fallback_types.cpp
    fallback_stages.cpp
    fallback_chain.cpp
    sufficiency.cpp
)

target_include_directories(llmap_mapper PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_mapper PUBLIC cxx_std_23)
target_link_libraries(llmap_mapper PUBLIC llmap_core)
