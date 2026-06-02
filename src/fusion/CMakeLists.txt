# LLmap fusion — Multi-Signal Fusion Engine (10 factors).
#
# The differentiating layer that fuses sequence + modification +
# coverage + expression + phasing + pseudogene + junction + barcode +
# MAPQ + length signals into a single multiplicative likelihood for
# Wave-Collapse. Each factor lives in a tiny private function inside
# likelihood_factors.cpp; the public entry point is ComputeFactors.

add_library(llmap_fusion STATIC
    likelihood_factors.cpp
    expression_db.cpp
)

target_include_directories(llmap_fusion PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_fusion PUBLIC cxx_std_23)
target_link_libraries(llmap_fusion PUBLIC llmap_core llmap_anchor)
