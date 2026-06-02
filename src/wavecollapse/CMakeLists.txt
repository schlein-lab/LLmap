# LLmap wavecollapse — Transcript-Mode spliced likelihood + kernel.
#
# Standalone module; does NOT modify the existing em_iterator (DNA-mode
# WaveCollapse stays bit-identical). The Transcript-Mode dispatcher
# (Block 6) selects these functions at the L(r|b) / K(b,b') call sites.

add_library(llmap_wavecollapse STATIC
    spliced_likelihood.cpp
)

target_include_directories(llmap_wavecollapse PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_wavecollapse PUBLIC cxx_std_23)
target_link_libraries(llmap_wavecollapse PUBLIC
    llmap_core llmap_anchor llmap_annot llmap_fusion)
