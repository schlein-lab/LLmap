# LLmap mapping — Transcript-Mode spliced-chain joiner.
#
# Operates on top of the classical chainer's output. The classical
# chainer (src/classical/) is NOT modified — DNA-mode pipelines that
# skip this joiner produce bit-identical output to before.

add_library(llmap_mapping STATIC
    chain_spliced.cpp
)

target_include_directories(llmap_mapping PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_mapping PUBLIC cxx_std_23)
target_link_libraries(llmap_mapping PUBLIC llmap_core)
