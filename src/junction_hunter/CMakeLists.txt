# LLmap junction_hunter — Mode-5 NAHR breakpoint detector.
#
# Two-stage pipeline:
#   Stage-1 ("coarse_locate") — k=51 anchor counting per NAHR-pair.
#   Stage-2 ("multi_k_consensus") — five-k consensus + outside-in
#     monotonicity + breakpoint detection.
#
# No minimap2 dependency at any stage (see feedback_minimap2_nahr_forbidden).

add_library(llmap_junction_hunter STATIC
    junction_hunter_types.cpp
    nahr_pair_loader.cpp
    anchor_table_builder.cpp
    read_tiler.cpp
    pair_kmer_index.cpp
    consensus_caller.cpp
)

target_include_directories(llmap_junction_hunter PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_junction_hunter PUBLIC cxx_std_23)
target_link_libraries(llmap_junction_hunter PUBLIC llmap_core)
