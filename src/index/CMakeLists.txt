# LLmap index — Transcript-Mode reverse k-mer index.
#
# Two translation units, kept under the 400-LOC soft cap:
#   transcript_kmer_index.cpp        — storage / lookup / hashing
#   transcript_kmer_index_build.cpp  — BuildFromAnchorStore (6 origin tables)

add_library(llmap_index STATIC
    transcript_kmer_index.cpp
    transcript_kmer_index_build.cpp
)

target_include_directories(llmap_index PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_index PUBLIC cxx_std_23)
target_link_libraries(llmap_index PUBLIC llmap_core llmap_annot)
