# Anchor module — multi-source anchor catalog for Transcript-Mode.
#
# Sources: GENCODE / MANE / RefSeq / IMGT/GENE-DB / IPD-IMGT-HLA (extern-
# ref only) / FANTOM5-CAT / ChessDB / HPRC per-Hap / BRANCH bubble /
# computed cluster.
#
# Depends on llmap_core (TranscriptKind enum) and zlib (gzipped GFF3 parse).

add_library(llmap_anchor STATIC
    anchor_record.cpp
    anchor_store.cpp
    anchor_loader_gencode.cpp
    anchor_loader_misc.cpp
    cluster_anchor_promotion.cpp
)

target_include_directories(llmap_anchor PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_anchor PUBLIC cxx_std_23)

# zlib for gzopen — already a transitive dep through htslib but we declare
# explicitly so this library is buildable standalone.
find_package(ZLIB REQUIRED)
target_link_libraries(llmap_anchor PUBLIC llmap_core ZLIB::ZLIB)
