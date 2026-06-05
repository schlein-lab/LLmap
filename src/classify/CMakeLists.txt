# Classify module — IGH-mode read classifier (canonical vs chimdup IGHG4).
#
# Public surface:
#   - classify/segdup_catalog_iface.h    (read-only catalog interface used by
#                                          upstream-anchor; defines the contract)
#   - classify/segdup_catalog_adapter.h  (adapts llmap::catalog::SegDupCatalog
#                                          to ISegDupCatalog)
#   - classify/upstream_anchor.h         (Classify() entry point)

add_library(llmap_classify STATIC
    upstream_anchor.cpp
    segdup_catalog_adapter.cpp
    taxbin.cpp
)

target_include_directories(llmap_classify PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_classify PUBLIC cxx_std_23)
target_link_libraries(llmap_classify PUBLIC llmap_core llmap_catalog)
