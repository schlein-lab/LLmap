# LLmap — SegDup Catalog module.
# Provides nlohmann/json-backed loaders for catalog/curated/*.json (T1)
# and catalog/bulk/*.jsonl (T2). Coordinate + locus-id lookup.

include(FetchContent)

# Try local config first (CMAKE_PREFIX_PATH / system / conda), then fall
# back to FetchContent from GitHub. CMake 3.24+ FIND_PACKAGE_ARGS makes
# this seamless.
find_package(nlohmann_json 3.10 QUIET)
if(NOT nlohmann_json_FOUND)
    FetchContent_Declare(
        nlohmann_json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
        URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
    )
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install OFF CACHE INTERNAL "")
    FetchContent_MakeAvailable(nlohmann_json)
endif()

add_library(llmap_catalog STATIC
    segdup_catalog.cpp
    segdup_catalog_parse.cpp
    bulk_loader.cpp
)

target_include_directories(llmap_catalog
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_catalog PUBLIC cxx_std_23)

target_link_libraries(llmap_catalog
    PUBLIC
        nlohmann_json::nlohmann_json
)

# Optional: expose the catalog source dir to consumers via a definition,
# so tests / CLI can locate the canonical T1 JSONs without env vars.
target_compile_definitions(llmap_catalog
    PUBLIC
        LLMAP_CATALOG_CURATED_DIR="${CMAKE_SOURCE_DIR}/catalog/curated"
        LLMAP_CATALOG_BULK_DIR="${CMAKE_SOURCE_DIR}/catalog/bulk"
)
