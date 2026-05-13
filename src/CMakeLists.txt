add_library(llmap_core STATIC
    core/alignment_record.cpp
    core/bucket_pyramid.cpp
    core/wave_state.cpp
)

target_include_directories(llmap_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_core PUBLIC cxx_std_23)

# CLI placeholder — will grow into full CLI surface in later phases
add_executable(llmap cli/llmap_main.cpp)
target_link_libraries(llmap PRIVATE llmap_core)
