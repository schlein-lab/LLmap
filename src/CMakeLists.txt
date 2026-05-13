add_library(llmap_core STATIC
    core/alignment_record.cpp
    core/bucket_pyramid.cpp
    core/wave_state.cpp
    synthetic/igh_locus_generator.cpp
)

target_include_directories(llmap_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_core PUBLIC cxx_std_23)

# AI module — foundation model embedders
add_library(llmap_ai STATIC
    ai/foundation_embedder.cpp
    ai/bucket_embedder.cpp
)

target_include_directories(llmap_ai PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_ai PUBLIC cxx_std_23)
target_link_libraries(llmap_ai PUBLIC llmap_core)

# Link ONNX Runtime if available
if(LLMAP_ONNXRUNTIME_AVAILABLE)
    target_compile_definitions(llmap_ai PRIVATE LLMAP_HAS_ONNXRUNTIME)
    if(TARGET onnxruntime::onnxruntime)
        target_link_libraries(llmap_ai PUBLIC onnxruntime::onnxruntime)
    elseif(ONNXRUNTIME_LIBRARY)
        target_link_libraries(llmap_ai PUBLIC ${ONNXRUNTIME_LIBRARY})
        target_include_directories(llmap_ai PUBLIC ${ONNXRUNTIME_INCLUDE_DIR})
    endif()
endif()

# CLI placeholder — will grow into full CLI surface in later phases
add_executable(llmap cli/llmap_main.cpp)
target_link_libraries(llmap PRIVATE llmap_core llmap_ai)
