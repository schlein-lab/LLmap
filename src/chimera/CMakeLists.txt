# LLmap chimera — ChimeraDetector + VdjLocusMask.

add_library(llmap_chimera STATIC
    chimera_detector.cpp
)

target_include_directories(llmap_chimera PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_chimera PUBLIC cxx_std_23)
target_link_libraries(llmap_chimera PUBLIC llmap_core)
