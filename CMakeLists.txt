cmake_minimum_required(VERSION 3.28)

project(LLmap
    VERSION 0.0.1
    DESCRIPTION "Lossless. LLM-Augmented. Wave-Particle. Sequence mapper."
    HOMEPAGE_URL "https://losslessmap.com"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Build options — GPU/AI optional but architecturally assumed
option(LLMAP_ENABLE_CUDA       "Build CUDA kernels (WaveCollapse on GPU)" ON)
option(LLMAP_ENABLE_FOUNDATION "Build foundation-model inference path"    ON)
option(LLMAP_ENABLE_CLAUDE     "Build Claude-agent integration"           ON)
option(LLMAP_ENABLE_TESTS      "Build unit and integration tests"         ON)
option(LLMAP_ENABLE_BENCH      "Build benchmark harness"                  ON)
option(LLMAP_USE_NATIVE_ARCH   "Compile with -march=native (faster, less portable)" ON)

# Warnings & sanitizers
add_compile_options(
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter
)

if(LLMAP_USE_NATIVE_ARCH AND CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-march=native -mtune=native)
endif()

# CUDA — enable language if requested and available
if(LLMAP_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        set(CMAKE_CUDA_ARCHITECTURES 80 86 89 90)  # Ampere, Ada, Hopper
        message(STATUS "CUDA enabled: ${CMAKE_CUDA_COMPILER}")
    else()
        message(WARNING "CUDA requested but no compiler found — falling back to CPU paths only")
        set(LLMAP_ENABLE_CUDA OFF)
    endif()
endif()

# Subdirectories
add_subdirectory(src)

if(LLMAP_ENABLE_TESTS)
    enable_testing()
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    )
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()

# Status summary
message(STATUS "")
message(STATUS "LLmap configuration:")
message(STATUS "  Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard    : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  CUDA            : ${LLMAP_ENABLE_CUDA}")
message(STATUS "  Foundation Model: ${LLMAP_ENABLE_FOUNDATION}")
message(STATUS "  Claude Agent    : ${LLMAP_ENABLE_CLAUDE}")
message(STATUS "  Tests           : ${LLMAP_ENABLE_TESTS}")
message(STATUS "  Benchmarks      : ${LLMAP_ENABLE_BENCH}")
message(STATUS "  Native arch     : ${LLMAP_USE_NATIVE_ARCH}")
message(STATUS "")
