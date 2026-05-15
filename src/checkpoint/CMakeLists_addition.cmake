# LLmap -- Layer 3 checkpoint library.
#
# Insert this snippet into src/CMakeLists.txt (alongside the other static
# libraries). The library has no extra third-party dependencies. The
# included SHA-256 implementation is self-contained, so no OpenSSL link
# is required.

add_library(llmap_checkpoint STATIC
    checkpoint/checkpoint_cache.cpp
    checkpoint/checkpoint_dispatcher.cpp
    checkpoint/checkpoint_prompts.cpp
    checkpoint/checkpoint_tools.cpp
)

target_include_directories(llmap_checkpoint PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(llmap_checkpoint PUBLIC cxx_std_23)
target_link_libraries(llmap_checkpoint PUBLIC
    llmap_core
    llmap_annot
    llmap_claude_agent
)

# Then add `llmap_checkpoint` to the `target_link_libraries(llmap PRIVATE ...)`
# list for the `llmap` executable.
