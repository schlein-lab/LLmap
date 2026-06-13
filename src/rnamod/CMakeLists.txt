# LLmap rnamod — RNA modification + splicing-state + AID-footprint layer.
#
# Three orthogonal annotators that ride on top of any successful
# alignment. None of them call modifications from raw signal; they
# import + preserve external tool output and turn it into structured
# AlignmentRecord payloads.

add_library(llmap_rnamod STATIC
    modification_db.cpp
    splicing_state.cpp
    aid_footprint.cpp
    damage_profile.cpp
)

target_include_directories(llmap_rnamod PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_rnamod PUBLIC cxx_std_23)
target_link_libraries(llmap_rnamod PUBLIC llmap_core)
