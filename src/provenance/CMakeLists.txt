# LLmap provenance — read provenance class + contamination spectrum.
#
# Every read carries exactly one provenance flag (Host or a confounder class);
# the spectrum aggregator enforces the lossless Σ-invariant. Dependency-light
# (stdlib + llmap_core) so the detector classes and the WaveCollapse EM can both
# feed it without dragging heavy deps.

add_library(llmap_provenance STATIC
    provenance_class.cpp
    contamination_spectrum.cpp
    mapping_confusion.cpp
    provenance_resolver.cpp
    pseudogene_catalog.cpp
    chimera_provenance.cpp
    te_family_catalog.cpp
    te_spread_mass.cpp
    mei_detector.cpp
    provenance_baseline.cpp
    te_population_prior.cpp
    numt_catalog.cpp
    exo_panel.cpp
)

target_include_directories(llmap_provenance PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_compile_features(llmap_provenance PUBLIC cxx_std_23)
# llmap_chimera: the Class-D chimera_provenance wrapper uses the Block-7 detector.
target_link_libraries(llmap_provenance PUBLIC llmap_core llmap_chimera)
