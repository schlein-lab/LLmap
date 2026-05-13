// LLmap — Synthetic IGH-locus generator preset configurations.

#include "synthetic/igh_locus_generator.h"

namespace llmap::synthetic::presets {

IGHLocusConfig canonical_only(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.0f;
    cfg.mosaic.canonical_depth = 30;
    return cfg;
}

IGHLocusConfig balanced_mosaic(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.5f;
    cfg.mosaic.canonical_depth = 15;
    cfg.mosaic.dup_depth = 15;
    return cfg;
}

IGHLocusConfig dup_fraction_5(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.05f;
    return cfg;
}

IGHLocusConfig dup_fraction_10(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.10f;
    return cfg;
}

IGHLocusConfig dup_fraction_30(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.30f;
    return cfg;
}

IGHLocusConfig dup_fraction_50(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 0.50f;
    return cfg;
}

IGHLocusConfig dup_fraction_100(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.mosaic.dup_fraction = 1.0f;
    return cfg;
}

IGHLocusConfig seq_identical_stress(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.seq_identical_fraction = 0.5f;
    cfg.mosaic.dup_fraction = 0.5f;
    cfg.n_psvs = 10;
    return cfg;
}

IGHLocusConfig tiny_test(std::uint64_t seed) {
    IGHLocusConfig cfg;
    cfg.seed = seed;
    cfg.locus_length = 5000;
    cfg.n_psvs = 5;
    cfg.read_length = 1000;
    cfg.read_length_stddev = 100.0f;
    cfg.mosaic.canonical_depth = 10;
    cfg.mosaic.dup_fraction = 0.3f;
    return cfg;
}

}  // namespace llmap::synthetic::presets
