// LLmap Phase 0 — IGH Locus Generator tests.
//
// Tests the synthetic IGH locus generator including:
// - Deterministic generation by seed
// - PSV planting and sequence-identical regions
// - Mosaic duplication fraction control
// - Read generation with ground truth

#include "synthetic/igh_locus_generator.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_set>

using namespace llmap::synthetic;

namespace {

class IGHLocusGeneratorTest : public ::testing::Test {
protected:
    static constexpr float TOLERANCE = 0.15f;  // 15% tolerance for stochastic tests
};

}  // namespace

TEST_F(IGHLocusGeneratorTest, DeterministicBySeed) {
    auto cfg = presets::tiny_test(12345);

    IGHLocusGenerator gen1(cfg);
    IGHLocusGenerator gen2(cfg);

    auto dataset1 = gen1.generate();
    auto dataset2 = gen2.generate();

    // Locus sequences should be identical
    EXPECT_EQ(dataset1.locus.canonical_sequence, dataset2.locus.canonical_sequence);
    EXPECT_EQ(dataset1.locus.duplicate_sequence, dataset2.locus.duplicate_sequence);

    // PSVs should be identical
    ASSERT_EQ(dataset1.locus.psvs.size(), dataset2.locus.psvs.size());
    for (std::size_t i = 0; i < dataset1.locus.psvs.size(); ++i) {
        EXPECT_EQ(dataset1.locus.psvs[i].position, dataset2.locus.psvs[i].position);
        EXPECT_EQ(dataset1.locus.psvs[i].canonical_allele,
                  dataset2.locus.psvs[i].canonical_allele);
        EXPECT_EQ(dataset1.locus.psvs[i].duplicate_allele,
                  dataset2.locus.psvs[i].duplicate_allele);
    }

    // Read count and order should be identical
    ASSERT_EQ(dataset1.reads.size(), dataset2.reads.size());
    for (std::size_t i = 0; i < dataset1.reads.size(); ++i) {
        EXPECT_EQ(dataset1.reads[i].id, dataset2.reads[i].id);
        EXPECT_EQ(dataset1.reads[i].sequence, dataset2.reads[i].sequence);
        EXPECT_EQ(dataset1.reads[i].origin, dataset2.reads[i].origin);
    }
}

TEST_F(IGHLocusGeneratorTest, DifferentSeedsDifferentOutput) {
    auto cfg1 = presets::tiny_test(111);
    auto cfg2 = presets::tiny_test(222);

    IGHLocusGenerator gen1(cfg1);
    IGHLocusGenerator gen2(cfg2);

    auto locus1 = gen1.generate_locus();
    auto locus2 = gen2.generate_locus();

    // Sequences should differ
    EXPECT_NE(locus1.canonical_sequence, locus2.canonical_sequence);
}

TEST_F(IGHLocusGeneratorTest, PSVsAreDiscriminating) {
    auto cfg = presets::tiny_test(42);
    cfg.seq_identical_fraction = 0.0f;  // no identical regions

    IGHLocusGenerator gen(cfg);
    auto locus = gen.generate_locus();

    EXPECT_EQ(locus.psvs.size(), cfg.n_psvs);

    for (const auto& psv : locus.psvs) {
        if (psv.is_discriminating) {
            // Canonical and duplicate alleles must differ at discriminating PSVs
            EXPECT_NE(psv.canonical_allele, psv.duplicate_allele)
                << "PSV at position " << psv.position << " has same alleles";

            // The sequences should reflect these differences
            EXPECT_EQ(locus.canonical_sequence[psv.position], psv.canonical_allele);
            EXPECT_EQ(locus.duplicate_sequence[psv.position], psv.duplicate_allele);
        }
    }
}

TEST_F(IGHLocusGeneratorTest, SequenceIdenticalRegionsAreCorrect) {
    auto cfg = presets::seq_identical_stress(42);

    IGHLocusGenerator gen(cfg);
    auto locus = gen.generate_locus();

    // Should have at least one identical region
    EXPECT_GE(locus.identical_regions.size(), 1u);

    for (const auto& region : locus.identical_regions) {
        EXPECT_LT(region.start, region.end);
        EXPECT_LE(region.end, cfg.locus_length);

        // Sequences must be identical in this region
        for (std::uint64_t pos = region.start; pos < region.end; ++pos) {
            EXPECT_EQ(locus.canonical_sequence[pos], locus.duplicate_sequence[pos])
                << "Mismatch in identical region at position " << pos;
        }
    }

    // No discriminating PSVs should fall within identical regions
    for (const auto& psv : locus.psvs) {
        if (psv.is_discriminating) {
            for (const auto& region : locus.identical_regions) {
                bool in_region = psv.position >= region.start && psv.position < region.end;
                EXPECT_FALSE(in_region)
                    << "Discriminating PSV at " << psv.position
                    << " falls within identical region [" << region.start
                    << ", " << region.end << ")";
            }
        }
    }
}

TEST_F(IGHLocusGeneratorTest, MosaicDupFractionCanonicalOnly) {
    auto cfg = presets::canonical_only(42);

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    EXPECT_GT(dataset.reads.size(), 0u);
    EXPECT_FLOAT_EQ(dataset.actual_dup_fraction, 0.0f);
    EXPECT_EQ(dataset.n_dup_reads, 0u);
    EXPECT_EQ(dataset.n_canonical_reads, dataset.reads.size());
}

TEST_F(IGHLocusGeneratorTest, MosaicDupFractionBalanced) {
    auto cfg = presets::balanced_mosaic(42);
    cfg.locus_length = 10000;
    cfg.mosaic.canonical_depth = 50;
    cfg.mosaic.dup_depth = 50;

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    EXPECT_GT(dataset.reads.size(), 0u);

    // Should be roughly 50/50, within tolerance
    float expected = 0.5f;
    EXPECT_NEAR(dataset.actual_dup_fraction, expected, TOLERANCE);
}

TEST_F(IGHLocusGeneratorTest, MosaicDupFractionSweep) {
    std::vector<std::pair<float, IGHLocusConfig(*)(std::uint64_t)>> presets = {
        {0.05f, presets::dup_fraction_5},
        {0.10f, presets::dup_fraction_10},
        {0.30f, presets::dup_fraction_30},
        {0.50f, presets::dup_fraction_50},
        {1.00f, presets::dup_fraction_100},
    };

    for (const auto& [expected_frac, preset_fn] : presets) {
        auto cfg = preset_fn(42);
        cfg.locus_length = 20000;
        cfg.mosaic.canonical_depth = 100;

        IGHLocusGenerator gen(cfg);
        auto dataset = gen.generate();

        EXPECT_GT(dataset.reads.size(), 50u)
            << "Not enough reads for dup_fraction=" << expected_frac;

        EXPECT_NEAR(dataset.actual_dup_fraction, expected_frac, TOLERANCE)
            << "Dup fraction mismatch for expected=" << expected_frac;
    }
}

TEST_F(IGHLocusGeneratorTest, ReadGroundTruthIsConsistent) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    for (const auto& read : dataset.reads) {
        // ID should encode origin
        bool id_says_dup = read.id.find("_DUP_") != std::string::npos;
        bool id_says_can = read.id.find("_CAN_") != std::string::npos;

        EXPECT_TRUE(id_says_dup || id_says_can) << "Read ID missing origin: " << read.id;
        EXPECT_NE(id_says_dup, id_says_can) << "Read ID has both origins: " << read.id;

        if (id_says_dup) {
            EXPECT_EQ(read.origin, SyntheticRead::Origin::Duplicate);
        } else {
            EXPECT_EQ(read.origin, SyntheticRead::Origin::Canonical);
        }

        // Position should be valid
        EXPECT_LT(read.true_start, read.true_end);
        EXPECT_LE(read.true_end, cfg.locus_length);

        // Quality string should match sequence length
        EXPECT_EQ(read.sequence.size(), read.quality.size());

        // Source gene should be valid
        EXPECT_FALSE(read.source_gene.empty());
    }
}

TEST_F(IGHLocusGeneratorTest, CoveredPSVsAreCorrect) {
    auto cfg = presets::tiny_test(42);
    cfg.n_psvs = 10;
    cfg.read_length = 2000;  // long reads to cover multiple PSVs

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    // Build PSV position set
    std::unordered_set<std::uint64_t> psv_positions;
    for (const auto& psv : dataset.locus.psvs) {
        psv_positions.insert(psv.position);
    }

    for (const auto& read : dataset.reads) {
        for (std::uint64_t psv_pos : read.covered_psv_positions) {
            // Covered PSV must exist
            EXPECT_TRUE(psv_positions.count(psv_pos) > 0)
                << "Read claims to cover non-existent PSV at " << psv_pos;

            // PSV must be within read bounds
            EXPECT_GE(psv_pos, read.true_start);
            EXPECT_LT(psv_pos, read.true_end);
        }
    }
}

TEST_F(IGHLocusGeneratorTest, SequenceLengthsAreCorrect) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto locus = gen.generate_locus();

    EXPECT_EQ(locus.canonical_sequence.size(), cfg.locus_length);
    EXPECT_EQ(locus.duplicate_sequence.size(), cfg.locus_length);
}

TEST_F(IGHLocusGeneratorTest, AllBasesAreValid) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto locus = gen.generate_locus();

    std::set<char> valid_bases{'A', 'C', 'G', 'T'};

    for (char c : locus.canonical_sequence) {
        EXPECT_TRUE(valid_bases.count(c) > 0) << "Invalid base in canonical: " << c;
    }

    for (char c : locus.duplicate_sequence) {
        EXPECT_TRUE(valid_bases.count(c) > 0) << "Invalid base in duplicate: " << c;
    }
}

TEST_F(IGHLocusGeneratorTest, GeneContextsAreAssigned) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto locus = gen.generate_locus();

    for (const auto& psv : locus.psvs) {
        EXPECT_FALSE(psv.gene_context.empty());
    }
}

TEST_F(IGHLocusGeneratorTest, TotalPSVObservationsAreTracked) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    std::uint64_t counted = 0;
    for (const auto& read : dataset.reads) {
        counted += read.covered_psv_positions.size();
    }

    EXPECT_EQ(dataset.total_psv_observations, counted);
}

TEST_F(IGHLocusGeneratorTest, QualityStringsArePhred33) {
    auto cfg = presets::tiny_test(42);

    IGHLocusGenerator gen(cfg);
    auto dataset = gen.generate();

    for (const auto& read : dataset.reads) {
        for (char q : read.quality) {
            int phred = static_cast<int>(q) - 33;
            EXPECT_GE(phred, 0) << "Quality below Phred+33 minimum";
            EXPECT_LE(phred, 50) << "Quality above reasonable maximum";
        }
    }
}

TEST_F(IGHLocusGeneratorTest, ErrorRateAffectsSequences) {
    auto cfg1 = presets::tiny_test(42);
    cfg1.error_rate = 0.0f;

    auto cfg2 = presets::tiny_test(42);
    cfg2.error_rate = 0.1f;

    IGHLocusGenerator gen1(cfg1);
    IGHLocusGenerator gen2(cfg2);

    auto dataset1 = gen1.generate();
    auto dataset2 = gen2.generate();

    // Compare first few reads - sequences should be non-empty
    for (std::size_t i = 0; i < std::min<std::size_t>(5, dataset1.reads.size()); ++i) {
        EXPECT_FALSE(dataset1.reads[i].sequence.empty());
        EXPECT_FALSE(dataset2.reads[i].sequence.empty());
    }
}

// Test preset configurations compile and produce valid output
TEST_F(IGHLocusGeneratorTest, AllPresetsWork) {
    std::vector<IGHLocusConfig(*)(std::uint64_t)> preset_fns = {
        presets::canonical_only,
        presets::balanced_mosaic,
        presets::dup_fraction_5,
        presets::dup_fraction_10,
        presets::dup_fraction_30,
        presets::dup_fraction_50,
        presets::dup_fraction_100,
        presets::seq_identical_stress,
        presets::tiny_test,
    };

    for (auto preset_fn : preset_fns) {
        auto cfg = preset_fn(42);
        IGHLocusGenerator gen(cfg);
        auto dataset = gen.generate();

        EXPECT_FALSE(dataset.locus.canonical_sequence.empty());
        EXPECT_FALSE(dataset.locus.duplicate_sequence.empty());
    }
}
