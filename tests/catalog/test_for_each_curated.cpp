// LLmap — Unit tests for SegDupCatalog::for_each_curated.
//
// The for_each_curated() API replaces the IGHG-trio hard-code that the
// IGH classify adapter previously carried. These tests pin down three
// things:
//   1. Iteration visits exactly the curated entries, no T2 bulk leakage.
//   2. The predicate-filtered overload composes correctly.
//   3. The classify::SegDupCatalogAdapter, after the refactor, finds a
//      non-IGHG haplotype_class — proving the IGHG-only hard-code is gone.
//
// We synthesise the catalog from in-memory entries via SegDupCatalog::add_entry
// so the tests stay independent of disk-resident catalog/curated/ contents.

#include <gtest/gtest.h>

#include "catalog/segdup_catalog.h"
#include "classify/segdup_catalog_adapter.h"

#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using llmap::catalog::GenomicCoords;
using llmap::catalog::SegDupCatalog;
using llmap::catalog::SegDupCatalogEntry;

// Build a minimal valid T1 entry with the fields the API actually reads.
// We deliberately use NON-IGHG identifiers so the classify-adapter test
// below can prove the old IGHG-trio hard-code is no longer in play.
SegDupCatalogEntry MakeCurated(std::string locus_id,
                                std::string haplotype_class,
                                std::string_view chrom,
                                std::int64_t start,
                                std::int64_t end) {
    SegDupCatalogEntry e;
    e.tier             = SegDupCatalogEntry::Tier::T1_Curated;
    e.locus_id         = std::move(locus_id);
    e.haplotype_class  = std::move(haplotype_class);
    e.structural_architecture = "tandem_duplication";
    GenomicCoords c;
    c.assembly = "GRCh38";
    c.chrom    = std::string(chrom);
    c.start    = start;
    c.end      = end;
    c.strand   = '+';
    e.coords_by_assembly.emplace(c.assembly, std::move(c));
    return e;
}

SegDupCatalogEntry MakeBulk(std::string locus_id,
                             std::string_view chrom,
                             std::int64_t start,
                             std::int64_t end) {
    SegDupCatalogEntry e = MakeCurated(std::move(locus_id), /*hap=*/"",
                                        chrom, start, end);
    e.tier             = SegDupCatalogEntry::Tier::T2_Bulk;
    e.haplotype_class.clear();
    return e;
}

// Plain for_each_curated visits exactly the curated entries.
TEST(ForEachCurated, VisitsCuratedOnly_SkipsBulk) {
    SegDupCatalog cat;
    cat.add_entry(MakeCurated("CLINICAL_NPHP1",
                               "nphp1_block_paralog",
                               "chr2", 110'000'000, 110'500'000));
    cat.add_entry(MakeCurated("CLINICAL_FCGR",
                               "fcgr_cluster_cnv",
                               "chr1", 161'000'000, 161'700'000));
    cat.add_entry(MakeBulk("BULK_UCSC_42",
                            "chr3", 50'000'000, 50'010'000));
    cat.add_entry(MakeCurated("Y_PALINDROME_P3",
                               "y_palindrome_p3",
                               "chrY", 21'000'000, 21'400'000));

    std::vector<std::string> visited;
    cat.for_each_curated([&](const SegDupCatalogEntry& e) {
        visited.push_back(e.locus_id);
    });

    ASSERT_EQ(visited.size(), 3u);
    std::unordered_set<std::string> as_set(visited.begin(), visited.end());
    EXPECT_TRUE(as_set.count("CLINICAL_NPHP1"));
    EXPECT_TRUE(as_set.count("CLINICAL_FCGR"));
    EXPECT_TRUE(as_set.count("Y_PALINDROME_P3"));
    EXPECT_FALSE(as_set.count("BULK_UCSC_42"))
        << "T2 bulk record must not appear in for_each_curated";
}

// Predicate-filtered overload composes correctly.
TEST(ForEachCurated, PredicateFiltered_OnlyMatchingEntries) {
    SegDupCatalog cat;
    cat.add_entry(MakeCurated("HLA_A_block",  "hla_a_class_i",  "chr6", 29'000'000, 29'050'000));
    cat.add_entry(MakeCurated("HLA_B_block",  "hla_b_class_i",  "chr6", 31'300'000, 31'400'000));
    cat.add_entry(MakeCurated("HLA_DRB1",     "hla_drb1_class_ii", "chr6", 32'500'000, 32'600'000));

    std::vector<std::string> visited;
    cat.for_each_curated(
        [](const SegDupCatalogEntry& e) {
            // Predicate: only Class I HLA entries
            return e.haplotype_class.find("class_i") != std::string::npos
                && e.haplotype_class.find("class_ii") == std::string::npos;
        },
        [&](const SegDupCatalogEntry& e) { visited.push_back(e.locus_id); });

    ASSERT_EQ(visited.size(), 2u);
    std::sort(visited.begin(), visited.end());
    EXPECT_EQ(visited[0], "HLA_A_block");
    EXPECT_EQ(visited[1], "HLA_B_block");
}

// Empty catalog → callback never fires, no error.
TEST(ForEachCurated, EmptyCatalog_NoInvocations) {
    SegDupCatalog cat;
    std::size_t calls = 0;
    cat.for_each_curated([&](const SegDupCatalogEntry&) { ++calls; });
    EXPECT_EQ(calls, 0u);
}

// The classify-adapter resolves a non-IGHG haplotype_class — proving the
// old IGHG-only hard-code (the comment block at segdup_catalog_adapter.cpp
// referenced 'k_known_ids = {IGHG4_chimdup_tandem, IGHG_canondup_nahr_block,
// IGHG4_chimdup_canonical_arch}') is no longer in play.
TEST(ForEachCurated, AdapterFindsNonIghgHaplotypeClass) {
    SegDupCatalog cat;
    cat.add_entry(MakeCurated("NCF1_cluster",
                               "ncf1_pseudo_converter",
                               "chr7", 74'000'000, 74'300'000));
    cat.add_entry(MakeCurated("NBPF_olduvai",
                               "nbpf_olduvai_cluster",
                               "chr1", 144'000'000, 145'500'000));

    llmap::classify::igh::SegDupCatalogAdapter adapter(cat, "GRCh38");

    auto hit_a = adapter.LookupByHaplotypeClass("ncf1_pseudo_converter");
    ASSERT_TRUE(hit_a.has_value())
        << "Adapter must find NCF1 haplotype_class via for_each_curated; "
           "previously the IGHG-trio whitelist would have rejected this.";
    EXPECT_EQ(hit_a->locus_id, "NCF1_cluster");

    auto hit_b = adapter.LookupByHaplotypeClass("nbpf_olduvai_cluster");
    ASSERT_TRUE(hit_b.has_value());
    EXPECT_EQ(hit_b->locus_id, "NBPF_olduvai");

    auto miss = adapter.LookupByHaplotypeClass("does_not_exist_class_xyz");
    EXPECT_FALSE(miss.has_value());
}

}  // namespace
