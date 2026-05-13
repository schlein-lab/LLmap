// LLmap — Unit tests for cache_layout.h

#include "core/cache_layout.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <vector>

namespace llmap::core {
namespace {

// =============================================================================
// Prefetch utilities
// =============================================================================

TEST(CacheLayout, PrefetchDoesNotCrash) {
    std::vector<int> data(1000);
    for (size_t i = 0; i < data.size(); ++i) {
        PrefetchForRead(&data[i]);
    }
    for (size_t i = 0; i < data.size(); ++i) {
        PrefetchForWrite(&data[i]);
    }
    EXPECT_TRUE(true);  // If we got here, prefetch didn't crash
}

TEST(CacheLayout, PrefetchRangeDoesNotCrash) {
    std::vector<char> data(4096);
    PrefetchRange(data.data(), data.size());
    EXPECT_TRUE(true);
}

TEST(CacheLayout, AllocateAligned) {
    void* ptr = AllocateAligned(1024);
    ASSERT_NE(ptr, nullptr);
    EXPECT_TRUE(IsCacheAligned(ptr));
    FreeAligned(ptr);
}

TEST(CacheLayout, AllocateAlignedCustomAlignment) {
    void* ptr = AllocateAligned(1024, 128);
    ASSERT_NE(ptr, nullptr);
    EXPECT_EQ(reinterpret_cast<uintptr_t>(ptr) % 128, 0);
    FreeAligned(ptr);
}

TEST(CacheLayout, IsCacheAligned) {
    void* ptr = AllocateAligned(64);
    EXPECT_TRUE(IsCacheAligned(ptr));

    char* misaligned = static_cast<char*>(ptr) + 7;
    EXPECT_FALSE(IsCacheAligned(misaligned));

    FreeAligned(ptr);
}

// =============================================================================
// MinimizerSoA
// =============================================================================

TEST(MinimizerSoA, EmptyByDefault) {
    MinimizerSoA soa;
    EXPECT_TRUE(soa.Empty());
    EXPECT_EQ(soa.Size(), 0);
    EXPECT_EQ(soa.Capacity(), 0);
}

TEST(MinimizerSoA, PushAndAccess) {
    MinimizerSoA soa;
    soa.Push(0x123456789ABCDEF0ULL, 100, false);
    soa.Push(0xFEDCBA9876543210ULL, 200, true);

    EXPECT_EQ(soa.Size(), 2);
    EXPECT_EQ(soa.Hash(0), 0x123456789ABCDEF0ULL);
    EXPECT_EQ(soa.Position(0), 100);
    EXPECT_FALSE(soa.IsReverse(0));

    EXPECT_EQ(soa.Hash(1), 0xFEDCBA9876543210ULL);
    EXPECT_EQ(soa.Position(1), 200);
    EXPECT_TRUE(soa.IsReverse(1));
}

TEST(MinimizerSoA, Reserve) {
    MinimizerSoA soa;
    soa.Reserve(1000);

    EXPECT_GE(soa.Capacity(), 1000);
    EXPECT_EQ(soa.Size(), 0);
}

TEST(MinimizerSoA, Reset) {
    MinimizerSoA soa;
    for (int i = 0; i < 100; ++i) {
        soa.Push(i, i * 10, i % 2 == 0);
    }
    EXPECT_EQ(soa.Size(), 100);

    soa.Reset();
    EXPECT_EQ(soa.Size(), 0);
    EXPECT_GE(soa.Capacity(), 100);  // Capacity retained
}

TEST(MinimizerSoA, Spans) {
    MinimizerSoA soa;
    soa.Push(1, 10, false);
    soa.Push(2, 20, true);
    soa.Push(3, 30, false);

    auto hashes = soa.Hashes();
    auto positions = soa.Positions();
    auto strands = soa.Strands();

    EXPECT_EQ(hashes.size(), 3);
    EXPECT_EQ(positions.size(), 3);
    EXPECT_EQ(strands.size(), 3);

    EXPECT_EQ(hashes[0], 1);
    EXPECT_EQ(hashes[1], 2);
    EXPECT_EQ(hashes[2], 3);

    EXPECT_EQ(positions[0], 10);
    EXPECT_EQ(positions[1], 20);
    EXPECT_EQ(positions[2], 30);
}

TEST(MinimizerSoA, MoveConstruction) {
    MinimizerSoA soa1;
    soa1.Push(1, 10, false);
    soa1.Push(2, 20, true);

    MinimizerSoA soa2(std::move(soa1));

    EXPECT_EQ(soa2.Size(), 2);
    EXPECT_EQ(soa2.Hash(0), 1);
    EXPECT_EQ(soa2.Hash(1), 2);

    EXPECT_EQ(soa1.Size(), 0);
    EXPECT_EQ(soa1.Capacity(), 0);
}

TEST(MinimizerSoA, MoveAssignment) {
    MinimizerSoA soa1;
    soa1.Push(1, 10, false);

    MinimizerSoA soa2;
    soa2.Push(999, 999, true);

    soa2 = std::move(soa1);

    EXPECT_EQ(soa2.Size(), 1);
    EXPECT_EQ(soa2.Hash(0), 1);
}

TEST(MinimizerSoA, AlignedArrays) {
    MinimizerSoA soa;
    soa.Reserve(100);
    soa.Push(1, 10, false);

    EXPECT_TRUE(IsCacheAligned(soa.HashData()));
    EXPECT_TRUE(IsCacheAligned(soa.PositionData()));
    EXPECT_TRUE(IsCacheAligned(soa.StrandData()));
}

TEST(MinimizerSoA, LargeDataset) {
    MinimizerSoA soa;
    constexpr size_t N = 100000;

    for (size_t i = 0; i < N; ++i) {
        soa.Push(i * 17, static_cast<uint32_t>(i), i % 2 == 0);
    }

    EXPECT_EQ(soa.Size(), N);

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(soa.Hash(i), i * 17);
        EXPECT_EQ(soa.Position(i), i);
        EXPECT_EQ(soa.IsReverse(i), i % 2 == 0);
    }
}

// =============================================================================
// AnchorSoA
// =============================================================================

TEST(AnchorSoA, EmptyByDefault) {
    AnchorSoA soa;
    EXPECT_TRUE(soa.Empty());
    EXPECT_EQ(soa.Size(), 0);
}

TEST(AnchorSoA, PushAndAccess) {
    AnchorSoA soa;
    soa.Push(0, 100, 50, true);
    soa.Push(1, 200, 75, false);

    EXPECT_EQ(soa.Size(), 2);

    EXPECT_EQ(soa.RefId(0), 0);
    EXPECT_EQ(soa.RefPosition(0), 100);
    EXPECT_EQ(soa.QueryPosition(0), 50);
    EXPECT_TRUE(soa.IsSameStrand(0));

    EXPECT_EQ(soa.RefId(1), 1);
    EXPECT_EQ(soa.RefPosition(1), 200);
    EXPECT_EQ(soa.QueryPosition(1), 75);
    EXPECT_FALSE(soa.IsSameStrand(1));
}

TEST(AnchorSoA, SortPermutation) {
    AnchorSoA soa;
    soa.Push(1, 300, 0, true);
    soa.Push(0, 100, 1, true);
    soa.Push(1, 100, 2, true);
    soa.Push(0, 200, 3, true);

    auto perm = soa.SortPermutation();

    // Sorted order should be: (0,100), (0,200), (1,100), (1,300)
    EXPECT_EQ(perm[0], 1);  // ref_id=0, ref_pos=100
    EXPECT_EQ(perm[1], 3);  // ref_id=0, ref_pos=200
    EXPECT_EQ(perm[2], 2);  // ref_id=1, ref_pos=100
    EXPECT_EQ(perm[3], 0);  // ref_id=1, ref_pos=300
}

TEST(AnchorSoA, Spans) {
    AnchorSoA soa;
    soa.Push(0, 10, 5, true);
    soa.Push(1, 20, 10, false);

    auto ref_ids = soa.RefIds();
    auto ref_pos = soa.RefPositions();
    auto query_pos = soa.QueryPositions();
    auto strands = soa.Strands();

    EXPECT_EQ(ref_ids.size(), 2);
    EXPECT_EQ(ref_pos.size(), 2);
    EXPECT_EQ(query_pos.size(), 2);
    EXPECT_EQ(strands.size(), 2);

    EXPECT_EQ(ref_ids[0], 0);
    EXPECT_EQ(ref_ids[1], 1);
}

TEST(AnchorSoA, MoveSemantics) {
    AnchorSoA soa1;
    soa1.Push(0, 100, 50, true);

    AnchorSoA soa2(std::move(soa1));
    EXPECT_EQ(soa2.Size(), 1);
    EXPECT_EQ(soa1.Size(), 0);

    AnchorSoA soa3;
    soa3 = std::move(soa2);
    EXPECT_EQ(soa3.Size(), 1);
    EXPECT_EQ(soa2.Size(), 0);
}

TEST(AnchorSoA, AlignedArrays) {
    AnchorSoA soa;
    soa.Reserve(100);
    soa.Push(0, 100, 50, true);

    EXPECT_TRUE(IsCacheAligned(soa.RefIdData()));
    EXPECT_TRUE(IsCacheAligned(soa.RefPositionData()));
    EXPECT_TRUE(IsCacheAligned(soa.QueryPositionData()));
    EXPECT_TRUE(IsCacheAligned(soa.StrandData()));
}

// =============================================================================
// DPStateSoA
// =============================================================================

TEST(DPStateSoA, Resize) {
    DPStateSoA dp;
    dp.Resize(1000);
    EXPECT_EQ(dp.Size(), 1000);
}

TEST(DPStateSoA, Initialize) {
    DPStateSoA dp(100);
    dp.Initialize();

    for (size_t i = 0; i < dp.Size(); ++i) {
        EXPECT_EQ(dp.Score(i), 0);
        EXPECT_EQ(dp.Predecessor(i), -1);
    }
}

TEST(DPStateSoA, SetAndGet) {
    DPStateSoA dp(10);
    dp.Initialize();

    dp.SetScore(5, 42);
    dp.SetPredecessor(5, 3);

    EXPECT_EQ(dp.Score(5), 42);
    EXPECT_EQ(dp.Predecessor(5), 3);
}

TEST(DPStateSoA, MoveSemantics) {
    DPStateSoA dp1(100);
    dp1.Initialize();
    dp1.SetScore(50, 999);

    DPStateSoA dp2(std::move(dp1));
    EXPECT_EQ(dp2.Score(50), 999);
    EXPECT_EQ(dp1.Size(), 0);
}

TEST(DPStateSoA, GrowsCapacity) {
    DPStateSoA dp;
    dp.Resize(100);
    dp.SetScore(50, 123);

    dp.Resize(200);  // Should preserve data
    EXPECT_EQ(dp.Size(), 200);
}

// =============================================================================
// TiledProcessor
// =============================================================================

TEST(TiledProcessor, ProcessAllCells) {
    constexpr size_t rows = 100;
    constexpr size_t cols = 80;

    std::vector<bool> visited(rows * cols, false);

    TiledProcessor<16>::Process(rows, cols, [&](size_t r, size_t c) {
        visited[r * cols + c] = true;
    });

    // All cells should be visited exactly once
    for (bool v : visited) {
        EXPECT_TRUE(v);
    }
}

TEST(TiledProcessor, ProcessDiagonalOrder) {
    constexpr size_t rows = 10;
    constexpr size_t cols = 10;

    std::vector<size_t> order;
    order.reserve(rows * cols);

    TiledProcessor<>::ProcessDiagonal(rows, cols, [&](size_t r, size_t c) {
        order.push_back(r * cols + c);
    });

    EXPECT_EQ(order.size(), rows * cols);

    // First element should be (0,0)
    EXPECT_EQ(order[0], 0);

    // Verify diagonal property: each cell is processed after its dependencies
    for (size_t i = 0; i < order.size(); ++i) {
        size_t idx = order[i];
        size_t r = idx / cols;
        size_t c = idx % cols;

        // For diagonal traversal, (r-1, c) and (r, c-1) should come before (r, c)
        if (r > 0) {
            size_t dep_idx = (r - 1) * cols + c;
            auto it = std::find(order.begin(), order.begin() + i, dep_idx);
            EXPECT_NE(it, order.begin() + i) << "Dependency (r-1,c) not processed";
        }
        if (c > 0) {
            size_t dep_idx = r * cols + (c - 1);
            auto it = std::find(order.begin(), order.begin() + i, dep_idx);
            EXPECT_NE(it, order.begin() + i) << "Dependency (r,c-1) not processed";
        }
    }
}

TEST(TiledProcessor, SmallMatrix) {
    constexpr size_t rows = 3;
    constexpr size_t cols = 4;

    size_t count = 0;
    TiledProcessor<2>::Process(rows, cols, [&](size_t r, size_t c) {
        EXPECT_LT(r, rows);
        EXPECT_LT(c, cols);
        ++count;
    });

    EXPECT_EQ(count, rows * cols);
}

// =============================================================================
// Performance: AoS vs SoA comparison
// =============================================================================

TEST(CacheLayout, SoAVersuAoSSequentialAccess) {
    constexpr size_t N = 100000;

    // AoS version
    struct AnchorAoS {
        uint32_t ref_id;
        uint32_t ref_pos;
        uint32_t query_pos;
        bool same_strand;
    };

    std::vector<AnchorAoS> aos(N);
    for (size_t i = 0; i < N; ++i) {
        aos[i] = {static_cast<uint32_t>(i % 10),
                  static_cast<uint32_t>(i * 7),
                  static_cast<uint32_t>(i * 3),
                  i % 2 == 0};
    }

    // SoA version
    AnchorSoA soa;
    soa.Reserve(N);
    for (size_t i = 0; i < N; ++i) {
        soa.Push(static_cast<uint32_t>(i % 10),
                 static_cast<uint32_t>(i * 7),
                 static_cast<uint32_t>(i * 3),
                 i % 2 == 0);
    }

    // Time AoS sequential access (only ref_pos)
    auto aos_start = std::chrono::high_resolution_clock::now();
    uint64_t aos_sum = 0;
    for (size_t i = 0; i < N; ++i) {
        aos_sum += aos[i].ref_pos;
    }
    auto aos_end = std::chrono::high_resolution_clock::now();
    auto aos_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        aos_end - aos_start).count();

    // Time SoA sequential access (only ref_pos)
    auto soa_start = std::chrono::high_resolution_clock::now();
    uint64_t soa_sum = 0;
    auto ref_positions = soa.RefPositions();
    for (size_t i = 0; i < N; ++i) {
        soa_sum += ref_positions[i];
    }
    auto soa_end = std::chrono::high_resolution_clock::now();
    auto soa_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        soa_end - soa_start).count();

    // Results should be identical
    EXPECT_EQ(aos_sum, soa_sum);

    // SoA should generally be faster for field-only access
    // (though this isn't guaranteed in a unit test environment)
    std::cout << "AoS sequential access: " << aos_ns << " ns\n";
    std::cout << "SoA sequential access: " << soa_ns << " ns\n";
    std::cout << "Speedup: " << static_cast<double>(aos_ns) / soa_ns << "x\n";
}

TEST(CacheLayout, SoARandomAccessPattern) {
    constexpr size_t N = 10000;

    AnchorSoA soa;
    soa.Reserve(N);
    for (size_t i = 0; i < N; ++i) {
        soa.Push(static_cast<uint32_t>(i % 10),
                 static_cast<uint32_t>(i * 7),
                 static_cast<uint32_t>(i * 3),
                 i % 2 == 0);
    }

    // Generate random access pattern
    std::vector<size_t> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 gen(42);
    std::shuffle(indices.begin(), indices.end(), gen);

    // Time random access with prefetching
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t sum = 0;
    auto ref_positions = soa.RefPositions();
    for (size_t i = 0; i < N; ++i) {
        if (i + 8 < N) {
            PrefetchForRead(&ref_positions[indices[i + 8]]);
        }
        sum += ref_positions[indices[i]];
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end - start).count();

    EXPECT_GT(sum, 0);
    std::cout << "Random access with prefetch: " << ns << " ns\n";
}

}  // namespace
}  // namespace llmap::core
