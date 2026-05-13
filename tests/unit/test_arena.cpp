#include "core/arena.h"

#include <gtest/gtest.h>

#include <thread>
#include <vector>

namespace llmap::core {
namespace {

// =============================================================================
// Arena Tests
// =============================================================================

TEST(ArenaTest, BasicAllocation) {
    Arena arena(1024);

    int* p1 = arena.Allocate<int>();
    ASSERT_NE(p1, nullptr);
    *p1 = 42;
    EXPECT_EQ(*p1, 42);

    double* p2 = arena.Allocate<double>();
    ASSERT_NE(p2, nullptr);
    *p2 = 3.14;
    EXPECT_EQ(*p2, 3.14);
}

TEST(ArenaTest, ArrayAllocation) {
    Arena arena;

    int* arr = arena.Allocate<int>(100);
    ASSERT_NE(arr, nullptr);

    for (int i = 0; i < 100; ++i) {
        arr[i] = i * 2;
    }

    for (int i = 0; i < 100; ++i) {
        EXPECT_EQ(arr[i], i * 2);
    }
}

TEST(ArenaTest, SpanAllocation) {
    Arena arena;

    auto span = arena.AllocateSpan<float>(50);
    EXPECT_EQ(span.size(), 50);

    for (size_t i = 0; i < span.size(); ++i) {
        span[i] = static_cast<float>(i);
    }

    for (size_t i = 0; i < span.size(); ++i) {
        EXPECT_EQ(span[i], static_cast<float>(i));
    }
}

TEST(ArenaTest, Construct) {
    Arena arena;

    struct Point {
        int x, y;
        Point(int x_, int y_) : x(x_), y(y_) {}
    };

    Point* p = arena.Construct<Point>(10, 20);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->x, 10);
    EXPECT_EQ(p->y, 20);
}

TEST(ArenaTest, ResetKeepsFirstBlock) {
    Arena arena(1024);

    arena.Allocate(512);
    EXPECT_EQ(arena.BlocksUsed(), 1);

    arena.Reset();
    EXPECT_EQ(arena.BlocksUsed(), 1);  // First block retained
    EXPECT_EQ(arena.BytesAllocated(), 0);

    int* p = arena.Allocate<int>();
    ASSERT_NE(p, nullptr);
    *p = 123;
    EXPECT_EQ(*p, 123);
}

TEST(ArenaTest, ResetReleasesExtraBlocks) {
    Arena arena(512);

    arena.Allocate(400);
    arena.Allocate(400);  // Forces second block
    EXPECT_GE(arena.BlocksUsed(), 2);

    arena.Reset();
    EXPECT_EQ(arena.BlocksUsed(), 1);
}

TEST(ArenaTest, Clear) {
    Arena arena(1024);

    arena.Allocate(512);
    EXPECT_EQ(arena.BlocksUsed(), 1);

    arena.Clear();
    EXPECT_EQ(arena.BlocksUsed(), 0);
    EXPECT_EQ(arena.BytesAllocated(), 0);
}

TEST(ArenaTest, LargeAllocation) {
    Arena arena(1024);

    char* large = arena.Allocate<char>(10000);
    ASSERT_NE(large, nullptr);
    std::memset(large, 'A', 10000);
    EXPECT_EQ(large[9999], 'A');
}

TEST(ArenaTest, Alignment) {
    Arena arena;

    char* c = arena.Allocate<char>();
    double* d = arena.Allocate<double>();

    EXPECT_EQ(reinterpret_cast<uintptr_t>(d) % alignof(double), 0);
}

TEST(ArenaTest, Statistics) {
    Arena arena(1024);

    EXPECT_EQ(arena.BytesAllocated(), 0);
    EXPECT_EQ(arena.BlocksUsed(), 0);

    arena.Allocate(100);
    EXPECT_GE(arena.BytesAllocated(), 100);
    EXPECT_EQ(arena.BlocksUsed(), 1);
}

// =============================================================================
// ScratchBuffer Tests
// =============================================================================

TEST(ScratchBufferTest, BasicOperations) {
    ScratchBuffer<int> buf;

    EXPECT_TRUE(buf.empty());
    EXPECT_EQ(buf.size(), 0);

    buf.push_back(1);
    buf.push_back(2);
    buf.push_back(3);

    EXPECT_EQ(buf.size(), 3);
    EXPECT_EQ(buf[0], 1);
    EXPECT_EQ(buf[1], 2);
    EXPECT_EQ(buf[2], 3);
}

TEST(ScratchBufferTest, FrontBack) {
    ScratchBuffer<int> buf;
    buf.push_back(10);
    buf.push_back(20);
    buf.push_back(30);

    EXPECT_EQ(buf.front(), 10);
    EXPECT_EQ(buf.back(), 30);
}

TEST(ScratchBufferTest, EmplaceBack) {
    struct Pair {
        int a, b;
        Pair(int x, int y) : a(x), b(y) {}
    };

    ScratchBuffer<Pair> buf;
    auto& p = buf.emplace_back(5, 10);

    EXPECT_EQ(p.a, 5);
    EXPECT_EQ(p.b, 10);
    EXPECT_EQ(buf.size(), 1);
}

TEST(ScratchBufferTest, ClearKeepsCapacity) {
    ScratchBuffer<int> buf;
    buf.reserve(100);
    size_t cap = buf.capacity();

    buf.push_back(1);
    buf.push_back(2);
    buf.clear();

    EXPECT_EQ(buf.size(), 0);
    EXPECT_EQ(buf.capacity(), cap);
}

TEST(ScratchBufferTest, Resize) {
    ScratchBuffer<int> buf;
    buf.resize(10);

    EXPECT_EQ(buf.size(), 10);

    for (size_t i = 0; i < 10; ++i) {
        buf[i] = static_cast<int>(i);
    }

    buf.resize(5);
    EXPECT_EQ(buf.size(), 5);
}

TEST(ScratchBufferTest, ResizeZero) {
    ScratchBuffer<int> buf;
    buf.resize_zero(10);

    EXPECT_EQ(buf.size(), 10);
    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(buf[i], 0);
    }
}

TEST(ScratchBufferTest, PopBack) {
    ScratchBuffer<int> buf;
    buf.push_back(1);
    buf.push_back(2);
    buf.push_back(3);

    buf.pop_back();
    EXPECT_EQ(buf.size(), 2);
    EXPECT_EQ(buf.back(), 2);
}

TEST(ScratchBufferTest, Iteration) {
    ScratchBuffer<int> buf;
    buf.push_back(1);
    buf.push_back(2);
    buf.push_back(3);

    int sum = 0;
    for (int x : buf) {
        sum += x;
    }
    EXPECT_EQ(sum, 6);
}

TEST(ScratchBufferTest, Span) {
    ScratchBuffer<int> buf;
    buf.push_back(10);
    buf.push_back(20);

    auto s = buf.span();
    EXPECT_EQ(s.size(), 2);
    EXPECT_EQ(s[0], 10);
    EXPECT_EQ(s[1], 20);
}

TEST(ScratchBufferTest, Release) {
    ScratchBuffer<int> buf;
    buf.push_back(1);
    buf.push_back(2);
    buf.push_back(3);

    auto vec = buf.release();
    EXPECT_EQ(vec.size(), 3);
    EXPECT_EQ(vec[0], 1);
    EXPECT_EQ(vec[1], 2);
    EXPECT_EQ(vec[2], 3);

    EXPECT_EQ(buf.size(), 0);  // Buffer cleared after release
}

TEST(ScratchBufferTest, MoveConstruct) {
    ScratchBuffer<int> buf1;
    buf1.push_back(1);
    buf1.push_back(2);

    ScratchBuffer<int> buf2(std::move(buf1));
    EXPECT_EQ(buf2.size(), 2);
    EXPECT_EQ(buf2[0], 1);
    EXPECT_EQ(buf2[1], 2);

    EXPECT_EQ(buf1.size(), 0);
    EXPECT_EQ(buf1.capacity(), 0);
}

TEST(ScratchBufferTest, MoveAssign) {
    ScratchBuffer<int> buf1;
    buf1.push_back(1);

    ScratchBuffer<int> buf2;
    buf2.push_back(100);
    buf2 = std::move(buf1);

    EXPECT_EQ(buf2.size(), 1);
    EXPECT_EQ(buf2[0], 1);
}

TEST(ScratchBufferTest, LargeGrowth) {
    ScratchBuffer<int> buf;
    for (int i = 0; i < 10000; ++i) {
        buf.push_back(i);
    }

    EXPECT_EQ(buf.size(), 10000);
    for (int i = 0; i < 10000; ++i) {
        EXPECT_EQ(buf[i], i);
    }
}

// =============================================================================
// ScratchSpace Tests
// =============================================================================

TEST(ScratchSpaceTest, LocalInstance) {
    auto& s1 = ScratchSpace::Local();
    auto& s2 = ScratchSpace::Local();
    EXPECT_EQ(&s1, &s2);  // Same thread, same instance
}

TEST(ScratchSpaceTest, ArenaAccess) {
    auto& space = ScratchSpace::Local();
    auto& arena = space.GetArena();

    int* p = arena.Allocate<int>();
    ASSERT_NE(p, nullptr);
    *p = 42;
    EXPECT_EQ(*p, 42);

    space.Reset();
}

TEST(ScratchSpaceTest, BufferAccess) {
    auto& space = ScratchSpace::Local();
    auto& buf = space.GetBuffer<int>();

    buf.push_back(1);
    buf.push_back(2);
    EXPECT_EQ(buf.size(), 2);

    auto& buf2 = space.GetBuffer<int>();
    EXPECT_EQ(buf2.size(), 0);  // Cleared on each GetBuffer call
}

TEST(ScratchSpaceTest, DifferentTypes) {
    auto& space = ScratchSpace::Local();

    auto& int_buf = space.GetBuffer<int>();
    int_buf.push_back(42);

    auto& float_buf = space.GetBuffer<float>();
    float_buf.push_back(3.14f);

    EXPECT_EQ(int_buf.size(), 1);  // Different type, different buffer
}

TEST(ScratchSpaceTest, ThreadIsolation) {
    auto& main_space = ScratchSpace::Local();
    ScratchSpace* thread_space = nullptr;

    std::thread t([&] {
        thread_space = &ScratchSpace::Local();
    });
    t.join();

    EXPECT_NE(&main_space, thread_space);  // Different threads, different instances
}

// =============================================================================
// ScratchGuard Tests
// =============================================================================

TEST(ScratchGuardTest, ResetOnExit) {
    auto& space = ScratchSpace::Local();

    {
        ScratchGuard guard(space);
        auto& arena = space.GetArena();
        arena.Allocate(1000);
        EXPECT_GT(arena.BytesAllocated(), 0);
    }

    EXPECT_EQ(space.GetArena().BytesAllocated(), 0);
}

// =============================================================================
// Performance-related Tests
// =============================================================================

TEST(ArenaTest, RepeatedAllocReset) {
    Arena arena(4096);

    for (int iter = 0; iter < 100; ++iter) {
        for (int i = 0; i < 50; ++i) {
            int* p = arena.Allocate<int>(10);
            ASSERT_NE(p, nullptr);
        }
        arena.Reset();
    }

    EXPECT_EQ(arena.BlocksUsed(), 1);  // Only one block after all resets
}

TEST(ScratchBufferTest, ReuseThroughClear) {
    ScratchBuffer<int> buf;

    for (int iter = 0; iter < 100; ++iter) {
        for (int i = 0; i < 1000; ++i) {
            buf.push_back(i);
        }
        size_t cap = buf.capacity();
        buf.clear();
        EXPECT_EQ(buf.capacity(), cap);  // Capacity preserved
    }

    EXPECT_GE(buf.capacity(), 1000);  // Should have grown to accommodate 1000
}

}  // namespace
}  // namespace llmap::core
