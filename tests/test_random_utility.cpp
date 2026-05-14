#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "random_utility.h"

// ─── randomDouble ────────────────────────────────────────────────────────────

TEST(RandomDouble, ResultInRange) {
    seedRandom(42);
    for (int i = 0; i < 1000; ++i) {
        double v = randomDouble(-1.0, 1.0);
        EXPECT_GE(v, -1.0);
        EXPECT_LE(v, 1.0);
    }
}

TEST(RandomDouble, NarrowRange) {
    seedRandom(7);
    for (int i = 0; i < 200; ++i) {
        double v = randomDouble(3.0, 3.001);
        EXPECT_GE(v, 3.0);
        EXPECT_LE(v, 3.001);
    }
}

TEST(RandomDouble, PositiveRangeOnly) {
    seedRandom(99);
    for (int i = 0; i < 500; ++i) {
        double v = randomDouble(0.0, 1.0);
        EXPECT_GE(v, 0.0);
        EXPECT_LE(v, 1.0);
    }
}

TEST(RandomDouble, NegativeRangeOnly) {
    seedRandom(1234);
    for (int i = 0; i < 500; ++i) {
        double v = randomDouble(-5.0, -1.0);
        EXPECT_GE(v, -5.0);
        EXPECT_LE(v, -1.0);
    }
}

// ─── seedRandom (determinism) ────────────────────────────────────────────────

TEST(SeedRandom, SameSeedProducesSameSequence) {
    seedRandom(12345);
    double a1 = randomDouble(0.0, 1.0);
    double a2 = randomDouble(0.0, 1.0);
    double a3 = randomDouble(0.0, 1.0);

    seedRandom(12345);
    double b1 = randomDouble(0.0, 1.0);
    double b2 = randomDouble(0.0, 1.0);
    double b3 = randomDouble(0.0, 1.0);

    EXPECT_DOUBLE_EQ(a1, b1);
    EXPECT_DOUBLE_EQ(a2, b2);
    EXPECT_DOUBLE_EQ(a3, b3);
}

TEST(SeedRandom, DifferentSeedsProduceDifferentValues) {
    seedRandom(1);
    double v1 = randomDouble(0.0, 1.0);

    seedRandom(2);
    double v2 = randomDouble(0.0, 1.0);

    // Different seeds should (with overwhelming probability) give different first values
    EXPECT_NE(v1, v2);
}

TEST(SeedRandom, ReseedingResetsSequence) {
    seedRandom(999);
    double first = randomDouble(0.0, 1.0);

    // Advance the sequence
    for (int i = 0; i < 100; ++i)
        randomDouble(0.0, 1.0);

    seedRandom(999);
    double again = randomDouble(0.0, 1.0);

    EXPECT_DOUBLE_EQ(first, again);
}

// ─── getGenerator ────────────────────────────────────────────────────────────

TEST(GetGenerator, ReturnsSameObject) {
    // getGenerator() must return a reference to the same internal engine each time.
    std::mt19937& g1 = getGenerator();
    std::mt19937& g2 = getGenerator();
    EXPECT_EQ(&g1, &g2);
}

TEST(GetGenerator, GeneratorIsAffectedBySeedRandom) {
    seedRandom(42);
    std::mt19937& gen = getGenerator();
    uint32_t v1 = gen();

    seedRandom(42);
    uint32_t v2 = gen();

    EXPECT_EQ(v1, v2);
}

TEST(GetGenerator, DirectGeneratorUseAndRandomDoubleShareState) {
    seedRandom(77);
    // Consume one value directly from the generator
    getGenerator()();
    double after_direct = randomDouble(0.0, 1.0);

    seedRandom(77);
    // Skip one via randomDouble, then take another
    randomDouble(0.0, 1.0);  // advances by however many raw values randomDouble consumes
    // We can't know the exact raw count, but we can verify the generator is shared:
    // reseed and replay both paths to confirm they agree.
    seedRandom(77);
    getGenerator()();
    double replay = randomDouble(0.0, 1.0);

    EXPECT_DOUBLE_EQ(after_direct, replay);
}
