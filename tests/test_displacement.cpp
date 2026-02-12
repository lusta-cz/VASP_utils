#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "poscar_file.h"
#include "random_utility.h"

static const std::string kNaClPath = std::string(TEST_DATA_DIR) + "/NaCl_conv_fcc.poscar";

class DisplacementTest : public ::testing::Test {
protected:
    POSCAR poscar;
    void SetUp() override {
        ASSERT_TRUE(poscar.readPOSCAR(kNaClPath));
    }
};

TEST_F(DisplacementTest, DisplacementMagnitude) {
    seedRandom(42);

    POSCAR original = poscar;
    poscar.toCartesian();
    original.toCartesian();

    poscar.displaceAtoms(1, 0.05);

    // Find which atom was displaced and check its magnitude
    const double amp = 0.05;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        double dx = poscar.coordinates[i].x - original.coordinates[i].x;
        double dy = poscar.coordinates[i].y - original.coordinates[i].y;
        double dz = poscar.coordinates[i].z - original.coordinates[i].z;
        double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        EXPECT_LE(dist, amp + 1e-12) << "Atom " << i << " displaced beyond amplitude";
    }
}

TEST_F(DisplacementTest, AllAtomsDisplaced) {
    seedRandom(123);

    auto original_coords = poscar.coordinates;
    poscar.displaceAtoms(poscar.total_atoms, 0.01);

    int changed = 0;
    constexpr double tol = 1e-15;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        if (std::abs(poscar.coordinates[i].x - original_coords[i].x) > tol ||
            std::abs(poscar.coordinates[i].y - original_coords[i].y) > tol ||
            std::abs(poscar.coordinates[i].z - original_coords[i].z) > tol) {
            ++changed;
        }
    }
    EXPECT_EQ(changed, poscar.total_atoms);
}

TEST_F(DisplacementTest, CoordinateSystemPreserved) {
    ASSERT_TRUE(poscar.is_direct);
    seedRandom(99);

    poscar.displaceAtoms(2, 0.01);
    EXPECT_TRUE(poscar.is_direct) << "Coordinate system should remain Direct after displacement";
}

TEST_F(DisplacementTest, DeterministicWithSeed) {
    seedRandom(777);
    POSCAR run1 = poscar;
    run1.displaceAtoms(4, 0.02);

    seedRandom(777);
    POSCAR run2 = poscar;
    run2.displaceAtoms(4, 0.02);

    for (size_t i = 0; i < run1.coordinates.size(); ++i) {
        EXPECT_DOUBLE_EQ(run1.coordinates[i].x, run2.coordinates[i].x) << "Atom " << i;
        EXPECT_DOUBLE_EQ(run1.coordinates[i].y, run2.coordinates[i].y) << "Atom " << i;
        EXPECT_DOUBLE_EQ(run1.coordinates[i].z, run2.coordinates[i].z) << "Atom " << i;
    }
}

TEST_F(DisplacementTest, ZeroAmplitude) {
    seedRandom(42);
    auto original_coords = poscar.coordinates;

    poscar.displaceAtoms(poscar.total_atoms, 0.0);

    constexpr double tol = 1e-15;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        EXPECT_NEAR(poscar.coordinates[i].x, original_coords[i].x, tol);
        EXPECT_NEAR(poscar.coordinates[i].y, original_coords[i].y, tol);
        EXPECT_NEAR(poscar.coordinates[i].z, original_coords[i].z, tol);
    }
}
