#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "poscar_file.h"

static const std::string kNaClPath = std::string(TEST_DATA_DIR) + "/NaCl_conv_fcc.poscar";

class CoordinateConversionTest : public ::testing::Test {
protected:
    POSCAR poscar;
    void SetUp() override {
        ASSERT_TRUE(poscar.readPOSCAR(kNaClPath));
    }
};

TEST_F(CoordinateConversionTest, DirectToCartesian) {
    ASSERT_TRUE(poscar.is_direct);
    poscar.toCartesian();
    EXPECT_FALSE(poscar.is_direct);

    const double a = 5.5881264354399347;
    const double half_a = a / 2.0;
    constexpr double tol = 1e-8;

    // Atom 0: Na at (0,0,0) direct → (0,0,0) Cartesian
    EXPECT_NEAR(poscar.coordinates[0].x, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[0].y, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[0].z, 0.0, tol);

    // Atom 1: Na at (0,0.5,0.5) direct → (0, a/2, a/2) Cartesian
    EXPECT_NEAR(poscar.coordinates[1].x, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[1].y, half_a, tol);
    EXPECT_NEAR(poscar.coordinates[1].z, half_a, tol);

    // Atom 3: Na at (0.5,0.5,0) direct → (a/2, a/2, 0) Cartesian
    EXPECT_NEAR(poscar.coordinates[3].x, half_a, tol);
    EXPECT_NEAR(poscar.coordinates[3].y, half_a, tol);
    EXPECT_NEAR(poscar.coordinates[3].z, 0.0, tol);
}

TEST_F(CoordinateConversionTest, CartesianToDirect) {
    // Save original Direct coordinates
    auto original_coords = poscar.coordinates;

    poscar.toCartesian();
    ASSERT_FALSE(poscar.is_direct);

    poscar.toDirect();
    ASSERT_TRUE(poscar.is_direct);

    constexpr double tol = 1e-9;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        EXPECT_NEAR(poscar.coordinates[i].x, original_coords[i].x, tol) << "Mismatch at atom " << i << " x";
        EXPECT_NEAR(poscar.coordinates[i].y, original_coords[i].y, tol) << "Mismatch at atom " << i << " y";
        EXPECT_NEAR(poscar.coordinates[i].z, original_coords[i].z, tol) << "Mismatch at atom " << i << " z";
    }
}

TEST_F(CoordinateConversionTest, RoundTripStability) {
    auto original_coords = poscar.coordinates;

    for (int cycle = 0; cycle < 5; ++cycle) {
        poscar.toCartesian();
        poscar.toDirect();
    }

    constexpr double tol = 1e-8;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        EXPECT_NEAR(poscar.coordinates[i].x, original_coords[i].x, tol)
            << "Diverged at atom " << i << " after 5 round trips";
        EXPECT_NEAR(poscar.coordinates[i].y, original_coords[i].y, tol);
        EXPECT_NEAR(poscar.coordinates[i].z, original_coords[i].z, tol);
    }
}

TEST_F(CoordinateConversionTest, FlagUpdates) {
    EXPECT_TRUE(poscar.is_direct);

    poscar.toCartesian();
    EXPECT_FALSE(poscar.is_direct);

    poscar.toDirect();
    EXPECT_TRUE(poscar.is_direct);
}

TEST_F(CoordinateConversionTest, AlreadyDirectIsNoop) {
    ASSERT_TRUE(poscar.is_direct);
    auto before = poscar.coordinates;

    poscar.toDirect();

    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].x, before[i].x);
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].y, before[i].y);
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].z, before[i].z);
    }
}

TEST_F(CoordinateConversionTest, AlreadyCartesianIsNoop) {
    poscar.toCartesian();
    ASSERT_FALSE(poscar.is_direct);
    auto before = poscar.coordinates;

    poscar.toCartesian();

    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].x, before[i].x);
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].y, before[i].y);
        EXPECT_DOUBLE_EQ(poscar.coordinates[i].z, before[i].z);
    }
}
