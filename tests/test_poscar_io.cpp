#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <string>

#include "poscar_file.h"

static const std::string kNaClPath = std::string(TEST_DATA_DIR) + "/NaCl_conv_fcc.poscar";

TEST(PoscarIO, ReadValidPOSCAR) {
    POSCAR poscar;
    ASSERT_TRUE(poscar.readPOSCAR(kNaClPath));

    EXPECT_EQ(poscar.comment, "Na4 Cl4");
    EXPECT_DOUBLE_EQ(poscar.scale, 1.0);
    EXPECT_TRUE(poscar.is_direct);
    EXPECT_EQ(poscar.total_atoms, 8);

    ASSERT_EQ(poscar.elements.size(), 2u);
    EXPECT_EQ(poscar.elements[0], "Na");
    EXPECT_EQ(poscar.elements[1], "Cl");

    ASSERT_EQ(poscar.num_atoms.size(), 2u);
    EXPECT_EQ(poscar.num_atoms[0], 4);
    EXPECT_EQ(poscar.num_atoms[1], 4);

    EXPECT_EQ(poscar.coordinates.size(), 8u);
}

TEST(PoscarIO, LatticeVectors) {
    POSCAR poscar;
    ASSERT_TRUE(poscar.readPOSCAR(kNaClPath));

    const double a = 5.5881264354399347;
    constexpr double tol = 1e-10;

    EXPECT_NEAR(poscar.lattice[0][0], a, tol);
    EXPECT_NEAR(poscar.lattice[1][1], a, tol);
    EXPECT_NEAR(poscar.lattice[2][2], a, tol);

    // Off-diagonal should be ~0
    EXPECT_NEAR(poscar.lattice[0][1], 0.0, tol);
    EXPECT_NEAR(poscar.lattice[0][2], 0.0, tol);
    EXPECT_NEAR(poscar.lattice[1][0], 0.0, tol);
    EXPECT_NEAR(poscar.lattice[1][2], 0.0, tol);
    EXPECT_NEAR(poscar.lattice[2][0], 0.0, tol);
    EXPECT_NEAR(poscar.lattice[2][1], 0.0, tol);
}

TEST(PoscarIO, CoordinateValues) {
    POSCAR poscar;
    ASSERT_TRUE(poscar.readPOSCAR(kNaClPath));

    constexpr double tol = 1e-10;

    // Na at (0, 0, 0)
    EXPECT_NEAR(poscar.coordinates[0].x, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[0].y, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[0].z, 0.0, tol);

    // Na at (0, 0.5, 0.5)
    EXPECT_NEAR(poscar.coordinates[1].x, 0.0, tol);
    EXPECT_NEAR(poscar.coordinates[1].y, 0.5, tol);
    EXPECT_NEAR(poscar.coordinates[1].z, 0.5, tol);

    // Cl at (0.5, 0.5, 0.5)
    EXPECT_NEAR(poscar.coordinates[7].x, 0.5, tol);
    EXPECT_NEAR(poscar.coordinates[7].y, 0.5, tol);
    EXPECT_NEAR(poscar.coordinates[7].z, 0.5, tol);
}

TEST(PoscarIO, ReadNonExistentFile) {
    POSCAR poscar;
    EXPECT_FALSE(poscar.readPOSCAR("nonexistent_file.poscar"));
}

TEST(PoscarIO, ReadWriteRoundTrip) {
    POSCAR original;
    ASSERT_TRUE(original.readPOSCAR(kNaClPath));

    const std::string tmpFile = std::string(TEST_DATA_DIR) + "/test_roundtrip_tmp.poscar";
    ASSERT_TRUE(original.writePOSCAR(tmpFile));

    POSCAR reloaded;
    ASSERT_TRUE(reloaded.readPOSCAR(tmpFile));
    std::remove(tmpFile.c_str());

    EXPECT_EQ(original.comment, reloaded.comment);
    EXPECT_DOUBLE_EQ(original.scale, reloaded.scale);
    EXPECT_EQ(original.is_direct, reloaded.is_direct);
    EXPECT_EQ(original.total_atoms, reloaded.total_atoms);
    EXPECT_EQ(original.elements, reloaded.elements);
    EXPECT_EQ(original.num_atoms, reloaded.num_atoms);

    constexpr double tol = 1e-9;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(original.lattice[i][j], reloaded.lattice[i][j], tol);
        }
    }

    ASSERT_EQ(original.coordinates.size(), reloaded.coordinates.size());
    for (size_t i = 0; i < original.coordinates.size(); ++i) {
        EXPECT_NEAR(original.coordinates[i].x, reloaded.coordinates[i].x, tol);
        EXPECT_NEAR(original.coordinates[i].y, reloaded.coordinates[i].y, tol);
        EXPECT_NEAR(original.coordinates[i].z, reloaded.coordinates[i].z, tol);
    }
}
