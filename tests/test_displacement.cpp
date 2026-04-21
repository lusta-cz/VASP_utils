#include <gtest/gtest.h>

#include <cmath>
#include <set>
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

    // Find which atom was displaced and check component-wise bounds
    const double amp = 0.05;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        double dx = poscar.coordinates[i].x - original.coordinates[i].x;
        double dy = poscar.coordinates[i].y - original.coordinates[i].y;
        double dz = poscar.coordinates[i].z - original.coordinates[i].z;
        EXPECT_LE(std::abs(dx), amp + 1e-12) << "Atom " << i << " dx beyond amplitude";
        EXPECT_LE(std::abs(dy), amp + 1e-12) << "Atom " << i << " dy beyond amplitude";
        EXPECT_LE(std::abs(dz), amp + 1e-12) << "Atom " << i << " dz beyond amplitude";
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

TEST_F(DisplacementTest, SpeciesRestrictedDisplacement) {
    seedRandom(12);
    POSCAR original = poscar;

    std::vector<size_t> cl_indices;
    const auto species = poscar.atomSpecies();
    for (size_t i = 0; i < species.size(); ++i) {
        if (species[i] == "Cl") {
            cl_indices.push_back(i);
        }
    }

    DisplacementOptions options;
    options.candidate_indices = cl_indices;
    options.ampx = 0.03;
    options.ampy = 0.03;
    options.ampz = 0.03;

    poscar.displaceAtoms(static_cast<int>(cl_indices.size()), options);

    const std::set<size_t> cl_set(cl_indices.begin(), cl_indices.end());
    constexpr double tol = 1e-15;
    for (size_t i = 0; i < poscar.coordinates.size(); ++i) {
        const bool changed = std::abs(poscar.coordinates[i].x - original.coordinates[i].x) > tol ||
                             std::abs(poscar.coordinates[i].y - original.coordinates[i].y) > tol ||
                             std::abs(poscar.coordinates[i].z - original.coordinates[i].z) > tol;
        if (cl_set.count(i) > 0) {
            EXPECT_TRUE(changed) << "Expected Cl atom " << i << " to move";
        } else {
            EXPECT_FALSE(changed) << "Expected non-Cl atom " << i << " to stay fixed";
        }
    }
}

TEST_F(DisplacementTest, WrapFractionalCoordinates) {
    POSCAR p;
    p.comment = "wrap";
    p.scale = 1.0;
    p.lattice[0][0] = 1.0;
    p.lattice[0][1] = 0.0;
    p.lattice[0][2] = 0.0;
    p.lattice[1][0] = 0.0;
    p.lattice[1][1] = 1.0;
    p.lattice[1][2] = 0.0;
    p.lattice[2][0] = 0.0;
    p.lattice[2][1] = 0.0;
    p.lattice[2][2] = 1.0;
    p.elements = {"H"};
    p.num_atoms = {1};
    p.total_atoms = 1;
    p.is_direct = true;
    p.coordinates = {{1.2, -0.1, 2.0}};

    DisplacementOptions options;
    options.wrap_direct = true;
    options.ampx = 0.0;
    options.ampy = 0.0;
    options.ampz = 0.0;
    p.displaceAtoms(1, options);

    EXPECT_GE(p.coordinates[0].x, 0.0);
    EXPECT_LT(p.coordinates[0].x, 1.0);
    EXPECT_GE(p.coordinates[0].y, 0.0);
    EXPECT_LT(p.coordinates[0].y, 1.0);
    EXPECT_GE(p.coordinates[0].z, 0.0);
    EXPECT_LT(p.coordinates[0].z, 1.0);
    EXPECT_NEAR(p.coordinates[0].x, 0.2, 1e-12);
    EXPECT_NEAR(p.coordinates[0].y, 0.9, 1e-12);
    EXPECT_NEAR(p.coordinates[0].z, 0.0, 1e-12);
}

TEST_F(DisplacementTest, ZeroNetDisplacement) {
    seedRandom(1234);
    POSCAR p = poscar;
    p.toCartesian();
    const auto original = p.coordinates;

    DisplacementOptions options;
    options.zero_net = true;
    options.ampx = 0.05;
    options.ampy = 0.05;
    options.ampz = 0.05;
    p.displaceAtoms(p.total_atoms, options);

    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (size_t i = 0; i < p.coordinates.size(); ++i) {
        sx += p.coordinates[i].x - original[i].x;
        sy += p.coordinates[i].y - original[i].y;
        sz += p.coordinates[i].z - original[i].z;
    }

    EXPECT_NEAR(sx, 0.0, 1e-12);
    EXPECT_NEAR(sy, 0.0, 1e-12);
    EXPECT_NEAR(sz, 0.0, 1e-12);
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
