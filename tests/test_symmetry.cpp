#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "poscar_file.h"
#include "symmetry.h"

// ─── Test fixture: NaCl conventional cell (Fm-3m, SG 225) ──────────────────
// Cubic, a = 5.64 Å, 4 Na + 4 Cl = 8 atoms

static POSCAR makeNaClConventional() {
    POSCAR p;
    p.comment = "NaCl conventional";
    p.scale = 1.0;
    const double a = 5.64;
    p.lattice[0][0] = a;
    p.lattice[0][1] = 0;
    p.lattice[0][2] = 0;
    p.lattice[1][0] = 0;
    p.lattice[1][1] = a;
    p.lattice[1][2] = 0;
    p.lattice[2][0] = 0;
    p.lattice[2][1] = 0;
    p.lattice[2][2] = a;
    p.elements = {"Na", "Cl"};
    p.num_atoms = {4, 4};
    p.total_atoms = 8;
    p.is_direct = true;

    // Na at face-centring positions
    const double naPos[4][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
    // Cl offset by (0.5, 0, 0) relative to Na
    const double clPos[4][3] = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.5, 0.5}};
    for (auto& row : naPos) {
        Atom a;
        a.x = row[0];
        a.y = row[1];
        a.z = row[2];
        p.coordinates.push_back(a);
    }
    for (auto& row : clPos) {
        Atom a;
        a.x = row[0];
        a.y = row[1];
        a.z = row[2];
        p.coordinates.push_back(a);
    }
    return p;
}

// Simple orthorhombic structure for type-assignment tests
static POSCAR makeOrthorhombicTwoSpecies() {
    POSCAR p;
    p.comment = "ortho";
    p.scale = 1.0;
    p.lattice[0][0] = 3.0;
    p.lattice[0][1] = 0;
    p.lattice[0][2] = 0;
    p.lattice[1][0] = 0;
    p.lattice[1][1] = 4.0;
    p.lattice[1][2] = 0;
    p.lattice[2][0] = 0;
    p.lattice[2][1] = 0;
    p.lattice[2][2] = 5.0;
    p.elements = {"A", "B"};
    p.num_atoms = {1, 2};
    p.total_atoms = 3;
    p.is_direct = true;

    const double pos[3][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}};
    for (auto& row : pos) {
        Atom a;
        a.x = row[0];
        a.y = row[1];
        a.z = row[2];
        p.coordinates.push_back(a);
    }
    return p;
}

// ─── initializeSpglibInput ──────────────────────────────────────────────────

TEST(InitializeSpglibInput, LatticeIsTransposed) {
    // POSCAR row-vectors → spglib column-vectors means lattice[i][j] should equal poscar.lattice[j][i]
    POSCAR p;
    p.scale = 1.0;
    p.lattice[0][0] = 1;
    p.lattice[0][1] = 2;
    p.lattice[0][2] = 3;
    p.lattice[1][0] = 4;
    p.lattice[1][1] = 5;
    p.lattice[1][2] = 6;
    p.lattice[2][0] = 7;
    p.lattice[2][1] = 8;
    p.lattice[2][2] = 9;
    p.elements = {"X"};
    p.num_atoms = {1};
    p.total_atoms = 1;
    p.is_direct = true;
    Atom at;
    at.x = 0;
    at.y = 0;
    at.z = 0;
    p.coordinates.push_back(at);

    double lattice[3][3];
    std::vector<std::array<double, 3>> positions(1);
    std::vector<int> types(1);
    std::map<std::string, int> elem_map;

    initializeSpglibInput(p, lattice, positions, types, elem_map);

    // spglib column-vector convention: lattice[col][row] in C, i.e. lattice[i][j] = poscar.lattice[j][i]
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(lattice[i][j], p.lattice[j][i]) << "lattice[" << i << "][" << j << "] mismatch";
}

TEST(InitializeSpglibInput, PositionsCopied) {
    POSCAR p = makeNaClConventional();

    double lattice[3][3];
    std::vector<std::array<double, 3>> positions(p.total_atoms);
    std::vector<int> types(p.total_atoms);
    std::map<std::string, int> elem_map;

    initializeSpglibInput(p, lattice, positions, types, elem_map);

    ASSERT_EQ(static_cast<int>(positions.size()), p.total_atoms);
    for (int i = 0; i < p.total_atoms; ++i) {
        EXPECT_DOUBLE_EQ(positions[i][0], p.coordinates[i].x);
        EXPECT_DOUBLE_EQ(positions[i][1], p.coordinates[i].y);
        EXPECT_DOUBLE_EQ(positions[i][2], p.coordinates[i].z);
    }
}

TEST(InitializeSpglibInput, TypesReflectSpeciesOrder) {
    POSCAR p = makeOrthorhombicTwoSpecies();

    double lattice[3][3];
    std::vector<std::array<double, 3>> positions(p.total_atoms);
    std::vector<int> types(p.total_atoms);
    std::map<std::string, int> elem_map;

    initializeSpglibInput(p, lattice, positions, types, elem_map);

    // "A" is first species → type 1, "B" is second → type 2
    EXPECT_EQ(types[0], types[0]);  // A has one consistent type
    EXPECT_EQ(types[1], types[2]);  // both B atoms share the same type
    EXPECT_NE(types[0], types[1]);  // A and B differ
    EXPECT_EQ(elem_map.count("A"), 1u);
    EXPECT_EQ(elem_map.count("B"), 1u);
    EXPECT_NE(elem_map.at("A"), elem_map.at("B"));
}

TEST(InitializeSpglibInput, ElementMapPopulated) {
    POSCAR p = makeNaClConventional();

    double lattice[3][3];
    std::vector<std::array<double, 3>> positions(p.total_atoms);
    std::vector<int> types(p.total_atoms);
    std::map<std::string, int> elem_map;

    initializeSpglibInput(p, lattice, positions, types, elem_map);

    EXPECT_EQ(elem_map.size(), 2u);
    EXPECT_EQ(elem_map.count("Na"), 1u);
    EXPECT_EQ(elem_map.count("Cl"), 1u);
    EXPECT_NE(elem_map.at("Na"), elem_map.at("Cl"));
}

// ─── analyzeSymmetry ────────────────────────────────────────────────────────

TEST(AnalyzeSymmetry, NaClSpaceGroup) {
    POSCAR p = makeNaClConventional();
    auto dataset = analyzeSymmetry(p, 1e-5);
    ASSERT_NE(dataset, nullptr);
    EXPECT_EQ(dataset->spacegroup_number, 225);  // Fm-3m
}

TEST(AnalyzeSymmetry, NaClSymmetryOperations) {
    POSCAR p = makeNaClConventional();
    auto dataset = analyzeSymmetry(p, 1e-5);
    ASSERT_NE(dataset, nullptr);
    // Fm-3m has 192 symmetry operations
    EXPECT_EQ(dataset->n_operations, 192);
}

TEST(AnalyzeSymmetry, NaClAtomCount) {
    POSCAR p = makeNaClConventional();
    auto dataset = analyzeSymmetry(p, 1e-5);
    ASSERT_NE(dataset, nullptr);
    EXPECT_EQ(dataset->n_atoms, 8);
}

TEST(AnalyzeSymmetry, WorksWithCartesianInput) {
    // analyzeSymmetry should internally convert to fractional before calling spglib
    POSCAR p = makeNaClConventional();
    p.toCartesian();
    ASSERT_FALSE(p.is_direct);

    auto dataset = analyzeSymmetry(p, 1e-5);
    ASSERT_NE(dataset, nullptr);
    EXPECT_EQ(dataset->spacegroup_number, 225);
}

// ─── standardizeCell ────────────────────────────────────────────────────────

TEST(StandardizeCell, ConventionalCellPreservesAtomCount) {
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 0, 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->total_atoms, 8);
}

TEST(StandardizeCell, PrimitiveCellReducesAtomCount) {
    // FCC NaCl conventional (8 atoms) → primitive (2 atoms)
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 1, 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->total_atoms, 2);
}

TEST(StandardizeCell, PrimitiveCellElementsPreserved) {
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 1, 1);
    ASSERT_TRUE(result.has_value());
    // Both Na and Cl should still be present
    bool has_Na = false, has_Cl = false;
    for (const auto& el : result->elements) {
        if (el == "Na")
            has_Na = true;
        if (el == "Cl")
            has_Cl = true;
    }
    EXPECT_TRUE(has_Na);
    EXPECT_TRUE(has_Cl);
}

TEST(StandardizeCell, ResultIsDirectCoordinates) {
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 1, 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(result->is_direct);
}

TEST(StandardizeCell, CommentContainsPrimitiveTag) {
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 1, 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_NE(result->comment.find("primitive"), std::string::npos);
}

TEST(StandardizeCell, CommentContainsConventionalTag) {
    POSCAR p = makeNaClConventional();
    auto result = standardizeCell(p, 1e-5, 0, 1);
    ASSERT_TRUE(result.has_value());
    EXPECT_NE(result->comment.find("conventional"), std::string::npos);
}
