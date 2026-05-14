#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "elastic.h"
#include "poscar_file.h"

// ─── Helpers ────────────────────────────────────────────────────────────────

static POSCAR makeCubicPOSCAR(double a = 4.0) {
    POSCAR p;
    p.comment = "test";
    p.scale = 1.0;
    p.lattice[0][0] = a;
    p.lattice[0][1] = 0;
    p.lattice[0][2] = 0;
    p.lattice[1][0] = 0;
    p.lattice[1][1] = a;
    p.lattice[1][2] = 0;
    p.lattice[2][0] = 0;
    p.lattice[2][1] = 0;
    p.lattice[2][2] = a;
    p.elements = {"Cu"};
    p.num_atoms = {1};
    p.total_atoms = 1;
    p.is_direct = true;
    Atom at;
    at.x = 0.0;
    at.y = 0.0;
    at.z = 0.0;
    p.coordinates.push_back(at);
    return p;
}

// ─── crystalSystemFromSpaceGroup ────────────────────────────────────────────

TEST(CrystalSystemFromSpaceGroup, TriclinicBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(1), CrystalSystem::Triclinic);
    EXPECT_EQ(crystalSystemFromSpaceGroup(2), CrystalSystem::Triclinic);
}

TEST(CrystalSystemFromSpaceGroup, MonoclinicBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(3), CrystalSystem::Monoclinic);
    EXPECT_EQ(crystalSystemFromSpaceGroup(15), CrystalSystem::Monoclinic);
}

TEST(CrystalSystemFromSpaceGroup, OrthorhombicBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(16), CrystalSystem::Orthorhombic);
    EXPECT_EQ(crystalSystemFromSpaceGroup(74), CrystalSystem::Orthorhombic);
}

TEST(CrystalSystemFromSpaceGroup, TetragonalBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(75), CrystalSystem::TetragonalI);
    EXPECT_EQ(crystalSystemFromSpaceGroup(88), CrystalSystem::TetragonalI);
    EXPECT_EQ(crystalSystemFromSpaceGroup(89), CrystalSystem::TetragonalII);
    EXPECT_EQ(crystalSystemFromSpaceGroup(142), CrystalSystem::TetragonalII);
}

TEST(CrystalSystemFromSpaceGroup, TrigonalBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(143), CrystalSystem::TrigonalI);
    EXPECT_EQ(crystalSystemFromSpaceGroup(148), CrystalSystem::TrigonalI);
    EXPECT_EQ(crystalSystemFromSpaceGroup(149), CrystalSystem::TrigonalII);
    EXPECT_EQ(crystalSystemFromSpaceGroup(167), CrystalSystem::TrigonalII);
}

TEST(CrystalSystemFromSpaceGroup, HexagonalBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(168), CrystalSystem::Hexagonal);
    EXPECT_EQ(crystalSystemFromSpaceGroup(194), CrystalSystem::Hexagonal);
}

TEST(CrystalSystemFromSpaceGroup, CubicBoundaries) {
    EXPECT_EQ(crystalSystemFromSpaceGroup(195), CrystalSystem::Cubic);
    EXPECT_EQ(crystalSystemFromSpaceGroup(225), CrystalSystem::Cubic);  // NaCl Fm-3m
    EXPECT_EQ(crystalSystemFromSpaceGroup(230), CrystalSystem::Cubic);
}

TEST(CrystalSystemFromSpaceGroup, OutOfRangeThrows) {
    EXPECT_THROW(crystalSystemFromSpaceGroup(0), std::out_of_range);
    EXPECT_THROW(crystalSystemFromSpaceGroup(231), std::out_of_range);
    EXPECT_THROW(crystalSystemFromSpaceGroup(-1), std::out_of_range);
}

// ─── crystalSystemName ──────────────────────────────────────────────────────

TEST(CrystalSystemName, AllSystemsReturnNonEmpty) {
    const CrystalSystem systems[] = {
        CrystalSystem::Triclinic,   CrystalSystem::Monoclinic,   CrystalSystem::Orthorhombic,
        CrystalSystem::TetragonalI, CrystalSystem::TetragonalII, CrystalSystem::TrigonalI,
        CrystalSystem::TrigonalII,  CrystalSystem::Hexagonal,    CrystalSystem::Cubic,
    };
    for (auto cs : systems)
        EXPECT_FALSE(crystalSystemName(cs).empty());
}

TEST(CrystalSystemName, KnownValues) {
    EXPECT_EQ(crystalSystemName(CrystalSystem::Cubic), "cubic");
    EXPECT_EQ(crystalSystemName(CrystalSystem::Hexagonal), "hexagonal");
    EXPECT_EQ(crystalSystemName(CrystalSystem::Triclinic), "triclinic");
    EXPECT_EQ(crystalSystemName(CrystalSystem::Monoclinic), "monoclinic");
    EXPECT_EQ(crystalSystemName(CrystalSystem::Orthorhombic), "orthorhombic");
}

// ─── nIndependentConstants ──────────────────────────────────────────────────

TEST(NIndependentConstants, CorrectCounts) {
    EXPECT_EQ(nIndependentConstants(CrystalSystem::Triclinic), 21);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::Monoclinic), 13);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::Orthorhombic), 9);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::TetragonalI), 7);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::TetragonalII), 6);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::TrigonalI), 7);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::TrigonalII), 6);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::Hexagonal), 5);
    EXPECT_EQ(nIndependentConstants(CrystalSystem::Cubic), 3);
}

// ─── energyStrainModes ──────────────────────────────────────────────────────

TEST(EnergyStrainModes, ModeCounts) {
    EXPECT_EQ(energyStrainModes(CrystalSystem::Cubic).size(), 3u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::Hexagonal).size(), 5u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::TrigonalII).size(), 6u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::TrigonalI).size(), 7u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::TetragonalII).size(), 6u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::TetragonalI).size(), 7u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::Orthorhombic).size(), 9u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::Monoclinic).size(), 13u);
    EXPECT_EQ(energyStrainModes(CrystalSystem::Triclinic).size(), 21u);
}

TEST(EnergyStrainModes, CubicVoigtVectors) {
    auto modes = energyStrainModes(CrystalSystem::Cubic);

    // "vol" mode: e1=e2=e3=1, shear=0
    EXPECT_EQ(modes[0].label, "vol");
    EXPECT_DOUBLE_EQ(modes[0].voigt[0], 1.0);
    EXPECT_DOUBLE_EQ(modes[0].voigt[1], 1.0);
    EXPECT_DOUBLE_EQ(modes[0].voigt[2], 1.0);
    EXPECT_DOUBLE_EQ(modes[0].voigt[3], 0.0);
    EXPECT_DOUBLE_EQ(modes[0].voigt[4], 0.0);
    EXPECT_DOUBLE_EQ(modes[0].voigt[5], 0.0);
    EXPECT_FALSE(modes[0].symmetric);

    // "e1-e2" mode: symmetric = true
    EXPECT_EQ(modes[1].label, "e1-e2");
    EXPECT_DOUBLE_EQ(modes[1].voigt[0], 1.0);
    EXPECT_DOUBLE_EQ(modes[1].voigt[1], -1.0);
    EXPECT_TRUE(modes[1].symmetric);

    // "e4" shear mode: symmetric = true
    EXPECT_EQ(modes[2].label, "e4");
    EXPECT_DOUBLE_EQ(modes[2].voigt[3], 1.0);
    EXPECT_TRUE(modes[2].symmetric);
}

TEST(EnergyStrainModes, TriclinicHasSixPureModes) {
    auto modes = energyStrainModes(CrystalSystem::Triclinic);
    // First 6 are pure (single-component) modes
    for (int i = 0; i < 6; ++i) {
        int nonzero = 0;
        for (int k = 0; k < 6; ++k)
            if (modes[i].voigt[k] != 0.0)
                nonzero++;
        EXPECT_EQ(nonzero, 1) << "Pure mode " << i << " should have exactly one nonzero Voigt component";
    }
}

TEST(EnergyStrainModes, TriclinicHas15MixedModes) {
    auto modes = energyStrainModes(CrystalSystem::Triclinic);
    // Modes 6–20 are mixed (two-component) modes
    for (size_t i = 6; i < modes.size(); ++i) {
        int nonzero = 0;
        for (int k = 0; k < 6; ++k)
            if (modes[i].voigt[k] != 0.0)
                nonzero++;
        EXPECT_EQ(nonzero, 2) << "Mixed mode " << i << " should have exactly two nonzero Voigt components";
    }
}

TEST(EnergyStrainModes, AllLabelsNonEmpty) {
    for (auto cs :
         {CrystalSystem::Cubic, CrystalSystem::Hexagonal, CrystalSystem::Orthorhombic, CrystalSystem::Triclinic}) {
        for (const auto& m : energyStrainModes(cs))
            EXPECT_FALSE(m.label.empty());
    }
}

// ─── stressStrainModes ──────────────────────────────────────────────────────

TEST(StressStrainModes, AlwaysSixModes) {
    EXPECT_EQ(stressStrainModes().size(), 6u);
}

TEST(StressStrainModes, VoigtBasisVectors) {
    auto modes = stressStrainModes();
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 6; ++k)
            EXPECT_DOUBLE_EQ(modes[i].voigt[k], (k == i) ? 1.0 : 0.0) << "mode " << i << " component " << k;
    }
}

TEST(StressStrainModes, Labels) {
    auto modes = stressStrainModes();
    const std::string expected[] = {"e1", "e2", "e3", "e4", "e5", "e6"};
    for (int i = 0; i < 6; ++i)
        EXPECT_EQ(modes[i].label, expected[i]);
}

// ─── applyStrain ────────────────────────────────────────────────────────────

TEST(ApplyStrain, ZeroAmplitudeLeavesLatticeUnchanged) {
    POSCAR p = makeCubicPOSCAR(4.0);
    std::array<double, 6> voigt = {1, 0, 0, 0, 0, 0};
    POSCAR out = applyStrain(p, voigt, 0.0);

    constexpr double tol = 1e-12;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(out.lattice[i][j], p.lattice[i][j], tol);
}

TEST(ApplyStrain, PureLongitudinalStrainE1) {
    POSCAR p = makeCubicPOSCAR(4.0);
    std::array<double, 6> voigt = {1, 0, 0, 0, 0, 0};  // e1 only
    double amp = 0.01;
    POSCAR out = applyStrain(p, voigt, amp);

    constexpr double tol = 1e-10;
    // a-axis stretched by (1 + amp)
    EXPECT_NEAR(out.lattice[0][0], 4.0 * (1.0 + amp), tol);
    // b and c axes unchanged
    EXPECT_NEAR(out.lattice[1][1], 4.0, tol);
    EXPECT_NEAR(out.lattice[2][2], 4.0, tol);
    // Off-diagonals remain zero
    EXPECT_NEAR(out.lattice[0][1], 0.0, tol);
    EXPECT_NEAR(out.lattice[0][2], 0.0, tol);
}

TEST(ApplyStrain, VolumetricStrainScalesAllAxes) {
    POSCAR p = makeCubicPOSCAR(4.0);
    std::array<double, 6> voigt = {1, 1, 1, 0, 0, 0};  // vol mode
    double amp = 0.01;
    POSCAR out = applyStrain(p, voigt, amp);

    constexpr double tol = 1e-10;
    EXPECT_NEAR(out.lattice[0][0], 4.0 * (1.0 + amp), tol);
    EXPECT_NEAR(out.lattice[1][1], 4.0 * (1.0 + amp), tol);
    EXPECT_NEAR(out.lattice[2][2], 4.0 * (1.0 + amp), tol);
}

TEST(ApplyStrain, ShearStrainE4IntroducesOffDiagonal) {
    // e4 = 2*ε_yz, so off-diagonal lattice element = amp/2 * a
    POSCAR p = makeCubicPOSCAR(4.0);
    std::array<double, 6> voigt = {0, 0, 0, 1, 0, 0};  // e4 (yz shear)
    double amp = 0.02;
    POSCAR out = applyStrain(p, voigt, amp);

    constexpr double tol = 1e-10;
    // Diagonal elements unchanged
    EXPECT_NEAR(out.lattice[0][0], 4.0, tol);
    EXPECT_NEAR(out.lattice[1][1], 4.0, tol);
    EXPECT_NEAR(out.lattice[2][2], 4.0, tol);
    // e4 = 2*ε_yz → ε_yz = amp/2; lattice[1] gets component along z, lattice[2] gets component along y
    // new_lattice[1][2] = old_lattice[1][1] * ε_yz = 4 * (amp/2)
    EXPECT_NEAR(out.lattice[1][2], 4.0 * (amp / 2.0), tol);
    EXPECT_NEAR(out.lattice[2][1], 4.0 * (amp / 2.0), tol);
}

TEST(ApplyStrain, AtomicCoordinatesPreserved) {
    POSCAR p = makeCubicPOSCAR(4.0);
    Atom extra;
    extra.x = 0.5;
    extra.y = 0.25;
    extra.z = 0.75;
    p.coordinates.push_back(extra);
    p.total_atoms = 2;

    std::array<double, 6> voigt = {1, 1, 1, 0, 0, 0};
    POSCAR out = applyStrain(p, voigt, 0.01);

    ASSERT_EQ(out.coordinates.size(), 2u);
    EXPECT_DOUBLE_EQ(out.coordinates[1].x, 0.5);
    EXPECT_DOUBLE_EQ(out.coordinates[1].y, 0.25);
    EXPECT_DOUBLE_EQ(out.coordinates[1].z, 0.75);
}

// ─── writeElasticLog / readElasticLog ───────────────────────────────────────

static ElasticDeformLog makeSampleLog() {
    ElasticDeformLog log;
    log.method = "energy";
    log.amplitude = 0.01;
    log.npoints = 5;
    log.symmetry_mode = "auto";
    log.crystal_system = CrystalSystem::Cubic;
    log.space_group_number = 225;
    log.space_group_symbol = "Fm-3m";
    log.point_group = "m-3m";
    log.n_independent = 3;
    log.volume = 64.0;
    log.reference_dir = "ref";

    ElasticStrainMode m;
    m.label = "vol";
    m.voigt = {1, 1, 1, 0, 0, 0};
    m.symmetric = false;
    log.modes.push_back(m);
    log.amplitudes.push_back({-0.01, 0.0, 0.01});
    log.dirs.push_back({"vol_m", "vol_0", "vol_p"});
    return log;
}

TEST(ElasticLog, RoundTrip) {
    const std::string tmp = std::string(std::tmpnam(nullptr));
    ElasticDeformLog orig = makeSampleLog();
    ASSERT_TRUE(writeElasticLog(tmp, orig));

    auto result = readElasticLog(tmp);
    ASSERT_TRUE(result.has_value());

    EXPECT_EQ(result->method, orig.method);
    EXPECT_DOUBLE_EQ(result->amplitude, orig.amplitude);
    EXPECT_EQ(result->npoints, orig.npoints);
    EXPECT_EQ(result->symmetry_mode, orig.symmetry_mode);
    EXPECT_EQ(result->crystal_system, orig.crystal_system);
    EXPECT_EQ(result->space_group_number, orig.space_group_number);
    EXPECT_EQ(result->space_group_symbol, orig.space_group_symbol);
    EXPECT_EQ(result->point_group, orig.point_group);
    EXPECT_EQ(result->n_independent, orig.n_independent);
    EXPECT_DOUBLE_EQ(result->volume, orig.volume);
    EXPECT_EQ(result->reference_dir, orig.reference_dir);

    ASSERT_EQ(result->modes.size(), 1u);
    EXPECT_EQ(result->modes[0].label, "vol");
    for (int k = 0; k < 6; ++k)
        EXPECT_DOUBLE_EQ(result->modes[0].voigt[k], orig.modes[0].voigt[k]);

    ASSERT_EQ(result->amplitudes[0].size(), 3u);
    EXPECT_NEAR(result->amplitudes[0][0], -0.01, 1e-9);
    EXPECT_NEAR(result->amplitudes[0][2], 0.01, 1e-9);

    ASSERT_EQ(result->dirs[0].size(), 3u);
    EXPECT_EQ(result->dirs[0][0], "vol_m");
    EXPECT_EQ(result->dirs[0][2], "vol_p");

    std::remove(tmp.c_str());
}

TEST(ElasticLog, NonExistentFileReturnsNullopt) {
    auto result = readElasticLog("/nonexistent/path/elastic.log");
    EXPECT_FALSE(result.has_value());
}

TEST(ElasticLog, EmptyMethodReturnsNullopt) {
    const std::string tmp = std::string(std::tmpnam(nullptr));
    // Write a file with N_MODES but no METHOD
    {
        std::ofstream f(tmp);
        f << "N_MODES=1\n";
        f << "MODE_1_LABEL=vol\n";
        f << "MODE_1_VOIGT=1 1 1 0 0 0\n";
        f << "MODE_1_AMPLITUDES=+0.010000\n";
        f << "MODE_1_DIRS=vol_p\n";
    }
    auto result = readElasticLog(tmp);
    EXPECT_FALSE(result.has_value());
    std::remove(tmp.c_str());
}

TEST(ElasticLog, ZeroModesReturnsNullopt) {
    const std::string tmp = std::string(std::tmpnam(nullptr));
    {
        std::ofstream f(tmp);
        f << "METHOD=energy\n";
        // N_MODES intentionally omitted (defaults to 0)
    }
    auto result = readElasticLog(tmp);
    EXPECT_FALSE(result.has_value());
    std::remove(tmp.c_str());
}

TEST(ElasticLog, AllCrystalSystemsRoundTrip) {
    const CrystalSystem systems[] = {
        CrystalSystem::Triclinic,   CrystalSystem::Monoclinic,   CrystalSystem::Orthorhombic,
        CrystalSystem::TetragonalI, CrystalSystem::TetragonalII, CrystalSystem::TrigonalI,
        CrystalSystem::TrigonalII,  CrystalSystem::Hexagonal,    CrystalSystem::Cubic,
    };
    for (auto cs : systems) {
        ElasticDeformLog log = makeSampleLog();
        log.crystal_system = cs;

        const std::string tmp = std::string(std::tmpnam(nullptr));
        ASSERT_TRUE(writeElasticLog(tmp, log));
        auto result = readElasticLog(tmp);
        ASSERT_TRUE(result.has_value());
        EXPECT_EQ(result->crystal_system, cs);
        std::remove(tmp.c_str());
    }
}
