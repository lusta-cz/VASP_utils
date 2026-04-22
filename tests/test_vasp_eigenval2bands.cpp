#include <gtest/gtest.h>

#include <cstdio>
#include <string>

#include "eigenval2bands.h"

static const std::string kTmpDir = std::string(TEST_DATA_DIR);
constexpr double kTol = 1e-10;

// ─── fixture file helpers ────────────────────────────────────────────────────

static std::string eigenval1Path() {
    return kTmpDir + "/test_ev2b_spin1_tmp.eigenval";
}
static std::string eigenval2Path() {
    return kTmpDir + "/test_ev2b_spin2_tmp.eigenval";
}
static std::string kpointsPath() {
    return kTmpDir + "/test_ev2b_kpoints_tmp.kpoints";
}
static std::string doscarPath() {
    return kTmpDir + "/test_ev2b_doscar_tmp.doscar";
}
static std::string outcarPath() {
    return kTmpDir + "/test_ev2b_outcar_tmp.outcar";
}

// nspin=1, 2 k-points, 3 bands
static void writeEigenval1(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("header line 1\n", f);
    std::fputs("0.1 0.1 0.1 0.5 2.0\n", f);
    std::fputs("0.02\n", f);
    std::fputs("CAR\n", f);
    std::fputs("test system\n", f);
    std::fputs("1 2 3\n", f);  // nions=1, nkpts=2, nbands=3
    std::fputs("\n", f);
    std::fputs("0.000000 0.000000 0.000000 0.250000\n", f);
    std::fputs("1 -13.5821\n", f);
    std::fputs("2  -4.7293\n", f);
    std::fputs("3   3.0000\n", f);
    std::fputs("\n", f);
    std::fputs("0.500000 0.000000 0.000000 0.250000\n", f);
    std::fputs("1 -13.1234\n", f);
    std::fputs("2  -4.5000\n", f);
    std::fputs("3   3.5000\n", f);
    std::fclose(f);
}

// nspin=2, 1 k-point, 2 bands — band lines carry: idx e_up e_dn occ_up occ_dn
static void writeEigenval2(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("header line 1\n", f);
    std::fputs("0.1 0.1 0.1 0.5 2.0\n", f);
    std::fputs("0.02\n", f);
    std::fputs("CAR\n", f);
    std::fputs("test system\n", f);
    std::fputs("1 1 2\n", f);  // nions=1, nkpts=1, nbands=2
    std::fputs("\n", f);
    std::fputs("0.000000 0.000000 0.000000 0.500000\n", f);
    std::fputs("1 -10.2345 -10.1234 1.00000 1.00000\n", f);
    std::fputs("2  -5.6789  -5.5678 0.50000 0.48000\n", f);
    std::fclose(f);
}

static void writeKpoints(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("K-path for band structure\n", f);
    std::fputs("10\n", f);         // 10 k-points per segment
    std::fputs("Line-mode\n", f);  // starts with 'L'
    std::fputs("reciprocal\n", f);
    std::fputs("0.000 0.000 0.000 ! Gamma\n", f);
    std::fputs("0.500 0.000 0.000 ! X\n", f);
    std::fclose(f);
}

static void writeDoscar(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("1 1 1 1 1\n", f);
    std::fputs("CAR\n", f);
    std::fputs("0.02\n", f);
    std::fputs("NELEM\n", f);
    std::fputs("comment\n", f);
    std::fputs("20.0000 -20.0000 301 5.1234 1.0000\n", f);  // Efermi = 5.1234
    std::fclose(f);
}

// Two E-fermi entries; parser must return the last one (7.8901)
static void writeOutcar(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs(" VASP output\n", f);
    std::fputs("   E-fermi :   3.1415   XC(G=0): -0.1234   alpha+bet: -0.01\n", f);
    std::fputs(" more SCF output\n", f);
    std::fputs("   E-fermi :   7.8901   XC(G=0): -0.1234   alpha+bet: -0.01\n", f);
    std::fputs(" final output\n", f);
    std::fclose(f);
}

// ─── Fixture ─────────────────────────────────────────────────────────────────

class EigenvalParseTest : public ::testing::Test {
protected:
    void SetUp() override {
        writeEigenval1(eigenval1Path().c_str());
        writeEigenval2(eigenval2Path().c_str());
        writeKpoints(kpointsPath().c_str());
        writeDoscar(doscarPath().c_str());
        writeOutcar(outcarPath().c_str());
    }
    void TearDown() override {
        std::remove(eigenval1Path().c_str());
        std::remove(eigenval2Path().c_str());
        std::remove(kpointsPath().c_str());
        std::remove(doscarPath().c_str());
        std::remove(outcarPath().c_str());
    }
};

// ─── parseEigenval: nspin=1 ───────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseEigenvalReturnsTrue) {
    EigenvalData data;
    EXPECT_TRUE(parseEigenval(eigenval1Path(), data));
}

TEST_F(EigenvalParseTest, ParseEigenvalTotalKpoints) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    EXPECT_EQ(data.total_kpoints, 2);
}

TEST_F(EigenvalParseTest, ParseEigenvalNbands) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    EXPECT_EQ(data.nbands, 3);
}

TEST_F(EigenvalParseTest, ParseEigenvalNspin1Detected) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    EXPECT_EQ(data.nspin, 1);
}

TEST_F(EigenvalParseTest, ParseEigenvalKpointVectorSize) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    EXPECT_EQ(data.kpoints.size(), 2u);
}

TEST_F(EigenvalParseTest, ParseEigenvalFirstKpointCoords) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    ASSERT_GE(data.kpoints.size(), 1u);
    EXPECT_NEAR(data.kpoints[0].x, 0.0, kTol);
    EXPECT_NEAR(data.kpoints[0].y, 0.0, kTol);
    EXPECT_NEAR(data.kpoints[0].z, 0.0, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalSecondKpointCoords) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    ASSERT_GE(data.kpoints.size(), 2u);
    EXPECT_NEAR(data.kpoints[1].x, 0.5, kTol);
    EXPECT_NEAR(data.kpoints[1].y, 0.0, kTol);
    EXPECT_NEAR(data.kpoints[1].z, 0.0, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalFirstKpointBandEnergies) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    ASSERT_GE(data.kpoints.size(), 1u);
    ASSERT_EQ(data.kpoints[0].energies_up.size(), 3u);
    EXPECT_NEAR(data.kpoints[0].energies_up[0], -13.5821, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[1], -4.7293, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[2], 3.0000, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalSecondKpointBandEnergies) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    ASSERT_GE(data.kpoints.size(), 2u);
    ASSERT_EQ(data.kpoints[1].energies_up.size(), 3u);
    EXPECT_NEAR(data.kpoints[1].energies_up[0], -13.1234, kTol);
    EXPECT_NEAR(data.kpoints[1].energies_up[1], -4.5000, kTol);
    EXPECT_NEAR(data.kpoints[1].energies_up[2], 3.5000, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalNspin1SpinDownEmpty) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    for (const auto& kp : data.kpoints) {
        EXPECT_TRUE(kp.energies_dn.empty());
    }
}

// ─── parseEigenval: nspin=2 ───────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseEigenvalNspin2Detected) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    EXPECT_EQ(data.nspin, 2);
}

TEST_F(EigenvalParseTest, ParseEigenvalNspin2Metadata) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    EXPECT_EQ(data.total_kpoints, 1);
    EXPECT_EQ(data.nbands, 2);
}

TEST_F(EigenvalParseTest, ParseEigenvalNspin2SpinUpEnergies) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    ASSERT_EQ(data.kpoints.size(), 1u);
    ASSERT_EQ(data.kpoints[0].energies_up.size(), 2u);
    EXPECT_NEAR(data.kpoints[0].energies_up[0], -10.2345, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[1], -5.6789, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalNspin2SpinDownEnergies) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    ASSERT_EQ(data.kpoints.size(), 1u);
    ASSERT_EQ(data.kpoints[0].energies_dn.size(), 2u);
    EXPECT_NEAR(data.kpoints[0].energies_dn[0], -10.1234, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_dn[1], -5.5678, kTol);
}

// ─── parseKpoints ─────────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseKpointsReturnsTrue) {
    int kpts = 0;
    EXPECT_TRUE(parseKpoints(kpointsPath(), kpts));
}

TEST_F(EigenvalParseTest, ParseKpointsPerSegmentValue) {
    int kpts = 0;
    ASSERT_TRUE(parseKpoints(kpointsPath(), kpts));
    EXPECT_EQ(kpts, 10);
}

// ─── parseFromDoscar ──────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseFromDoscarReturnsTrue) {
    double e_fermi = 0.0;
    EXPECT_TRUE(parseFromDoscar(doscarPath(), e_fermi));
}

TEST_F(EigenvalParseTest, ParseFromDoscarFermiLevel) {
    double e_fermi = 0.0;
    ASSERT_TRUE(parseFromDoscar(doscarPath(), e_fermi));
    EXPECT_NEAR(e_fermi, 5.1234, kTol);
}

// ─── parseFromOutcar ──────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseFromOutcarReturnsTrue) {
    double e_fermi = 0.0;
    EXPECT_TRUE(parseFromOutcar(outcarPath(), e_fermi));
}

TEST_F(EigenvalParseTest, ParseFromOutcarFermiLevel) {
    double e_fermi = 0.0;
    ASSERT_TRUE(parseFromOutcar(outcarPath(), e_fermi));
    EXPECT_NEAR(e_fermi, 7.8901, kTol);
}

TEST_F(EigenvalParseTest, ParseFromOutcarUsesLastOccurrence) {
    // OUTCAR has E-fermi = 3.1415 first, then 7.8901 last — only the last counts
    double e_fermi = 0.0;
    ASSERT_TRUE(parseFromOutcar(outcarPath(), e_fermi));
    EXPECT_NEAR(e_fermi, 7.8901, kTol);
    EXPECT_GT(e_fermi, 5.0);
}
