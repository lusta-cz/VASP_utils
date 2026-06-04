#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <string>

#include "eigenval2bands.h"

static const std::string kTmpDir = std::string(TEST_DATA_DIR);
constexpr double kTol = 1e-10;

// ─── Fixture File Helpers ────────────────────────────────────────────────────

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
    std::fputs(" 1 2 3\n", f);  // NIONS, NKPTS, NBANDS
    // kpt 1
    std::fputs("\n", f);
    std::fputs("  0.0 0.0 0.0  1.0\n", f);
    std::fputs(" 1  -10.1234\n", f);
    std::fputs(" 2   -5.5678\n", f);
    std::fputs(" 3    0.1234\n", f);
    // kpt 2
    std::fputs("\n", f);
    std::fputs("  0.5 0.0 0.0  1.0\n", f);
    std::fputs(" 1  -9.0000\n", f);
    std::fputs(" 2  -4.0000\n", f);
    std::fputs(" 3   1.5000\n", f);
    std::fclose(f);
}

// nspin=2, 1 k-point, 2 bands
static void writeEigenval2(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("header line 1\n", f);
    std::fputs("0.1 0.1 0.1 0.5 2.0\n", f);
    std::fputs("0.02\n", f);
    std::fputs("CAR\n", f);
    std::fputs("test system spin-polarized\n", f);
    std::fputs(" 1 1 2\n", f);  // NIONS, NKPTS, NBANDS
    // kpt 1
    std::fputs("\n", f);
    std::fputs("  0.0 0.0 0.0  1.0\n", f);
    // Band 1: Index, E_up, E_down, Occ_up, Occ_down
    std::fputs(" 1  -1.0000  -1.5000  1.0000  1.0000\n", f);
    // Band 2: Index, E_up, E_down, Occ_up, Occ_down
    std::fputs(" 2   2.0000   1.8000  0.0000  0.0000\n", f);
    std::fclose(f);
}

static void writeKpoints(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("High symmetry path test\n", f);
    std::fputs("10\n", f);  // 10 k-points per segment
    std::fputs("Line-mode\n", f);
    std::fputs("Reciprocal\n", f);
    std::fputs("  0.0 0.0 0.0 ! Gamma\n", f);
    std::fputs("  0.5 0.0 0.0 ! X\n", f);
    std::fclose(f);
}

static void writeDoscar(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs("l1\nl2\nl3\nl4\nl5\n", f);
    std::fputs("  10.0  0.0  1000  5.1234  1.000\n", f);  // Efermi is 4th value
    std::fclose(f);
}

static void writeOutcar(const char* path) {
    FILE* f = std::fopen(path, "w");
    std::fputs(" some lines before\n", f);
    std::fputs(" E-fermi :   4.9876     \n", f);
    std::fclose(f);
}

// ─── Test Fixture ────────────────────────────────────────────────────────────

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

// ─── parseEigenval ───────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseEigenvalSpin1Metadata) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    EXPECT_EQ(data.total_kpoints, 2);
    EXPECT_EQ(data.nbands, 3);
    EXPECT_EQ(data.nspin, 1);
}

TEST_F(EigenvalParseTest, ParseEigenvalSpin1Values) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));
    ASSERT_EQ(data.kpoints.size(), 2);

    EXPECT_NEAR(data.kpoints[0].x, 0.0, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[0], -10.1234, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[1], -5.5678, kTol);
}

TEST_F(EigenvalParseTest, ParseEigenvalSpin2Metadata) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    EXPECT_EQ(data.total_kpoints, 1);
    EXPECT_EQ(data.nbands, 2);
    EXPECT_EQ(data.nspin, 2);
}

TEST_F(EigenvalParseTest, ParseEigenvalSpin2Values) {
    EigenvalData data;
    ASSERT_TRUE(parseEigenval(eigenval2Path(), data));
    ASSERT_EQ(data.kpoints.size(), 1);

    EXPECT_NEAR(data.kpoints[0].energies_up[0], -1.0, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_dn[0], -1.5, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_up[1], 2.0, kTol);
    EXPECT_NEAR(data.kpoints[0].energies_dn[1], 1.8, kTol);
}

// ─── parseKpoints ─────────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, ParseKpointsReturnsTrue) {
    BZPath bz_path;
    EXPECT_TRUE(parseKpoints(kpointsPath(), bz_path));
}

TEST_F(EigenvalParseTest, ParseKpointsStructData) {
    BZPath bz_path;
    ASSERT_TRUE(parseKpoints(kpointsPath(), bz_path));
    EXPECT_EQ(bz_path.kpts_per_seg, 10);
    EXPECT_EQ(bz_path.num_segments, 1);
    ASSERT_EQ(bz_path.segments.size(), 1);
    EXPECT_EQ(bz_path.segments[0].start_label, "Gamma");
    EXPECT_EQ(bz_path.segments[0].end_label, "X");
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
    EXPECT_NEAR(e_fermi, 4.9876, kTol);
}

// ─── writeKpointsLog ──────────────────────────────────────────────────────────

TEST_F(EigenvalParseTest, WriteKpointsLogExecution) {
    BZPath bz_path;
    EigenvalData data;
    ASSERT_TRUE(parseKpoints(kpointsPath(), bz_path));
    ASSERT_TRUE(parseEigenval(eigenval1Path(), data));

    // Artificially expand data to fit the kpoints per segment requirements for the test
    data.total_kpoints = 10;
    data.kpoints.resize(10);
    for (int i = 0; i < 10; ++i) {
        data.kpoints[i].x = 0.05 * i;
        data.kpoints[i].y = 0.0;
        data.kpoints[i].z = 0.0;
    }

    double log_cumulative_dist = 0.0;
    std::string log_file = kTmpDir + "/test_kpoints.log";
    EXPECT_TRUE(writeKpointsLog(log_file, bz_path, data, log_cumulative_dist));
    EXPECT_GT(log_cumulative_dist, 0.0);
    std::remove(log_file.c_str());
}
