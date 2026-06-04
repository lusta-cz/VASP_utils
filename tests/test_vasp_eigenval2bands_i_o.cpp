#include <gtest/gtest.h>

#include <cstdio>
#include <string>

#include "eigenval2bands.h"

static const std::string kTmpDir = std::string(TEST_DATA_DIR);

// ─── Non-existent file tests ──────────────────────────────────────────────────

TEST(EigenvalIO, ParseEigenvalNonExistentFile) {
    EigenvalData data;
    EXPECT_FALSE(parseEigenval("nonexistent_eigenval.dat", data));
}

TEST(EigenvalIO, ParseKpointsNonExistentFile) {
    BZPath bz_path;
    EXPECT_FALSE(parseKpoints("nonexistent_kpoints.dat", bz_path));
}

TEST(EigenvalIO, ParseFromDoscarNonExistentFile) {
    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromDoscar("nonexistent_doscar.dat", e_fermi));
}

TEST(EigenvalIO, ParseFromOutcarNonExistentFile) {
    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromOutcar("nonexistent_outcar.dat", e_fermi));
}

// ─── EIGENVAL error cases ─────────────────────────────────────────────────────

TEST(EigenvalIO, ParseEigenvalEmptyFile) {
    const std::string path = kTmpDir + "/test_io_eigenval_empty_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

// ─── KPOINTS error cases ──────────────────────────────────────────────────────

TEST(EigenvalIO, ParseKpointsEmptyFile) {
    const std::string path = kTmpDir + "/test_io_kpoints_empty_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fclose(f);

    BZPath bz_path;
    EXPECT_FALSE(parseKpoints(path, bz_path));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsInvalidMode) {
    const std::string path = kTmpDir + "/test_io_kpoints_badmode_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("Invalid Mode Test\n", f);
    std::fputs("20\n", f);
    std::fputs("Cartesian\n", f);  // Expected "Line-mode"
    std::fclose(f);

    BZPath bz_path;
    EXPECT_FALSE(parseKpoints(path, bz_path));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsOddLinesBoundaryViolation) {
    const std::string path = kTmpDir + "/test_io_kpoints_oddlines_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("Odd lines check\n", f);
    std::fputs("15\n", f);
    std::fputs("Line-mode\n", f);
    std::fputs("Reciprocal\n", f);
    std::fputs("  0.0 0.0 0.0 ! Gamma\n", f);  // 1 coordinate line (odd)
    std::fclose(f);

    BZPath bz_path;
    EXPECT_FALSE(parseKpoints(path, bz_path));
    std::remove(path.c_str());
}

// ─── DOSCAR error cases ───────────────────────────────────────────────────────

TEST(EigenvalIO, ParseFromDoscarEmptyFile) {
    const std::string path = kTmpDir + "/test_io_doscar_empty_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromDoscar(path, e_fermi));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseFromDoscarShortHeader) {
    const std::string path = kTmpDir + "/test_io_doscar_short_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("line1\n", f);
    std::fputs("line2\n", f);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromDoscar(path, e_fermi));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseFromDoscarMalformedLine6) {
    const std::string path = kTmpDir + "/test_io_doscar_badl6_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("l1\nl2\nl3\nl4\nl5\n", f);
    std::fputs("only_one_value\n", f);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromDoscar(path, e_fermi));
    std::remove(path.c_str());
}

// ─── OUTCAR error cases ───────────────────────────────────────────────────────

TEST(EigenvalIO, ParseFromOutcarNoEFermi) {
    const std::string path = kTmpDir + "/test_io_outcar_nofermi_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs(" VASP output\n", f);
    std::fputs(" no fermi level here...\n", f);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromOutcar(path, e_fermi));
    std::remove(path.c_str());
}

// ─── LOGGER validation ────────────────────────────────────────────────────────

TEST(EigenvalIO, WriteKpointsLogInvalidPath) {
    BZPath bz_path;
    EigenvalData data;
    double final_dist = 0.0;
    // Attempting to write to a non-existent subdirectory path
    EXPECT_FALSE(writeKpointsLog("/invalid_dir_xyz/kpoints.log", bz_path, data, final_dist));
}
