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
    int kpts = 0;
    EXPECT_FALSE(parseKpoints("nonexistent_kpoints.dat", kpts));
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

TEST(EigenvalIO, ParseEigenvalTruncatedBeforeLine6) {
    const std::string path = kTmpDir + "/test_io_eigenval_trunc_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("header1\n", f);
    std::fputs("header2\n", f);
    std::fputs("header3\n", f);
    // Only 3 lines, parser needs 5 header + 1 data
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalMalformedLine6) {
    const std::string path = kTmpDir + "/test_io_eigenval_badl6_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("not integers here\n", f);
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalZeroKPoints) {
    const std::string path = kTmpDir + "/test_io_eigenval_zerokpt_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 0 3\n", f);  // nkpts=0
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalNegativeKPoints) {
    const std::string path = kTmpDir + "/test_io_eigenval_negkpt_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 -1 3\n", f);  // nkpts=-1
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalZeroBands) {
    const std::string path = kTmpDir + "/test_io_eigenval_zeroband_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 1 0\n", f);  // nbands=0
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalMissingKPointSeparator) {
    // File ends immediately after line 6 — no blank separator follows
    const std::string path = kTmpDir + "/test_io_eigenval_nosep_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 1 1\n", f);
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalMissingKPointCoords) {
    const std::string path = kTmpDir + "/test_io_eigenval_nocoords_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 1 1\n", f);
    std::fputs("\n", f);  // blank separator present, but no coords line
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalMalformedKPointLine) {
    const std::string path = kTmpDir + "/test_io_eigenval_badkpt_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 1 2\n", f);
    std::fputs("\n", f);
    std::fputs("not_a_number x y z\n", f);  // kx can't be parsed as double
    std::fputs("1 -5.0\n", f);
    std::fputs("2 -3.0\n", f);
    std::fclose(f);

    EigenvalData data;
    EXPECT_FALSE(parseEigenval(path, data));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseEigenvalMissingBandData) {
    const std::string path = kTmpDir + "/test_io_eigenval_noband_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("h1\nh2\nh3\nh4\nh5\n", f);
    std::fputs("1 1 2\n", f);
    std::fputs("\n", f);
    std::fputs("0.0 0.0 0.0 0.5\n", f);
    // Band lines missing — file ends here
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

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsMissingLine3) {
    // Only 2 lines; parser needs line 3 to confirm line mode
    const std::string path = kTmpDir + "/test_io_kpoints_trunc_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("K-Path\n", f);
    std::fputs("10\n", f);
    std::fclose(f);

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsNotLineMode) {
    const std::string path = kTmpDir + "/test_io_kpoints_notline_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("K-Path\n", f);
    std::fputs("10\n", f);
    std::fputs("Monkhorst-Pack\n", f);  // doesn't start with 'L'
    std::fclose(f);

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsPerSegmentOne) {
    // kpts_per_seg must be > 1
    const std::string path = kTmpDir + "/test_io_kpoints_one_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("K-Path\n", f);
    std::fputs("1\n", f);
    std::fputs("Line-mode\n", f);
    std::fclose(f);

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsPerSegmentZero) {
    const std::string path = kTmpDir + "/test_io_kpoints_zero_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("K-Path\n", f);
    std::fputs("0\n", f);
    std::fputs("Line-mode\n", f);
    std::fclose(f);

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseKpointsMalformedLine2) {
    const std::string path = kTmpDir + "/test_io_kpoints_badl2_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fputs("K-Path\n", f);
    std::fputs("not_an_integer\n", f);
    std::fputs("Line-mode\n", f);
    std::fclose(f);

    int kpts = 0;
    EXPECT_FALSE(parseKpoints(path, kpts));
    std::remove(path.c_str());
}

// ─── DOSCAR error cases ───────────────────────────────────────────────────────

TEST(EigenvalIO, ParseFromDoscarTruncated) {
    // Fewer than 6 lines — parser needs lines 1-5 (header) + line 6 (data)
    const std::string path = kTmpDir + "/test_io_doscar_trunc_tmp.dat";
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
    // Line 6 needs at least 4 numeric values: Emax Emin NEDOS Efermi
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
    std::fputs(" no fermi level here\n", f);
    std::fputs(" just regular output\n", f);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromOutcar(path, e_fermi));
    std::remove(path.c_str());
}

TEST(EigenvalIO, ParseFromOutcarEmptyFile) {
    const std::string path = kTmpDir + "/test_io_outcar_empty_tmp.dat";
    FILE* f = std::fopen(path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    std::fclose(f);

    double e_fermi = 0.0;
    EXPECT_FALSE(parseFromOutcar(path, e_fermi));
    std::remove(path.c_str());
}
