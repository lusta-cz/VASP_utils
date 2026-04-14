// Tests for getBravaisKPath against Setyawan-Curtarolo (2010) tables,
// verified numerically against SeeK-path / HPKOT reference.

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>

#include "kpath.h"
#include "poscar_file.h"

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static constexpr double kTol = 1e-9;

// Build a minimal SpglibDataset with only the fields used by detectBravais.
static SpglibDataset makeDataset(int sg, const char* intl) {
    SpglibDataset ds{};
    ds.spacegroup_number = sg;
    std::strncpy(ds.international_symbol, intl, sizeof(ds.international_symbol) - 1);
    return ds;
}

// Build a POSCAR whose lattice has the given row vectors (row i = lattice[i][*]).
static POSCAR makePOSCAR(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21,
                         double a22) {
    POSCAR p{};
    p.is_direct = true;
    p.total_atoms = 1;
    p.lattice[0][0] = a00;
    p.lattice[0][1] = a01;
    p.lattice[0][2] = a02;
    p.lattice[1][0] = a10;
    p.lattice[1][1] = a11;
    p.lattice[1][2] = a12;
    p.lattice[2][0] = a20;
    p.lattice[2][1] = a21;
    p.lattice[2][2] = a22;
    return p;
}

// Find a KPoint by label; returns nullptr if absent.
static const KPoint* findPt(const KPath& kp, const std::string& label) {
    for (const auto& pt : kp.points)
        if (pt.label == label)
            return &pt;
    return nullptr;
}

// Assert that a named point has the expected coordinates.
#define EXPECT_KPT(kpath, label, ex, ey, ez)         \
    do {                                             \
        const KPoint* p_ = findPt((kpath), (label)); \
        ASSERT_NE(p_, nullptr) << "missing " label;  \
        EXPECT_NEAR(p_->kx, (ex), kTol);             \
        EXPECT_NEAR(p_->ky, (ey), kTol);             \
        EXPECT_NEAR(p_->kz, (ez), kTol);             \
    } while (0)

// ─────────────────────────────────────────────────────────────────────────────
// cP  (Setyawan Table 2)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, cP_BravaisLabel) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(221, "Pm-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "cP");
}

TEST(KPath, cP_Points) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(221, "Pm-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "M", 0.5, 0.5, 0.0);
    EXPECT_KPT(*kp, "R", 0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "X", 0.0, 0.5, 0.0);
}

TEST(KPath, cP_Segments) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(221, "Pm-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->segments.size(), 7u);
}

// ─────────────────────────────────────────────────────────────────────────────
// cF  (Setyawan Table 3)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, cF_BravaisLabel) {
    auto conv = makePOSCAR(5, 0, 0, 0, 5, 0, 0, 0, 5);
    auto ds = makeDataset(225, "Fm-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "cF");
}

TEST(KPath, cF_Points) {
    auto conv = makePOSCAR(5, 0, 0, 0, 5, 0, 0, 0, 5);
    auto ds = makeDataset(225, "Fm-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "K", 0.375, 0.375, 0.75);
    EXPECT_KPT(*kp, "L", 0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "U", 0.625, 0.25, 0.625);
    EXPECT_KPT(*kp, "W", 0.5, 0.25, 0.75);
    EXPECT_KPT(*kp, "X", 0.5, 0.0, 0.5);
}

// ─────────────────────────────────────────────────────────────────────────────
// cI  (Setyawan Table 4)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, cI_Points) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(229, "Im-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "cI");
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "H", 0.5, -0.5, 0.5);
    EXPECT_KPT(*kp, "N", 0.0, 0.0, 0.5);
    EXPECT_KPT(*kp, "P", 0.25, 0.25, 0.25);
}

// ─────────────────────────────────────────────────────────────────────────────
// tP  (Setyawan Table 5)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, tP_Points) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 6);  // a=b=4, c=6
    auto ds = makeDataset(99, "P4mm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "tP");
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "A", 0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "M", 0.5, 0.5, 0.0);
    EXPECT_KPT(*kp, "R", 0.0, 0.5, 0.5);
    EXPECT_KPT(*kp, "X", 0.0, 0.5, 0.0);
    EXPECT_KPT(*kp, "Z", 0.0, 0.0, 0.5);
}

// ─────────────────────────────────────────────────────────────────────────────
// tI1 (c <= a)  verified vs SeeK-path
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, tI1_BravaisLabel) {
    // a=4, c=3 → c < a → tI1
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 3);
    auto ds = makeDataset(139, "I4/mmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "tI1");
}

TEST(KPath, tI1_ZPoint) {
    // H = (1 + c²/a²)/4 = (1 + 9/16)/4 = 25/64
    const double H = (1.0 + 9.0 / 16.0) / 4.0;  // 0.390625
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 3);
    auto ds = makeDataset(139, "I4/mmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_KPT(*kp, "Z", H, H, -H);
    EXPECT_KPT(*kp, "Z1", -H, 1.0 - H, H);
    EXPECT_KPT(*kp, "M", -0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "N", 0.0, 0.5, 0.0);
    EXPECT_KPT(*kp, "P", 0.25, 0.25, 0.25);
    EXPECT_KPT(*kp, "X", 0.0, 0.0, 0.5);
}

TEST(KPath, tI1_Equal_ca) {
    // c == a edge case → tI1
    const double H = 0.5;  // (1+1)/4 = 0.5
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(139, "I4/mmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "tI1");
    EXPECT_KPT(*kp, "Z", H, H, -H);
}

// ─────────────────────────────────────────────────────────────────────────────
// tI2 (c > a)  verified vs SeeK-path (Indium: a≈3.25, c≈4.95)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, tI2_BravaisLabel) {
    // a=3, c=4 → c > a → tI2
    auto conv = makePOSCAR(3, 0, 0, 0, 3, 0, 0, 0, 4);
    auto ds = makeDataset(139, "I4/mmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "tI2");
}

TEST(KPath, tI2_Points) {
    // a=3, c=4: H=(1+9/16)/4=25/64, Z=9/32
    const double H = (1.0 + 9.0 / 16.0) / 4.0;  // 0.390625
    const double Z = 9.0 / 32.0;                // 0.28125
    auto conv = makePOSCAR(3, 0, 0, 0, 3, 0, 0, 0, 4);
    auto ds = makeDataset(139, "I4/mmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_KPT(*kp, "Z", 0.5, 0.5, -0.5);
    EXPECT_KPT(*kp, "N", 0.0, 0.5, 0.0);
    EXPECT_KPT(*kp, "P", 0.25, 0.25, 0.25);
    EXPECT_KPT(*kp, "X", 0.0, 0.0, 0.5);
    EXPECT_KPT(*kp, "Sigma", -H, H, H);
    EXPECT_KPT(*kp, "Sigma1", H, 1.0 - H, -H);
    EXPECT_KPT(*kp, "R", -Z, Z, 0.5);
    EXPECT_KPT(*kp, "Y1", 0.5, 0.5, -Z);
}

// ─────────────────────────────────────────────────────────────────────────────
// oP  (Setyawan Table 7)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, oP_Points) {
    auto conv = makePOSCAR(3, 0, 0, 0, 4, 0, 0, 0, 5);
    auto ds = makeDataset(47, "Pmmm");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "oP");
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "R", 0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "S", 0.5, 0.5, 0.0);
    EXPECT_KPT(*kp, "T", 0.0, 0.5, 0.5);
    EXPECT_KPT(*kp, "U", 0.5, 0.0, 0.5);
    EXPECT_KPT(*kp, "X", 0.5, 0.0, 0.0);
    EXPECT_KPT(*kp, "Y", 0.0, 0.5, 0.0);
    EXPECT_KPT(*kp, "Z", 0.0, 0.0, 0.5);
}

// ─────────────────────────────────────────────────────────────────────────────
// hP  (Setyawan Table 10)
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, hP_Points) {
    // Hexagonal conventional: a1=(a,0,0), a2=(-a/2, a*sqrt(3)/2, 0), a3=(0,0,c)
    // Store as row vectors in POSCAR (lattice[i] = i-th vector)
    const double a = 3.21, c = 5.21;
    auto conv = makePOSCAR(a, 0, 0, -a / 2, a * std::sqrt(3.0) / 2, 0, 0, 0, c);
    auto ds = makeDataset(194, "P6_3/mmc");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "hP");
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "A", 0.0, 0.0, 0.5);
    EXPECT_KPT(*kp, "H", 1.0 / 3.0, 1.0 / 3.0, 0.5);
    EXPECT_KPT(*kp, "K", 1.0 / 3.0, 1.0 / 3.0, 0.0);
    EXPECT_KPT(*kp, "L", 0.5, 0.0, 0.5);
    EXPECT_KPT(*kp, "M", 0.5, 0.0, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// hR1 (cos α_r > 0)  — Bi-like: a_h=4.546, c_h=11.862
// Verified vs SeeK-path: all shared coordinates are identical.
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, hR1_BravaisLabel) {
    const double a_h = 4.546, c_h = 11.862;
    auto conv = makePOSCAR(a_h, 0, 0, -a_h / 2, a_h * std::sqrt(3.0) / 2, 0, 0, 0, c_h);
    auto ds = makeDataset(166, "R-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "hR1");
}

TEST(KPath, hR1_ParametricPoints) {
    const double a_h = 4.546, c_h = 11.862;
    const double a2 = a_h * a_h, c2 = c_h * c_h;
    const double cos_alpha = (2 * c2 - 3 * a2) / (2 * c2 + 6 * a2);
    const double eta = (1.0 + 4.0 * cos_alpha) / (2.0 + 4.0 * cos_alpha);
    const double nu = 0.75 - eta / 2.0;

    auto conv = makePOSCAR(a_h, 0, 0, -a_h / 2, a_h * std::sqrt(3.0) / 2, 0, 0, 0, c_h);
    auto ds = makeDataset(166, "R-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());

    // Fixed points
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "Z", 0.5, 0.5, 0.5);
    EXPECT_KPT(*kp, "L", 0.5, 0.0, 0.0);
    EXPECT_KPT(*kp, "L1", 0.0, 0.0, -0.5);
    EXPECT_KPT(*kp, "F", 0.5, 0.5, 0.0);
    // eta/nu-dependent points (match SeeK-path exactly)
    EXPECT_KPT(*kp, "X", nu, 0.0, -nu);
    EXPECT_KPT(*kp, "Q", 1.0 - nu, nu, 0.0);
    EXPECT_KPT(*kp, "B", eta, 0.5, 1.0 - eta);
    EXPECT_KPT(*kp, "P", eta, nu, nu);
    EXPECT_KPT(*kp, "P2", nu, nu, eta - 1.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// hR2 (cos α_r < 0)  — a_h == c_h gives cos = -1/8
// Verified vs SeeK-path: all shared coordinates are identical.
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, hR2_BravaisLabel) {
    const double a_h = 4.0, c_h = 4.0;
    auto conv = makePOSCAR(a_h, 0, 0, -a_h / 2, a_h * std::sqrt(3.0) / 2, 0, 0, 0, c_h);
    auto ds = makeDataset(166, "R-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "hR2");
}

TEST(KPath, hR2_ParametricPoints) {
    const double a_h = 4.0, c_h = 4.0;
    const double a2 = a_h * a_h, c2 = c_h * c_h;
    const double cos_alpha = (2 * c2 - 3 * a2) / (2 * c2 + 6 * a2);  // -0.125
    const double eta = (1.0 + cos_alpha) / (2.0 * (1.0 - cos_alpha));
    const double nu = 0.75 - eta / 2.0;

    auto conv = makePOSCAR(a_h, 0, 0, -a_h / 2, a_h * std::sqrt(3.0) / 2, 0, 0, 0, c_h);
    auto ds = makeDataset(166, "R-3m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());

    // Fixed points
    EXPECT_KPT(*kp, "Gamma", 0.0, 0.0, 0.0);
    EXPECT_KPT(*kp, "Z", 0.5, -0.5, 0.5);
    EXPECT_KPT(*kp, "L", 0.5, 0.0, 0.0);
    EXPECT_KPT(*kp, "F", 0.5, -0.5, 0.0);
    // eta/nu-dependent points (match SeeK-path exactly)
    EXPECT_KPT(*kp, "Q", eta, eta, eta);
    EXPECT_KPT(*kp, "Q1", 1.0 - eta, -eta, -eta);
    EXPECT_KPT(*kp, "P", 1.0 - nu, -nu, 1.0 - nu);
    EXPECT_KPT(*kp, "P1", nu, nu - 1.0, nu - 1.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Unknown / unsupported type returns nullopt
// ─────────────────────────────────────────────────────────────────────────────

TEST(KPath, UnknownSG_ReturnsNullopt) {
    auto conv = makePOSCAR(4, 0, 0, 0, 4, 0, 0, 0, 4);
    auto ds = makeDataset(0, "???");  // invalid SG — not in any range
    auto kp = getBravaisKPath(conv, ds);
    EXPECT_FALSE(kp.has_value());
}

// SG 10 P2/m is monoclinic P — now fully supported.
TEST(KPath, mP_NowSupported) {
    const double b = 4.0, c = 5.0;
    const double alpha = 70.0 * M_PI / 180.0;
    // b vector along y, c vector at angle alpha from b
    auto conv = makePOSCAR(3, 0, 0, 0, b, 0, 0, c * std::cos(alpha), c * std::sin(alpha));
    auto ds = makeDataset(10, "P2/m");
    auto kp = getBravaisKPath(conv, ds);
    ASSERT_TRUE(kp.has_value());
    EXPECT_EQ(kp->bravais_label, "mP");
}
