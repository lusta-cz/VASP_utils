// k-path tables based on Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010).
// K-point coordinates are in fractional units of the PRIMITIVE reciprocal lattice.
// Use the output POSCAR_prim (primitive=1 from spglib) as input for VASP.

#include "kpath.h"

#include <cmath>
#include <string>

#include "poscar_file.h"

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Bravais lattice classification
// ─────────────────────────────────────────────────────────────────────────────

enum class BravaisType { cP, cF, cI, tP, tI, oP, oF, oI, oC, hP, hR, mP, mC, aP, Unknown };

static BravaisType detectBravais(int sg, const char* intl_sym) {
    // Extract centering letter: first non-space character of international symbol.
    char centering = '\0';
    for (int i = 0; intl_sym[i] != '\0'; ++i) {
        if (intl_sym[i] != ' ') {
            centering = intl_sym[i];
            break;
        }
    }

    if (sg >= 195 && sg <= 230) {
        if (centering == 'P')
            return BravaisType::cP;
        if (centering == 'F')
            return BravaisType::cF;
        if (centering == 'I')
            return BravaisType::cI;
    } else if (sg >= 75 && sg <= 142) {
        if (centering == 'P')
            return BravaisType::tP;
        if (centering == 'I')
            return BravaisType::tI;
    } else if (sg >= 16 && sg <= 74) {
        if (centering == 'P')
            return BravaisType::oP;
        if (centering == 'F')
            return BravaisType::oF;
        if (centering == 'I')
            return BravaisType::oI;
        if (centering == 'C' || centering == 'A' || centering == 'B')
            return BravaisType::oC;
    } else if (sg >= 168 && sg <= 194) {
        return BravaisType::hP;
    } else if (sg >= 143 && sg <= 167) {
        if (centering == 'R')
            return BravaisType::hR;
        return BravaisType::hP;  // P centering in this range is also hexagonal
    } else if (sg >= 3 && sg <= 15) {
        if (centering == 'P')
            return BravaisType::mP;
        if (centering == 'C' || centering == 'A' || centering == 'B' || centering == 'I')
            return BravaisType::mC;
    } else if (sg >= 1 && sg <= 2) {
        return BravaisType::aP;
    }
    return BravaisType::Unknown;
}

// ─────────────────────────────────────────────────────────────────────────────
// Lattice parameter helpers
// ─────────────────────────────────────────────────────────────────────────────

static double vecnorm(const double v[3]) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double vecdot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void veccross(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// Angle between two vectors in radians.
static double vecangle(const double a[3], const double b[3]) {
    double cosval = vecdot(a, b) / (vecnorm(a) * vecnorm(b));
    // Clamp to [-1, 1] to guard against floating-point noise.
    if (cosval > 1.0)
        cosval = 1.0;
    if (cosval < -1.0)
        cosval = -1.0;
    return std::acos(cosval);
}

// Compute reciprocal lattice vectors (without 2π factor) from POSCAR row vectors.
// b1 = (a2 × a3) / V,  b2 = (a3 × a1) / V,  b3 = (a1 × a2) / V
static void reciprocal_lattice(const POSCAR& conv, double b1[3], double b2[3], double b3[3]) {
    const double* a1 = conv.lattice[0];
    const double* a2 = conv.lattice[1];
    const double* a3 = conv.lattice[2];

    double a2xa3[3], a3xa1[3], a1xa2[3];
    veccross(a2, a3, a2xa3);
    veccross(a3, a1, a3xa1);
    veccross(a1, a2, a1xa2);

    const double V = vecdot(a1, a2xa3);
    for (int i = 0; i < 3; ++i) {
        b1[i] = a2xa3[i] / V;
        b2[i] = a3xa1[i] / V;
        b3[i] = a1xa2[i] / V;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// K-path builders  (Setyawan-Curtarolo Tables 2-14)
// ─────────────────────────────────────────────────────────────────────────────

static KPath kpath_cP() {
    KPath p;
    p.bravais_label = "cP";
    p.points = {{"Gamma", 0, 0, 0}, {"M", 0.5, 0.5, 0}, {"R", 0.5, 0.5, 0.5}, {"X", 0, 0.5, 0}};
    p.segments = {{"Gamma", "X"}, {"X", "M"}, {"M", "Gamma"}, {"Gamma", "R"}, {"R", "X"}, {"X", "M"}, {"M", "R"}};
    return p;
}

static KPath kpath_cF() {
    KPath p;
    p.bravais_label = "cF";
    p.points = {{"Gamma", 0, 0, 0}, {"K", 0.375, 0.375, 0.75}, {"L", 0.5, 0.5, 0.5}, {"U", 0.625, 0.25, 0.625},
                {"W", 0.5, 0.25, 0.75}, {"X", 0.5, 0, 0.5}};
    p.segments = {{"Gamma", "X"}, {"X", "W"}, {"W", "K"}, {"K", "Gamma"}, {"Gamma", "L"}, {"L", "U"},
                  {"U", "W"}, {"W", "L"}, {"L", "K"}, {"K", "U"}, {"U", "X"}};
    return p;
}

static KPath kpath_cI() {
    KPath p;
    p.bravais_label = "cI";
    p.points = {{"Gamma", 0, 0, 0}, {"H", 0.5, -0.5, 0.5}, {"N", 0, 0, 0.5}, {"P", 0.25, 0.25, 0.25}};
    p.segments = {{"Gamma", "H"}, {"H", "N"}, {"N", "Gamma"}, {"Gamma", "P"}, {"P", "H"}, {"P", "N"}};
    return p;
}

static KPath kpath_tP() {
    KPath p;
    p.bravais_label = "tP";
    p.points = {{"Gamma", 0, 0, 0}, {"A", 0.5, 0.5, 0.5}, {"M", 0.5, 0.5, 0},
                {"R", 0, 0.5, 0.5}, {"X", 0, 0.5, 0},     {"Z", 0, 0, 0.5}};
    p.segments = {{"Gamma", "X"}, {"X", "M"}, {"M", "Gamma"}, {"Gamma", "Z"}, {"Z", "R"},
                  {"R", "A"}, {"A", "Z"}, {"X", "R"}, {"M", "A"}};
    return p;
}

// tI: two subcases based on c vs a of the conventional tetragonal cell.
// Follows HPKOT (Hinuma et al. 2017) / SeeK-path convention.
// a = |lattice[0]|, c = |lattice[2]|.
//
// tI1 (c <= a):  H = (1 + c²/a²) / 4
// tI2 (c  > a):  H = (1 + a²/c²) / 4,   Z = a² / (2c²)
static KPath kpath_tI(const POSCAR& conv) {
    const double a = vecnorm(conv.lattice[0]);
    const double c = vecnorm(conv.lattice[2]);
    const double a2 = a * a;
    const double c2 = c * c;

    KPath p;
    if (c <= a) {
        // tI1
        const double H = (1.0 + c2 / a2) / 4.0;
        p.bravais_label = "tI1";
        p.points = {{"Gamma", 0, 0, 0},   {"M", -0.5, 0.5, 0.5}, {"N", 0, 0.5, 0},      {"P", 0.25, 0.25, 0.25},
                    {"X", 0, 0, 0.5}, {"Z", H, H, -H},       {"Z1", -H, 1.0 - H, H}};
        p.segments = {{"Gamma", "X"}, {"X", "M"}, {"M", "Gamma"}, {"Gamma", "Z"}, {"Z1", "M"}, {"X", "P"}, {"P", "N"}, {"N", "Gamma"}};
    } else {
        // tI2
        const double H = (1.0 + a2 / c2) / 4.0;
        const double Z = a2 / (2.0 * c2);
        p.bravais_label = "tI2";
        p.points = {{"Gamma", 0, 0, 0},          {"Z", 0.5, 0.5, -0.5}, {"N", 0, 0.5, 0},
                    {"P", 0.25, 0.25, 0.25}, {"X", 0, 0, 0.5},      {"Sigma", -H, H, H},
                    {"Sigma1", H, 1.0 - H, -H},   {"Y", -Z, Z, 0.5},     {"Y1", 0.5, 0.5, -Z}};
        p.segments = {{"Gamma", "X"}, {"X", "P"},  {"P", "N"}, {"N", "Gamma"}, {"Gamma", "Z"},
                      {"Z", "Sigma"}, {"Sigma1", "Gamma"}, {"X", "Y"}, {"Y1", "Z"}};
    }
    return p;
}

static KPath kpath_oP() {
    KPath p;
    p.bravais_label = "oP";
    p.points = {{"Gamma", 0, 0, 0},     {"R", 0.5, 0.5, 0.5}, {"S", 0.5, 0.5, 0}, {"T", 0, 0.5, 0.5},
                {"U", 0.5, 0, 0.5}, {"X", 0.5, 0, 0},     {"Y", 0, 0.5, 0},   {"Z", 0, 0, 0.5}};
    p.segments = {{"Gamma", "X"}, {"X", "S"}, {"S", "Y"}, {"Y", "Gamma"}, {"Gamma", "Z"}, {"Z", "U"},
                  {"U", "R"}, {"R", "T"}, {"T", "Z"}, {"Y", "T"}, {"U", "X"}, {"S", "R"}};
    return p;
}

// oF: face-centred orthorhombic — three subcases based on 1/a² vs 1/b²+1/c².
static KPath kpath_oF(const POSCAR& conv) {
    const double a = vecnorm(conv.lattice[0]);
    const double b = vecnorm(conv.lattice[1]);
    const double c = vecnorm(conv.lattice[2]);
    const double a2 = a * a, b2 = b * b, c2 = c * c;

    KPath p;
    if (1.0 / a2 >= 1.0 / b2 + 1.0 / c2) {
        // oF1
        const double zeta = (1.0 + a2 / b2 - a2 / c2) / 4.0;
        const double eta = (1.0 + a2 / b2 + a2 / c2) / 4.0;
        p.bravais_label = "oF1";
        p.points = {{"Gamma", 0, 0, 0},
                    {"A", 0.5, 0.5 + zeta, zeta},
                    {"A1", 0.5, 0.5 - zeta, 1.0 - zeta},
                    {"L", 0.5, 0.5, 0.5},
                    {"T", 1.0, 0.5, 0.5},
                    {"X", 0.0, eta, eta},
                    {"X1", 1.0, 1.0 - eta, 1.0 - eta},
                    {"Y", 0.5, 0.0, 0.5},
                    {"Z", 0.5, 0.5, 0.0}};
        p.segments = {{"Gamma", "Y"},  {"Y", "T"},  {"T", "Z"}, {"Z", "Gamma"}, {"Gamma", "X"}, {"X", "A1"},
                      {"A1", "Y"}, {"T", "X1"}, {"X", "A"}, {"A", "Z"}, {"L", "Gamma"}};
    } else (1.0 / a2 < 1.0 / b2 + 1.0 / c2) {
        // oF2
        const double phi = (1.0 + c2 / b2 - c2 / a2) / 4.0;
        const double eta = (1.0 + a2 / b2 - a2 / c2) / 4.0;
        const double delta = (1.0 + b2 / a2 - b2 / c2) / 4.0;
        p.bravais_label = "oF2";
        p.points = {{"Gamma", 0, 0, 0},
                    {"C", 0.5, 0.5 - eta, 1.0 - eta},
                    {"C1", 0.5, 0.5 + eta, eta},
                    {"D", 0.5 - delta, 0.5, 1.0 - delta},
                    {"D1", 0.5 + delta, 0.5, delta},
                    {"L", 0.5, 0.5, 0.5},
                    {"H", 1.0 - phi, 0.5 - phi, 0.5},
                    {"H1", phi, 0.5 + phi, 0.5},
                    {"X", 0.0, 0.5, 0.5},
                    {"Y", 0.5, 0.0, 0.5},
                    {"Z", 0.5, 0.5, 0.0}};
        p.segments = {{"Gamma", "Y"},  {"Y", "C"}, {"C", "D"},  {"D", "X"},  {"X", "Gamma"}, {"Gamma", "Z"}, {"Z", "D1"},
                      {"D1", "H"}, {"H", "C"}, {"C1", "Z"}, {"X", "H1"}, {"H", "Y"}, {"L", "Gamma"}};
    } else {
        // oF3 (degenerate: 1/a² == 1/b² + 1/c²)
        const double zeta = (1.0 + a2 / b2 - a2 / c2) / 4.0;
        const double eta = (1.0 + a2 / b2 + a2 / c2) / 4.0;
        p.bravais_label = "oF3";
        p.points = {{"Gamma", 0, 0, 0},
                    {"A", 0.5, 0.5 + zeta, zeta},
                    {"A1", 0.5, 0.5 - zeta, 1.0 - zeta},
                    {"L", 0.5, 0.5, 0.5},
                    {"T", 1.0, 0.5, 0.5},
                    {"X", 0.0, eta, eta},
                    {"X1", 1.0, 1.0 - eta, 1.0 - eta},
                    {"Y", 0.5, 0.0, 0.5},
                    {"Z", 0.5, 0.5, 0.0}};
        p.segments = {{"Gamma", "Y"},  {"Y", "T"},  {"T", "Z"}, {"Z", "Gamma"}, {"Gamma", "X"},
                      {"X", "A1"}, {"A1", "Y"}, {"X", "A"}, {"A", "Z"}, {"L", "Gamma"}};
    }
    return p;
}

// oI: body-centred orthorhombic.
static KPath kpath_oI(const POSCAR& conv) {
    const double a = vecnorm(conv.lattice[0]);
    const double b = vecnorm(conv.lattice[1]);
    const double c = vecnorm(conv.lattice[2]);
    const double a2 = a * a, b2 = b * b, c2 = c * c;
    const double zeta = (1.0 + a2 / c2) / 4.0;
    const double eta = (1.0 + b2 / c2) / 4.0;
    const double delta = (b2 - a2) / (4.0 * c2);
    const double mu = (a2 + b2) / (4.0 * c2);

    KPath p;
    p.bravais_label = "oI";
    p.points = {{"Gamma", 0, 0, 0},
                {"L", -mu, mu, 0.5 - delta},
                {"L1", mu, -mu, 0.5 + delta},
                {"L2", 0.5 - delta, 0.5 + delta, -mu},
                {"R", 0.0, 0.5, 0.0},
                {"S", 0.5, 0.0, 0.0},
                {"T", 0.0, 0.0, 0.5},
                {"W", 0.25, 0.25, 0.25},
                {"X", -zeta, zeta, zeta},
                {"X1", zeta, 1.0 - zeta, -zeta},
                {"Y", eta, -eta, eta},
                {"Y1", 1.0 - eta, eta, -eta},
                {"Z", 0.5, 0.5, -0.5}};
    p.segments = {{"Gamma", "X"}, {"X", "L"}, {"L", "T"}, {"T", "W"}, {"W", "R"},  {"R", "X1"}, {"X1", "Z"},
                  {"Z", "Gamma"}, {"Gamma", "Y"}, {"Y", "S"}, {"S", "W"}, {"L1", "Y"}, {"Y1", "Z"}};
    return p;
}

// oC: C-centred (or A-centred) orthorhombic.
static KPath kpath_oC(const POSCAR& conv) {
    const double a = vecnorm(conv.lattice[0]);
    const double b = vecnorm(conv.lattice[1]);
    const double zeta = (1.0 + a * a / (b * b)) / 4.0;

    KPath p;
    p.bravais_label = "oC";
    p.points = {{"Gamma", 0, 0, 0},         {"A", zeta, zeta, 0.5},         {"A1", -zeta, 1.0 - zeta, 0.5},
                {"R", 0.0, 0.5, 0.5},   {"S", 0.0, 0.5, 0.0},           {"T", -0.5, 0.5, 0.5},
                {"X", zeta, zeta, 0.0}, {"X1", -zeta, 1.0 - zeta, 0.0}, {"Y", -0.5, 0.5, 0.0},
                {"Z", 0.0, 0.0, 0.5}};
    p.segments = {{"Gamma", "X"}, {"X", "S"},  {"S", "R"},   {"R", "A"},  {"A", "Z"}, {"Z", "Gamma"},
                  {"Gamma", "Y"}, {"Y", "X1"}, {"X1", "A1"}, {"A1", "T"}, {"T", "Y"}, {"Z", "T"}};
    return p;
}

static KPath kpath_hP() {
    KPath p;
    p.bravais_label = "hP";
    p.points = {{"Gamma", 0, 0, 0},     {"A", 0, 0, 0.5}, {"H", 1.0 / 3.0, 1.0 / 3.0, 0.5}, {"K", 1.0 / 3.0, 1.0 / 3.0, 0},
                {"L", 0.5, 0, 0.5}, {"M", 0.5, 0, 0}};
    p.segments = {{"Gamma", "M"}, {"M", "K"}, {"K", "Gamma"}, {"Gamma", "A"}, {"A", "L"},
                  {"L", "H"}, {"H", "A"}, {"L", "M"}, {"K", "H"}};
    return p;
}

// hR: two subcases based on the rhombohedral angle alpha_r.
// Using the conventional hexagonal description:
//   cos(alpha_r) = (2*c_h^2 - 3*a_h^2) / (2*c_h^2 + 6*a_h^2)
//
// hR1: alpha_r < pi/2  (cos > 0)   Setyawan Table 8
// hR2: alpha_r > pi/2  (cos < 0)   Setyawan Table 9
static KPath kpath_hR(const POSCAR& conv) {
    const double a_h = vecnorm(conv.lattice[0]);
    const double c_h = vecnorm(conv.lattice[2]);
    const double a2 = a_h * a_h;
    const double c2 = c_h * c_h;
    const double cos_alpha = (2.0 * c2 - 3.0 * a2) / (2.0 * c2 + 6.0 * a2);

    KPath p;
    if (cos_alpha > 0.0) {
        // hR1
        const double eta = (1.0 + 4.0 * cos_alpha) / (2.0 + 4.0 * cos_alpha);
        const double nu = 0.75 - eta / 2.0;
        p.bravais_label = "hR1";
        p.points = {{"Gamma", 0, 0, 0},
                    {"B", eta, 0.5, 1.0 - eta},
                    {"B1", 0.5, 1.0 - eta, eta - 1.0},
                    {"F", 0.5, 0.5, 0},
                    {"L", 0.5, 0, 0},
                    {"L1", 0, 0, -0.5},
                    {"P", eta, nu, nu},
                    {"P1", 1.0 - nu, 1.0 - nu, 1.0 - eta},
                    {"P2", nu, nu, eta - 1.0},
                    {"Q", 1.0 - nu, nu, 0},
                    {"X", nu, 0, -nu},
                    {"Z", 0.5, 0.5, 0.5}};
        p.segments = {{"Gamma", "L"}, {"L", "B1"}, {"B", "Z"},  {"Z", "Gamma"}, {"Gamma", "X"},
                      {"Q", "F"}, {"F", "P1"}, {"P1", "Z"}, {"L", "P"}};
    } else {
        // hR2
        const double eta = (1.0 + cos_alpha) / (2.0 * (1.0 - cos_alpha));
        const double nu = 0.75 - eta / 2.0;
        p.bravais_label = "hR2";
        p.points = {{"Gamma", 0, 0, 0},
                    {"F", 0.5, -0.5, 0},
                    {"L", 0.5, 0, 0},
                    {"P", 1.0 - nu, -nu, 1.0 - nu},
                    {"P1", nu, nu - 1.0, nu - 1.0},
                    {"Q", eta, eta, eta},
                    {"Q1", 1.0 - eta, -eta, -eta},
                    {"Z", 0.5, -0.5, 0.5}};
        p.segments = {{"Gamma", "P"},  {"P", "Z"},   {"Z", "Q"},  {"Q", "Gamma"}, {"Gamma", "F"},
                      {"F", "P1"}, {"P1", "Q1"}, {"Q1", "L"}, {"L", "Z"}};
    }
    return p;
}

// mP: primitive monoclinic.
// spglib returns the conventional cell in b-unique setting: lattice[1] is the unique axis.
// S-C formula uses the two non-unique axis lengths and the acute monoclinic angle between them.
// Mapping: S-C "b" <- |lattice[0]|, S-C "c" <- |lattice[2]|,
//          S-C angle <- acute(angle(lattice[0], lattice[2])).
// spglib's beta (a-c angle) is always >= 90°; we force it to acute for the S-C formula.
// Setyawan Table 10.
static KPath kpath_mP(const POSCAR& conv) {
    const double b = vecnorm(conv.lattice[0]);  // S-C "b": first non-unique
    const double c = vecnorm(conv.lattice[2]);  // S-C "c": second non-unique
    // S-C uses the acute monoclinic angle; spglib b-unique gives beta >= 90°.
    const double beta_raw = vecangle(conv.lattice[0], conv.lattice[2]);
    const double beta = (beta_raw > M_PI / 2.0) ? (M_PI - beta_raw) : beta_raw;
    const double eta = (1.0 - b * std::cos(beta) / c) / (2.0 * std::sin(beta) * std::sin(beta));
    const double nu = 0.5 - eta * c * std::cos(beta) / b;

    KPath p;
    p.bravais_label = "mP";
    p.points = {
        {"Gamma", 0, 0, 0},         {"A", 0.5, 0.5, 0.0},      {"C", 0.0, 0.5, 0.5},       {"D", 0.5, 0.0, 0.5},
        {"D1", 0.5, 0.0, -0.5}, {"E", 0.5, 0.5, 0.5},      {"H", 0.0, eta, 1.0 - nu},  {"H1", 0.0, 1.0 - eta, nu},
        {"H2", 0.0, eta, -nu},  {"M", 0.5, eta, 1.0 - nu}, {"M1", 0.5, 1.0 - eta, nu}, {"M2", 0.5, eta, -nu},
        {"X", 0.0, 0.5, 0.0},   {"Y", 0.0, 0.0, 0.5},      {"Y1", 0.0, 0.0, -0.5},     {"Z", 0.5, 0.0, 0.0}};
    p.segments = {{"Gamma", "Y"}, {"Y", "H"},  {"H", "C"}, {"C", "E"}, {"E", "M1"}, {"M1", "A"},
                  {"A", "X"}, {"X", "H1"}, {"M", "D"}, {"D", "Z"}, {"Y", "D"}};
    return p;
}

// mC: C-centred monoclinic — five subcases depending on reciprocal lattice angle kgamma
// and the parameter b*cos(alpha)/c + b²sin²(alpha)/a².
// spglib returns the conventional cell in b-unique setting: lattice[1] is the unique axis.
// S-C formula mapping: S-C "a" (unique) <- |lattice[1]|,
//                      S-C "b" <- |lattice[0]|, S-C "c" <- |lattice[2]|,
//                      S-C angle <- acute(angle(lattice[0], lattice[2])).
// kgamma is the gamma angle of the PRIMITIVE cell's reciprocal lattice.
// Setyawan Tables 11a-11e.
static KPath kpath_mC(const POSCAR& conv) {
    const double a = vecnorm(conv.lattice[1]);  // S-C "a": unique axis
    const double b = vecnorm(conv.lattice[0]);  // S-C "b": first non-unique
    const double c = vecnorm(conv.lattice[2]);  // S-C "c": second non-unique
    // S-C uses the acute monoclinic angle; spglib b-unique gives beta >= 90°.
    const double beta_raw = vecangle(conv.lattice[0], conv.lattice[2]);
    const double beta = (beta_raw > M_PI / 2.0) ? (M_PI - beta_raw) : beta_raw;

    // Compute kgamma from the primitive cell's reciprocal lattice.
    // C-centering: a_p = (a_conv + b_conv)/2, b_p = (-a_conv + b_conv)/2, c_p = c_conv.
    double ap[3], bp_[3], cp_[3];
    for (int i = 0; i < 3; ++i) {
        ap[i] = (conv.lattice[0][i] + conv.lattice[1][i]) * 0.5;
        bp_[i] = (conv.lattice[0][i] - conv.lattice[1][i]) * 0.5;  // (a-b)/2, not (-a+b)/2
        cp_[i] = conv.lattice[2][i];
    }
    double bpxcp[3], cpxap[3];
    veccross(bp_, cp_, bpxcp);
    veccross(cp_, ap, cpxap);
    const double Vp = vecdot(ap, bpxcp);
    double b1p[3], b2p[3];
    for (int i = 0; i < 3; ++i) {
        b1p[i] = bpxcp[i] / Vp;
        b2p[i] = cpxap[i] / Vp;
    }
    const double kgamma_deg = vecangle(b1p, b2p) * 180.0 / M_PI;

    KPath p;

    if (kgamma_deg > 90.0) {
        // mC1
        const double zeta = (2.0 - b * std::cos(beta) / c) / (4.0 * std::sin(beta) * std::sin(beta));
        const double eta = 0.5 + 2.0 * zeta * c * std::cos(beta) / b;
        const double psi = 0.75 - a * a / (4.0 * b * b * std::sin(beta) * std::sin(beta));
        const double phi = psi + (0.75 - psi) * b * std::cos(beta) / c;
        p.bravais_label = "mC1";
        p.points = {{"Gamma", 0, 0, 0},
                    {"N", 0.5, 0.0, 0.0},
                    {"N1", 0.0, -0.5, 0.0},
                    {"F", 1.0 - zeta, 1.0 - zeta, 1.0 - eta},
                    {"F1", zeta, zeta, eta},
                    {"F2", -zeta, -zeta, 1.0 - eta},
                    {"I", phi, 1.0 - phi, 0.5},
                    {"I1", 1.0 - phi, phi - 1.0, 0.5},
                    {"L", 0.5, 0.5, 0.5},
                    {"M", 0.5, 0.0, 0.5},
                    {"X", 1.0 - psi, psi - 1.0, 0.0},
                    {"X1", psi, 1.0 - psi, 0.0},
                    {"X2", psi - 1.0, -psi, 0.0},
                    {"Y", 0.5, 0.5, 0.0},
                    {"Y1", -0.5, -0.5, 0.0},
                    {"Z", 0.0, 0.0, 0.5}};
        p.segments = {{"Gamma", "Y"},  {"Y", "F"},  {"F", "L"}, {"L", "I"}, {"I1", "Z"},
                      {"Z", "F1"}, {"Y", "X1"}, {"X", "Gamma"}, {"Gamma", "N"}, {"M", "G"}};
    } else if (kgamma_deg < 90.0) {
        const double test = b * std::cos(beta) / c + b * b * std::sin(beta) * std::sin(beta) / (a * a);
        if (test < 1.0) {
            // mC3
            const double mu = (1.0 + b * b / (a * a)) / 4.0;
            const double delta = b * c * std::cos(beta) / (2.0 * a * a);
            const double zeta = mu - 0.25 + (1.0 - b * std::cos(beta) / c) / (4.0 * std::sin(beta) * std::sin(beta));
            const double eta = 0.5 + 2.0 * zeta * c * std::cos(beta) / b;
            const double phi = 1.0 + zeta - 2.0 * mu;
            const double psi = eta - 2.0 * delta;
            p.bravais_label = "mC3";
            p.points = {{"Gamma", 0, 0, 0},
                        {"F", 1.0 - phi, 1.0 - phi, 1.0 - psi},
                        {"F1", phi, phi - 1.0, psi},
                        {"F2", 1.0 - phi, -phi, 1.0 - psi},
                        {"H", zeta, zeta, eta},
                        {"H1", 1.0 - zeta, -zeta, 1.0 - eta},
                        {"H2", -zeta, -zeta, 1.0 - eta},
                        {"I", 0.5, -0.5, 0.5},
                        {"M", 0.5, 0.0, 0.5},
                        {"N", 0.5, 0.0, 0.0},
                        {"N1", 0.0, -0.5, 0.0},
                        {"X", 0.5, -0.5, 0.0},
                        {"Y", mu, mu, delta},
                        {"Y1", 1.0 - mu, -mu, -delta},
                        {"Y2", -mu, -mu, -delta},
                        {"Y3", mu, mu - 1.0, delta},
                        {"Z", 0.0, 0.0, 0.5}};
            p.segments = {{"Gamma", "Y"},   {"Y", "F"},  {"F", "H"}, {"H", "Z"}, {"Z", "I"}, {"I", "F1"},
                          {"H1", "Y1"}, {"Y1", "X"}, {"X", "Gamma"}, {"Gamma", "N"}, {"M", "Gamma"}};
        } else if (test > 1.0) {
            // mC5
            const double zeta =
                (b * b / (a * a) + (1.0 - b * std::cos(beta) / c) / (std::sin(beta) * std::sin(beta))) / 4.0;
            const double eta = 0.5 + 2.0 * zeta * c * std::cos(beta) / b;
            const double mu = eta / 2.0 + b * b / (4.0 * a * a) - b * c * std::cos(beta) / (2.0 * a * a);
            const double nu = 2.0 * mu - zeta;
            const double rho = 1.0 - zeta * a * a / (b * b);
            const double omega =
                (4.0 * nu - 1.0 - b * b * std::sin(beta) * std::sin(beta) / (a * a)) * c / (2.0 * b * std::cos(beta));
            const double delta = zeta * c * std::cos(beta) / b + omega / 2.0 - 0.25;
            p.bravais_label = "mC5";
            p.points = {{"Gamma", 0, 0, 0},
                        {"F", nu, nu, omega},
                        {"F1", 1.0 - nu, 1.0 - nu, 1.0 - omega},
                        {"F2", nu, nu - 1.0, omega},
                        {"H", zeta, zeta, eta},
                        {"H1", 1.0 - zeta, -zeta, 1.0 - eta},
                        {"H2", -zeta, -zeta, 1.0 - eta},
                        {"I", rho, 1.0 - rho, 0.5},
                        {"I1", 1.0 - rho, rho - 1.0, 0.5},
                        {"L", 0.5, 0.5, 0.5},
                        {"M", 0.5, 0.0, 0.5},
                        {"N", 0.5, 0.0, 0.0},
                        {"N1", 0.0, -0.5, 0.0},
                        {"X", 0.5, -0.5, 0.0},
                        {"Y", mu, mu, delta},
                        {"Y1", 1.0 - mu, -mu, -delta},
                        {"Y2", -mu, -mu, -delta},
                        {"Y3", mu, mu - 1.0, delta},
                        {"Z", 0.0, 0.0, 0.5}};
            p.segments = {{"Gamma", "Y"},  {"Y", "F"},   {"F", "L"},  {"L", "I"}, {"I1", "Z"}, {"Z", "H"},
                          {"H", "F1"}, {"H1", "Y1"}, {"Y1", "X"}, {"X", "Gamma"}, {"Gamma", "N"},  {"M", "Gamma"}};
        } else {
            // mC4 (degenerate)
            const double mu = (1.0 + b * b / (a * a)) / 4.0;
            const double delta = b * c * std::cos(beta) / (2.0 * a * a);
            const double zeta = mu - 0.25 + (1.0 - b * std::cos(beta) / c) / (4.0 * std::sin(beta) * std::sin(beta));
            const double eta = 0.5 + 2.0 * zeta * c * std::cos(beta) / b;
            const double phi = 1.0 + zeta - 2.0 * mu;
            const double psi = eta - 2.0 * delta;
            p.bravais_label = "mC4";
            p.points = {{"Gamma", 0, 0, 0},
                        {"F", 1.0 - phi, 1.0 - phi, 1.0 - psi},
                        {"F1", phi, phi - 1.0, psi},
                        {"F2", 1.0 - phi, -phi, 1.0 - psi},
                        {"H", zeta, zeta, eta},
                        {"H1", 1.0 - zeta, -zeta, 1.0 - eta},
                        {"H2", -zeta, -zeta, 1.0 - eta},
                        {"I", 0.5, -0.5, 0.5},
                        {"M", 0.5, 0.0, 0.5},
                        {"N", 0.5, 0.0, 0.0},
                        {"N1", 0.0, -0.5, 0.0},
                        {"X", 0.5, -0.5, 0.0},
                        {"Y", mu, mu, delta},
                        {"Y1", 1.0 - mu, -mu, -delta},
                        {"Y2", -mu, -mu, -delta},
                        {"Y3", mu, mu - 1.0, delta},
                        {"Z", 0.0, 0.0, 0.5}};
            p.segments = {{"Gamma", "Y"},   {"Y", "F"},  {"F", "H"}, {"H", "Z"}, {"Z", "I"},
                          {"H1", "Y1"}, {"Y1", "X"}, {"X", "Gamma"}, {"Gamma", "N"}, {"M", "Gamma"}};
        }
    } else {
        // mC2 (kgamma == 90°)
        const double zeta = (2.0 - b * std::cos(beta) / c) / (4.0 * std::sin(beta) * std::sin(beta));
        const double eta = 0.5 + 2.0 * zeta * c * std::cos(beta) / b;
        const double psi = 0.75 - a * a / (4.0 * b * b * std::sin(beta) * std::sin(beta));
        const double phi = psi + (0.75 - psi) * b * std::cos(beta) / c;
        p.bravais_label = "mC2";
        p.points = {{"Gamma", 0, 0, 0},
                    {"N", 0.5, 0.0, 0.0},
                    {"N1", 0.0, -0.5, 0.0},
                    {"F", 1.0 - zeta, 1.0 - zeta, 1.0 - eta},
                    {"F1", zeta, zeta, eta},
                    {"F2", -zeta, -zeta, 1.0 - eta},
                    {"F3", 1.0 - zeta, -zeta, 1.0 - eta},
                    {"I", phi, 1.0 - phi, 0.5},
                    {"I1", 1.0 - phi, phi - 1.0, 0.5},
                    {"L", 0.5, 0.5, 0.5},
                    {"M", 0.5, 0.0, 0.5},
                    {"X", 1.0 - psi, psi - 1.0, 0.0},
                    {"X1", psi, 1.0 - psi, 0.0},
                    {"X2", psi - 1.0, -psi, 0.0},
                    {"Y", 0.5, 0.5, 0.0},
                    {"Y1", -0.5, -0.5, 0.0},
                    {"Z", 0.0, 0.0, 0.5}};
        p.segments = {{"Gamma", "Y"}, {"Y", "F"}, {"F", "L"}, {"L", "I"}, {"I1", "Z"}, {"Z", "F1"}, {"N", "Gamma"}, {"Gamma", "M"}};
    }
    return p;
}

// aP: triclinic — two subcases based on reciprocal lattice angles of the
// input primitive cell.
//
// The caller must pass the primitive cell in the orientation that should be
// used for classification (i.e. the as-read POSCAR, not a spglib
// re-standardized version).  spglib may return a different P-1 basis than the
// one committed to by the user or upstream pymatgen, changing the angles.
//
// TRI1a: kalpha>90°, kbeta>90°, kgamma≥90°  (pymatgen name "TRI1a")
// TRI1b: otherwise                           (pymatgen name "TRI1b")
static KPath kpath_aP(const POSCAR& prim) {
    double a2xa3[3], a3xa1[3], a1xa2[3];
    veccross(prim.lattice[1], prim.lattice[2], a2xa3);
    veccross(prim.lattice[2], prim.lattice[0], a3xa1);
    veccross(prim.lattice[0], prim.lattice[1], a1xa2);
    const double V = vecdot(prim.lattice[0], a2xa3);
    double b1[3], b2[3], b3[3];
    for (int i = 0; i < 3; ++i) {
        b1[i] = a2xa3[i] / V;
        b2[i] = a3xa1[i] / V;
        b3[i] = a1xa2[i] / V;
    }
    const double ka = vecangle(b2, b3) * 180.0 / M_PI;
    const double kb = vecangle(b1, b3) * 180.0 / M_PI;
    const double kg = vecangle(b1, b2) * 180.0 / M_PI;
    const bool is_tri1a = (ka > 90.0 && kb > 90.0 && kg >= 90.0);

    KPath p;
    if (is_tri1a) {
        p.bravais_label = "aP1";
        p.points = {{"Gamma", 0, 0, 0},       {"L", 0.5, 0.5, 0.0}, {"M", 0.0, 0.5, 0.5}, {"N", 0.5, 0.0, 0.5},
                    {"R", 0.5, 0.5, 0.5}, {"X", 0.5, 0.0, 0.0}, {"Y", 0.0, 0.5, 0.0}, {"Z", 0.0, 0.0, 0.5}};
        p.segments = {{"X", "Gamma"}, {"Gamma", "Y"}, {"L", "Gamma"}, {"Gamma", "Z"}, {"N", "Gamma"}, {"Gamma", "M"}, {"R", "Gamma"}};
    } else {
        p.bravais_label = "aP2";
        p.points = {{"G", 0, 0, 0},        {"L", 0.5, -0.5, 0.0}, {"M", 0.0, 0.0, 0.5}, {"N", -0.5, -0.5, 0.5},
                    {"R", 0.0, -0.5, 0.5}, {"X", 0.0, -0.5, 0.0}, {"Y", 0.5, 0.0, 0.0}, {"Z", -0.5, 0.0, 0.5}};
        p.segments = {{"X", "Gamma"}, {"Gamma", "Y"}, {"L", "Gamma"}, {"Gamma", "Z"}, {"N", "Gamma"}, {"Gamma", "M"}, {"R", "Gamma"}};
    }
    return p;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Public interface
// ─────────────────────────────────────────────────────────────────────────────

std::optional<KPath> getBravaisKPath(const POSCAR& std_conv, const POSCAR& std_prim, const SpglibDataset& dataset) {
    const BravaisType bt = detectBravais(dataset.spacegroup_number, dataset.international_symbol);

    switch (bt) {
        case BravaisType::cP:
            return kpath_cP();
        case BravaisType::cF:
            return kpath_cF();
        case BravaisType::cI:
            return kpath_cI();
        case BravaisType::tP:
            return kpath_tP();
        case BravaisType::tI:
            return kpath_tI(std_conv);
        case BravaisType::oP:
            return kpath_oP();
        case BravaisType::oF:
            return kpath_oF(std_conv);
        case BravaisType::oI:
            return kpath_oI(std_conv);
        case BravaisType::oC:
            return kpath_oC(std_conv);
        case BravaisType::hP:
            return kpath_hP();
        case BravaisType::hR:
            return kpath_hR(std_conv);
        case BravaisType::mP:
            return kpath_mP(std_conv);
        case BravaisType::mC:
            return kpath_mC(std_conv);
        case BravaisType::aP:
            return kpath_aP(std_prim);  // primitive orientation matches pymatgen
        default:
            return std::nullopt;
    }
}
