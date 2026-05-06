#include <lapacke.h>

#include <algorithm>
#include <array>
#include <CLI/CLI.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "elastic.h"
#include "error_parse.h"

namespace fs = std::filesystem;

// ── Types ─────────────────────────────────────────────────────────────────────

// Stress in Voigt order: [0]=σ_xx [1]=σ_yy [2]=σ_zz [3]=σ_yz [4]=σ_xz [5]=σ_xy  (kBar)
// Sign: tensile positive (mechanics convention).
using StressVoigt = std::array<double, 6>;

// 6×6 elastic constant matrix, C[i][j] in Voigt notation (0-based).
using Cmat = std::array<std::array<double, 6>, 6>;

// ── VASP output parsers ───────────────────────────────────────────────────────

// Returns the final E0 = energy(sigma->0) (eV) from the last ionic step.
// Line format: "   N F= -.16281084E+03 E0= -.16281084E+03  d E =..."
std::optional<double> readOszicar(const fs::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        return std::nullopt;

    std::optional<double> last_e0;
    std::string line;
    while (std::getline(f, line)) {
        const auto pos = line.find("E0=");
        if (pos == std::string::npos)
            continue;
        std::istringstream iss(line.substr(pos + 3));
        double e0;
        if (iss >> e0)
            last_e0 = e0;
    }
    return last_e0;
}

// Returns the final energy(sigma->0) (eV) from OUTCAR.
// Line format: "  energy  without entropy=  -163.281  energy(sigma->0) =  -163.281"
std::optional<double> readOutcarEnergy(const fs::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        return std::nullopt;

    std::optional<double> last_e;
    std::string line;
    while (std::getline(f, line)) {
        const auto pos = line.find("energy(sigma->0) =");
        if (pos == std::string::npos)
            continue;
        std::istringstream iss(line.substr(pos + 18));
        double e;
        if (iss >> e)
            last_e = e;
    }
    return last_e;
}

// Returns the final stress tensor (kBar) from the last ionic step in OUTCAR.
//
// OUTCAR "in kB" column order: XX YY ZZ XY YZ ZX
// Returned Voigt order:        xx yy zz yz xz xy  (indices 0–5)
//
// Note: OUTCAR header reads "FORCE on cell = -STRESS". The "in kB" line therefore
// gives -σ (compressive positive). We negate to recover the Cauchy stress (tensile positive).
std::optional<StressVoigt> readOutcarStress(const fs::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        return std::nullopt;

    std::optional<StressVoigt> last_stress;
    std::string line;
    while (std::getline(f, line)) {
        if (line.find("in kB") == std::string::npos)
            continue;
        std::istringstream iss(line.substr(line.find("in kB") + 5));
        double xx, yy, zz, xy, yz, zx;
        if (!(iss >> xx >> yy >> zz >> xy >> yz >> zx))
            continue;
        // Negate (FORCE on cell = -STRESS) and reorder to Voigt
        last_stress = {-xx, -yy, -zz, -yz, -zx, -xy};
    }
    return last_stress;
}

// Returns the final free energy e_fr_energy (eV) from vasprun.xml.
std::optional<double> readVasprunEnergy(const fs::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        return std::nullopt;

    std::optional<double> last_energy;
    std::string line;
    while (std::getline(f, line)) {
        if (line.find("e_fr_energy") == std::string::npos)
            continue;
        const auto start = line.find('>');
        const auto end = line.rfind('<');
        if (start == std::string::npos || end == std::string::npos || end <= start)
            continue;
        std::istringstream iss(line.substr(start + 1, end - start - 1));
        double e;
        if (iss >> e)
            last_energy = e;
    }
    return last_energy;
}

// Returns the final stress tensor (kBar) from vasprun.xml.
// vasprun.xml writes the actual Cauchy stress (tensile positive), no negation needed.
//
// Block layout:
//   <varray name="stress" >
//     <v>  σ_xx  σ_xy  σ_xz </v>
//     <v>  σ_yx  σ_yy  σ_yz </v>
//     <v>  σ_zx  σ_zy  σ_zz </v>
//   </varray>
std::optional<StressVoigt> readVasprunStress(const fs::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        return std::nullopt;

    std::optional<StressVoigt> last_stress;
    std::string line;
    while (std::getline(f, line)) {
        if (line.find(R"(name="stress")") == std::string::npos)
            continue;

        std::array<std::array<double, 3>, 3> mat{};
        bool ok = true;
        for (int row = 0; row < 3 && ok; ++row) {
            if (!std::getline(f, line)) {
                ok = false;
                break;
            }
            const auto s = line.find('>');
            const auto e = line.rfind('<');
            if (s == std::string::npos || e == std::string::npos || e <= s) {
                ok = false;
                break;
            }
            std::istringstream iss(line.substr(s + 1, e - s - 1));
            if (!(iss >> mat[row][0] >> mat[row][1] >> mat[row][2]))
                ok = false;
        }
        if (!ok)
            continue;

        // mat[i][j] = σ_{ij}, i,j ∈ {x=0, y=1, z=2}
        // Voigt: [σ_xx, σ_yy, σ_zz, σ_yz, σ_xz, σ_xy]
        last_stress = {mat[0][0], mat[1][1], mat[2][2], mat[1][2], mat[0][2], mat[0][1]};
    }
    return last_stress;
}

// ── Per-mode collected data ───────────────────────────────────────────────────

struct ModeData {
    ElasticStrainMode mode;
    std::vector<double> amplitudes;
    std::vector<double> energies;       // ΔE = E(δ) − E_ref  [eV]    (energy method)
    std::vector<StressVoigt> stresses;  // σ(δ) − σ_ref       [kBar]  (stress method)
};

// ── Fetch helpers ─────────────────────────────────────────────────────────────

static std::optional<double> fetchEnergy(const fs::path& dir, const std::string& source) {
    if (source == "oszicar")
        return readOszicar(dir / "OSZICAR");
    if (source == "outcar")
        return readOutcarEnergy(dir / "OUTCAR");
    return readVasprunEnergy(dir / "vasprun.xml");
}

static std::optional<StressVoigt> fetchStress(const fs::path& dir, const std::string& source) {
    if (source == "outcar")
        return readOutcarStress(dir / "OUTCAR");
    return readVasprunStress(dir / "vasprun.xml");
}

// ── Polynomial fit result ─────────────────────────────────────────────────────

struct PolyFitResult {
    double a2;             // quadratic coefficient (eV)
    double a4;             // quartic coefficient (eV)
    double rmse;           // root mean square error on residuals (eV)
    double r2;             // coefficient of determination
    double anharmonicity;  // |a4| * delta_max^2 / |a2|; >0.1 → significant anharmonicity
};

// ── Born stability result ─────────────────────────────────────────────────────

struct StabilityResult {
    bool stable;
    std::vector<double> eigenvalues;  // ascending eigenvalues of C (same units as C)
};

// ── VRH result ────────────────────────────────────────────────────────────────

struct VRHResult {
    double K_V, G_V;   // Voigt bulk and shear moduli
    double K_R, G_R;   // Reuss bulk and shear moduli
    double K_H, G_H;   // Hill (arithmetic mean) bulk and shear moduli
    double E_H, nu_H;  // Hill Young's modulus and Poisson ratio
    double A_U;        // Universal anisotropy: 5*(G_V/G_R) + (K_V/K_R) - 6
    double pugh;       // Pugh ratio K_H/G_H (>1.75 → ductile tendency)
};

// ── Linear algebra (LAPACKE) ──────────────────────────────────────────────────

// Invert a 6×6 matrix by solving M·X = I via LAPACKE_dgesv (no dgetri needed).
static std::optional<Cmat> invert6(Cmat M) {
    std::array<double, 36> flat_M, flat_I;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            flat_M[i * 6 + j] = M[i][j];
            flat_I[i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }

    std::array<lapack_int, 6> ipiv{};
    if (LAPACKE_dgesv(LAPACK_ROW_MAJOR, 6, 6, flat_M.data(), 6, ipiv.data(), flat_I.data(), 6) != 0)
        return std::nullopt;

    Cmat inv{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            inv[i][j] = flat_I[i * 6 + j];
    return inv;
}

// Condition number of an m×n matrix via SVD (LAPACKE_dgesvd). Returns -1 on failure.
static double conditionNumber(std::vector<double> A_flat, int m, int n) {
    const int k = std::min(m, n);
    std::vector<double> s(k);
    std::vector<double> superb(std::max(k - 1, 1));
    if (LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'N', m, n, A_flat.data(), n, s.data(), nullptr, 1, nullptr, 1,
                       superb.data()) != 0)
        return -1.0;
    if (s[k - 1] < 1e-15 * s[0])
        return std::numeric_limits<double>::infinity();
    return s[0] / s[k - 1];
}

// ── Stress-strain fitting with weighted LS ────────────────────────────────────

// Fit C_ij (kBar) from stress-strain data via weighted least squares (LAPACKE_dgels).
// Weight: w_k = delta_max / |delta_k| (smaller strain → higher weight for accuracy).
// Symmetrizes C_ij = C_ji afterwards.
static std::optional<Cmat> fitStressStrain(const std::vector<ModeData>& mode_data) {
    int n_data = 0;
    double delta_max = 0.0;
    for (const auto& md : mode_data) {
        n_data += static_cast<int>(md.amplitudes.size());
        for (double a : md.amplitudes)
            delta_max = std::max(delta_max, std::abs(a));
    }
    if (n_data < 6 || delta_max < 1e-15)
        return std::nullopt;

    std::vector<double> A_flat(n_data * 6, 0.0);
    std::vector<double> B_flat(n_data * 6, 0.0);
    int row = 0;
    for (const auto& md : mode_data) {
        for (size_t a = 0; a < md.amplitudes.size(); ++a) {
            const double amp = md.amplitudes[a];
            const double sqrtw = std::sqrt(delta_max / std::max(std::abs(amp), 1e-15));
            for (int j = 0; j < 6; ++j)
                A_flat[row * 6 + j] = sqrtw * md.mode.voigt[j] * amp;
            for (int i = 0; i < 6; ++i)
                B_flat[row * 6 + i] = sqrtw * md.stresses[a][i];
            ++row;
        }
    }

    const double cond = conditionNumber(A_flat, n_data, 6);
    if (cond < 0.0 || cond > 1e10)
        std::cerr << "Warning: ill-conditioned stress-strain system (cond = " << std::scientific << std::setprecision(2)
                  << cond << ")\n";

    if (LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', n_data, 6, 6, A_flat.data(), 6, B_flat.data(), 6) != 0)
        return std::nullopt;

    // First 6 rows of B_flat hold the solution: B_flat[j*6+i] = X[j][i] = C[i][j]
    Cmat C{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            C[i][j] = B_flat[j * 6 + i];

    // Symmetrize: C_ij = C_ji
    for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j)
            C[i][j] = C[j][i] = 0.5 * (C[i][j] + C[j][i]);

    return C;
}

// ── Energy-strain polynomial fitting ─────────────────────────────────────────

// Fit polynomial ΔE(δ) via weighted least squares and return diagnostics.
// symmetric=true  → fit a2*δ² + a4*δ⁴
// symmetric=false → fit a1*δ + a2*δ² + a3*δ³ + a4*δ⁴
// Weight: w_k = delta_max / |delta_k|
static std::optional<PolyFitResult> fitPolynomial(const std::vector<double>& amps, const std::vector<double>& dE,
                                                  bool symmetric) {
    const int m = static_cast<int>(amps.size());
    const int ncols = symmetric ? 2 : 4;
    const int a2col = symmetric ? 0 : 1;
    const int a4col = symmetric ? 1 : 3;

    if (m < ncols)
        return std::nullopt;

    double delta_max = 0.0;
    for (double a : amps)
        delta_max = std::max(delta_max, std::abs(a));

    std::vector<double> A_flat(m * ncols);
    std::vector<double> b(m);
    for (int k = 0; k < m; ++k) {
        const double d = amps[k];
        const double sqrtw = std::sqrt(delta_max / std::max(std::abs(d), 1e-15));
        b[k] = sqrtw * dE[k];
        if (symmetric) {
            A_flat[k * ncols + 0] = sqrtw * d * d;
            A_flat[k * ncols + 1] = sqrtw * d * d * d * d;
        } else {
            A_flat[k * ncols + 0] = sqrtw * d;
            A_flat[k * ncols + 1] = sqrtw * d * d;
            A_flat[k * ncols + 2] = sqrtw * d * d * d;
            A_flat[k * ncols + 3] = sqrtw * d * d * d * d;
        }
    }

    if (LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, ncols, 1, A_flat.data(), ncols, b.data(), 1) != 0)
        return std::nullopt;

    // Coefficients are in b[0..ncols-1] after solve
    const double a2 = b[a2col];
    const double a4 = b[a4col];

    // Compute unweighted RMSE and R² by evaluating polynomial at each point
    double mean_dE = 0.0;
    for (int k = 0; k < m; ++k)
        mean_dE += dE[k];
    mean_dE /= m;

    double ss_res = 0.0, ss_tot = 0.0;
    for (int k = 0; k < m; ++k) {
        const double d = amps[k];
        double fitted;
        if (symmetric) {
            fitted = b[0] * d * d + b[1] * d * d * d * d;
        } else {
            fitted = b[0] * d + b[1] * d * d + b[2] * d * d * d + b[3] * d * d * d * d;
        }
        const double res = dE[k] - fitted;
        ss_res += res * res;
        ss_tot += (dE[k] - mean_dE) * (dE[k] - mean_dE);
    }

    const double rmse = std::sqrt(ss_res / m);
    const double r2 = (ss_tot > 1e-30) ? (1.0 - ss_res / ss_tot) : 1.0;
    const double anharmonicity = (std::abs(a2) > 1e-30) ? std::abs(a4) * delta_max * delta_max / std::abs(a2) : 0.0;

    return PolyFitResult{a2, a4, rmse, r2, anharmonicity};
}

// ── Crystal system helpers ────────────────────────────────────────────────────

// Independent C_ij pairs (i≤j, 0-based) for each crystal system.
static std::vector<std::pair<int, int>> independentPairs(CrystalSystem cs) {
    using P = std::pair<int, int>;
    switch (cs) {
        case CrystalSystem::Cubic:
            return {P{0, 0}, P{0, 1}, P{3, 3}};
        case CrystalSystem::Hexagonal:
            return {P{0, 0}, P{0, 1}, P{0, 2}, P{2, 2}, P{3, 3}};
        case CrystalSystem::TetragonalII:
            return {P{0, 0}, P{0, 1}, P{0, 2}, P{2, 2}, P{3, 3}, P{5, 5}};
        case CrystalSystem::TetragonalI:
            return {P{0, 0}, P{0, 1}, P{0, 2}, P{2, 2}, P{3, 3}, P{5, 5}, P{0, 5}};
        case CrystalSystem::TrigonalII:
            return {P{0, 0}, P{0, 1}, P{0, 2}, P{2, 2}, P{3, 3}, P{0, 3}};
        case CrystalSystem::TrigonalI:
            return {P{0, 0}, P{0, 1}, P{0, 2}, P{2, 2}, P{3, 3}, P{0, 3}, P{0, 4}};
        case CrystalSystem::Orthorhombic:
            return {P{0, 0}, P{1, 1}, P{2, 2}, P{0, 1}, P{0, 2}, P{1, 2}, P{3, 3}, P{4, 4}, P{5, 5}};
        case CrystalSystem::Monoclinic:
            return {P{0, 0}, P{1, 1}, P{2, 2}, P{0, 1}, P{0, 2}, P{1, 2}, P{3, 3},
                    P{4, 4}, P{5, 5}, P{0, 4}, P{1, 4}, P{2, 4}, P{3, 5}};
        case CrystalSystem::Triclinic:
        default: {
            std::vector<P> v;
            for (int i = 0; i < 6; ++i)
                for (int j = i; j < 6; ++j)
                    v.push_back(P{i, j});
            return v;
        }
    }
}

// Fill symmetry-equivalent C elements from the independent ones.
static void applySymmetry(Cmat& C, CrystalSystem cs) {
    switch (cs) {
        case CrystalSystem::Cubic:
            C[1][1] = C[2][2] = C[0][0];
            C[0][2] = C[1][2] = C[0][1];
            C[4][4] = C[5][5] = C[3][3];
            break;
        case CrystalSystem::Hexagonal:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[5][5] = 0.5 * (C[0][0] - C[0][1]);
            break;
        case CrystalSystem::TetragonalII:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            break;
        case CrystalSystem::TetragonalI:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[1][5] = -C[0][5];
            break;
        case CrystalSystem::TrigonalII:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[5][5] = 0.5 * (C[0][0] - C[0][1]);
            C[1][3] = -C[0][3];
            C[4][5] = C[0][3];
            break;
        case CrystalSystem::TrigonalI:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[5][5] = 0.5 * (C[0][0] - C[0][1]);
            C[1][3] = -C[0][3];
            C[4][5] = C[0][3];
            C[1][4] = -C[0][4];
            C[3][5] = C[0][4];
            break;
        default:
            break;
    }
    // Enforce C_ij = C_ji
    for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j)
            C[j][i] = C[i][j];
}

// ── Energy-strain fitting ─────────────────────────────────────────────────────

// Fit C_ij (kBar) from energy-strain data.
// Prints per-mode polynomial fit diagnostics.
static std::optional<Cmat> fitEnergyStrain(const std::vector<ModeData>& mode_data, double volume_ang3,
                                           CrystalSystem cs) {
    const auto pairs = independentPairs(cs);
    const int n = static_cast<int>(pairs.size());

    if (static_cast<int>(mode_data.size()) != n)
        return std::nullopt;

    std::vector<double> a2(n);
    for (int m = 0; m < n; ++m) {
        auto res = fitPolynomial(mode_data[m].amplitudes, mode_data[m].energies, mode_data[m].mode.symmetric);
        if (!res)
            return std::nullopt;
        a2[m] = res->a2;

        std::cout << "    Mode " << std::setw(2) << (m + 1) << " [" << std::setw(12) << std::left
                  << mode_data[m].mode.label << std::right << "]  R²=" << std::fixed << std::setprecision(6) << res->r2
                  << "  RMSE=" << std::scientific << std::setprecision(2) << res->rmse << " eV"
                  << "  anharm=" << std::fixed << std::setprecision(4) << res->anharmonicity;
        if (res->anharmonicity > 0.1)
            std::cout << "  ** HIGH";
        std::cout << "\n";
    }

    // Build n×n coefficient matrix and rhs
    std::vector<double> A_flat(n * n), rhs(n);
    for (int m = 0; m < n; ++m) {
        const auto& v = mode_data[m].mode.voigt;
        for (int k = 0; k < n; ++k) {
            const auto [i, j] = pairs[k];
            A_flat[m * n + k] = (i == j ? 1.0 : 2.0) * v[i] * v[j];
        }
        rhs[m] = 2.0 * a2[m] / volume_ang3;
    }

    const double cond = conditionNumber(std::vector<double>(A_flat), n, n);
    if (cond < 0.0 || cond > 1e10)
        std::cerr << "Warning: ill-conditioned energy-strain system (cond = " << std::scientific << std::setprecision(2)
                  << cond << ")\n";

    std::vector<lapack_int> ipiv(n);
    if (LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A_flat.data(), n, ipiv.data(), rhs.data(), 1) != 0)
        return std::nullopt;

    Cmat C{};
    for (int k = 0; k < n; ++k) {
        const auto [i, j] = pairs[k];
        C[i][j] = rhs[k];
    }
    applySymmetry(C, cs);

    constexpr double eVA3_to_kBar = 1602.1766208;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            C[i][j] *= eVA3_to_kBar;

    return C;
}

// ── VRH polycrystalline averages ──────────────────────────────────────────────

// All inputs and outputs in the same units as C (GPa when called with C_GPa).
static std::optional<VRHResult> computeVRH(const Cmat& C) {
    const double K_V = (1.0 / 9.0) * (C[0][0] + C[1][1] + C[2][2] + 2.0 * (C[0][1] + C[0][2] + C[1][2]));
    const double G_V = (1.0 / 15.0) * (C[0][0] + C[1][1] + C[2][2] - (C[0][1] + C[0][2] + C[1][2]) +
                                       3.0 * (C[3][3] + C[4][4] + C[5][5]));

    auto S_opt = invert6(C);
    if (!S_opt)
        return std::nullopt;
    const Cmat& S = *S_opt;

    const double K_R = 1.0 / (S[0][0] + S[1][1] + S[2][2] + 2.0 * (S[0][1] + S[0][2] + S[1][2]));
    const double G_R = 15.0 / (4.0 * (S[0][0] + S[1][1] + S[2][2]) - 4.0 * (S[0][1] + S[0][2] + S[1][2]) +
                               3.0 * (S[3][3] + S[4][4] + S[5][5]));

    const double K_H = 0.5 * (K_V + K_R);
    const double G_H = 0.5 * (G_V + G_R);
    const double E_H = 9.0 * K_H * G_H / (3.0 * K_H + G_H);
    const double nu_H = (3.0 * K_H - 2.0 * G_H) / (2.0 * (3.0 * K_H + G_H));

    // Universal anisotropy index (Ranganathan & Ostoja-Starzewski, PRL 2008)
    const double A_U =
        (G_R > 1e-10 && K_R > 1e-10) ? 5.0 * (G_V / G_R) + (K_V / K_R) - 6.0 : std::numeric_limits<double>::quiet_NaN();
    const double pugh = (G_H > 1e-10) ? K_H / G_H : std::numeric_limits<double>::quiet_NaN();

    return VRHResult{K_V, G_V, K_R, G_R, K_H, G_H, E_H, nu_H, A_U, pugh};
}

// ── Born mechanical stability ─────────────────────────────────────────────────

// Check Born stability: all eigenvalues of the 6×6 C matrix must be positive.
// Uses LAPACKE_dsyev (symmetric eigenvalue decomposition).
static StabilityResult checkBornStability(const Cmat& C_GPa) {
    std::array<double, 36> flat;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            flat[i * 6 + j] = C_GPa[i][j];

    std::array<double, 6> w{};
    const int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', 6, flat.data(), 6, w.data());

    StabilityResult res;
    res.eigenvalues.assign(w.begin(), w.end());
    if (info != 0) {
        res.stable = false;
        return res;
    }
    res.stable = (w[0] > 0.0);
    return res;
}

// ── Output ────────────────────────────────────────────────────────────────────

static void writeCij(std::ostream& out, const Cmat& C, const ElasticDeformLog& log, const StabilityResult& stab,
                     const std::optional<VRHResult>& vrh) {
    out << "# Elastic constants C_ij (GPa)\n";
    out << "# Method: " << log.method << "-strain\n";
    out << "# Crystal: " << crystalSystemName(log.crystal_system) << " (" << log.space_group_symbol << " #"
        << log.space_group_number << ")\n";
    out << "# Volume: " << std::fixed << std::setprecision(4) << log.volume << " Ang^3\n";
    out << "#\n";
    out << "#        C_1j      C_2j      C_3j      C_4j      C_5j      C_6j\n";

    out << std::fixed << std::setprecision(3);
    for (int i = 0; i < 6; ++i) {
        out << "C_i" << (i + 1);
        for (int j = 0; j < 6; ++j)
            out << std::setw(10) << C[i][j];
        out << "\n";
    }

    out << "#\n# Born mechanical stability: " << (stab.stable ? "STABLE" : "UNSTABLE") << "\n";
    out << "#  Eigenvalues (GPa):";
    out << std::fixed << std::setprecision(3);
    for (double ev : stab.eigenvalues)
        out << "  " << std::setw(9) << ev;
    out << "\n";

    if (vrh) {
        out << "#\n# Voigt-Reuss-Hill polycrystalline averages (GPa)\n";
        out << std::fixed << std::setprecision(2);
        out << "#  K_Voigt = " << std::setw(8) << vrh->K_V << "   G_Voigt = " << std::setw(8) << vrh->G_V << "\n";
        out << "#  K_Reuss = " << std::setw(8) << vrh->K_R << "   G_Reuss = " << std::setw(8) << vrh->G_R << "\n";
        out << "#  K_Hill  = " << std::setw(8) << vrh->K_H << "   G_Hill  = " << std::setw(8) << vrh->G_H << "\n";
        out << "#  E_Hill  = " << std::setw(8) << vrh->E_H << "   nu_Hill = " << std::setw(8) << vrh->nu_H << "\n";
        if (!std::isnan(vrh->A_U))
            out << "#  Univ. anisotropy A_U   = " << std::fixed << std::setprecision(4) << vrh->A_U
                << "  (0 = isotropic)\n";
        if (!std::isnan(vrh->pugh))
            out << "#  Pugh ratio K_H/G_H     = " << std::fixed << std::setprecision(3) << vrh->pugh
                << "  (>1.75 → ductile tendency)\n";
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    CLI::App app{"Extract elastic constants from VASP calculations (vasp_elastic_fit)"};

    std::string logFile{"elastic_deform.log"};
    std::string dataSource;
    std::string outputFile{"elastic_constants.dat"};
    bool terminal_only{false};
    bool averages{false};

    app.add_option("--log,-l", logFile, "Path to elastic_deform.log manifest")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--source,-s", dataSource,
                   "Data source: oszicar (default for energy), outcar (default for stress, "
                   "or energy via sigma->0), vasprun")
        ->check(CLI::IsMember({"oszicar", "outcar", "vasprun"}));
    app.add_option("--output,-o", outputFile, "Output file for C_ij table")->capture_default_str();
    app.add_flag("--terminal,-t", terminal_only, "Print C_ij to terminal only; do not write output file");
    app.add_flag("--averages,-a", averages, "Compute and print Voigt-Reuss-Hill polycrystalline averages");

    CLI11_PARSE(app, argc, argv);

    // ── Read manifest ────────────────────────────────────────────────────────
    auto log_opt = readElasticLog(logFile);
    if (!log_opt) {
        std::cerr << "Error: cannot read manifest: " << logFile << "\n";
        return 1;
    }
    const ElasticDeformLog& log = *log_opt;

    fs::path log_dir = fs::path(logFile).parent_path();
    if (log_dir.empty())
        log_dir = ".";

    if (dataSource.empty())
        dataSource = (log.method == "energy") ? "oszicar" : "outcar";

    std::cout << "=== Elastic Constants Extraction ===\n";
    std::cout << "Method:      " << log.method << "-strain\n";
    std::cout << "Crystal:     " << crystalSystemName(log.crystal_system) << " (#" << log.space_group_number << " "
              << log.space_group_symbol << ")\n";
    std::cout << "Independent: " << log.n_independent << " constants\n";
    std::cout << "Source:      " << dataSource << "\n";
    std::cout << "Volume:      " << std::fixed << std::setprecision(4) << log.volume << " Ang^3\n\n";

    // ── Read reference data ──────────────────────────────────────────────────
    const fs::path ref_dir = log_dir / log.reference_dir;
    std::optional<double> ref_energy;
    std::optional<StressVoigt> ref_stress;

    if (log.method == "energy") {
        ref_energy = fetchEnergy(ref_dir, dataSource);
        if (!ref_energy) {
            std::cerr << "Error: cannot read reference energy from " << ref_dir << "\n";
            return 1;
        }
        std::cout << "Reference E0 = " << std::fixed << std::setprecision(8) << *ref_energy << " eV\n\n";
    } else {
        ref_stress = fetchStress(ref_dir, dataSource);
        if (!ref_stress) {
            std::cerr << "Warning: cannot read reference stress from " << ref_dir
                      << "; assuming zero residual stress.\n\n";
            ref_stress = StressVoigt{};
        } else {
            const double P_kBar = -((*ref_stress)[0] + (*ref_stress)[1] + (*ref_stress)[2]) / 3.0;
            std::cout << "Reference stress (kBar):" << "  xx=" << std::fixed << std::setprecision(3) << (*ref_stress)[0]
                      << "  yy=" << (*ref_stress)[1] << "  zz=" << (*ref_stress)[2] << "\n";
            std::cout << "Residual pressure: " << std::fixed << std::setprecision(3) << P_kBar * 0.1 << " GPa";
            if (std::abs(P_kBar * 0.1) > 0.5)
                std::cout << "  ** Born-Huang correction will be applied";
            std::cout << "\n\n";
        }
    }

    // ── Collect data from deformed structures ────────────────────────────────
    std::vector<ModeData> mode_data;
    mode_data.reserve(log.modes.size());
    bool any_error = false;

    for (size_t m = 0; m < log.modes.size(); ++m) {
        ModeData md;
        md.mode = log.modes[m];

        for (size_t a = 0; a < log.amplitudes[m].size(); ++a) {
            const fs::path dir = log_dir / log.dirs[m][a];

            if (log.method == "energy") {
                auto e = fetchEnergy(dir, dataSource);
                if (!e) {
                    std::cerr << "Error: cannot read energy from " << dir << "\n";
                    any_error = true;
                } else {
                    md.amplitudes.push_back(log.amplitudes[m][a]);
                    md.energies.push_back(*e - *ref_energy);
                }
            } else {
                auto s = fetchStress(dir, dataSource);
                if (!s) {
                    std::cerr << "Error: cannot read stress from " << dir << "\n";
                    any_error = true;
                } else {
                    for (int i = 0; i < 6; ++i)
                        (*s)[i] -= (*ref_stress)[i];
                    md.amplitudes.push_back(log.amplitudes[m][a]);
                    md.stresses.push_back(*s);
                }
            }
        }

        const int expected = static_cast<int>(log.amplitudes[m].size());
        const int got = static_cast<int>(md.amplitudes.size());
        std::cout << "  Mode " << std::setw(2) << (m + 1) << " [" << std::setw(12) << std::left << log.modes[m].label
                  << std::right << "]:  " << got << "/" << expected << " structures";
        if (got < expected)
            std::cout << "  ** INCOMPLETE";
        std::cout << "\n";

        mode_data.push_back(std::move(md));
    }

    if (any_error)
        std::cerr << "\nWarning: some structures could not be read. Results may be incomplete.\n";

    // ── Mirror symmetric modes ───────────────────────────────────────────────
    // poscar_elastic_deform generates only positive amplitudes for symmetric modes.
    // Mirror: energy ΔE(−δ)=ΔE(+δ); stress σ(−δ)=−σ(+δ).
    int n_mirrored = 0;
    for (auto& md : mode_data) {
        if (!md.mode.symmetric)
            continue;
        bool all_positive = std::all_of(md.amplitudes.begin(), md.amplitudes.end(), [](double a) { return a > 0.0; });
        if (!all_positive)
            continue;

        const int n = static_cast<int>(md.amplitudes.size());
        for (int k = n - 1; k >= 0; --k) {
            md.amplitudes.push_back(-md.amplitudes[k]);
            if (log.method == "energy") {
                md.energies.push_back(md.energies[k]);
            } else {
                StressVoigt neg_s;
                for (int i = 0; i < 6; ++i)
                    neg_s[i] = -md.stresses[k][i];
                md.stresses.push_back(neg_s);
            }
            ++n_mirrored;
        }
    }
    if (n_mirrored > 0)
        std::cout << "\nMirrored " << n_mirrored << " data points for symmetric modes.\n";

    // ── Fit C_ij ─────────────────────────────────────────────────────────────
    std::cout << "\nFitting C_ij ...\n";
    std::optional<Cmat> C_kBar_opt;
    if (log.method == "stress") {
        C_kBar_opt = fitStressStrain(mode_data);
    } else {
        C_kBar_opt = fitEnergyStrain(mode_data, log.volume, log.crystal_system);
    }
    if (!C_kBar_opt) {
        std::cerr << "Error: fitting failed (singular system or insufficient data).\n";
        return 1;
    }

    // Convert kBar → GPa
    Cmat C_GPa{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            C_GPa[i][j] = (*C_kBar_opt)[i][j] * 0.1;

    // ── Pressure correction (stress method, Born-Huang) ──────────────────────
    // C_true[i][j] = C_meas[i][j] + P * delta[i][j]
    // delta: +1 on diagonal, -1 for off-diagonal normal pairs (i,j both < 3), 0 otherwise
    if (log.method == "stress" && ref_stress) {
        const double P_kBar = -((*ref_stress)[0] + (*ref_stress)[1] + (*ref_stress)[2]) / 3.0;
        const double P_GPa = P_kBar * 0.1;
        if (std::abs(P_GPa) > 0.5) {
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    int delta = 0;
                    if (i == j)
                        delta = 1;
                    else if (i < 3 && j < 3)
                        delta = -1;
                    C_GPa[i][j] += P_GPa * delta;
                }
            }
            std::cout << "Born-Huang pressure correction applied (" << std::fixed << std::setprecision(2) << P_GPa
                      << " GPa).\n";
        }
    }

    // ── Kleinman / Cauchy analysis (cubic only) ──────────────────────────────
    if (log.crystal_system == CrystalSystem::Cubic) {
        const double C11 = C_GPa[0][0];
        const double C12 = C_GPa[0][1];
        const double C44 = C_GPa[3][3];
        const double denom = 7.0 * C11 + 2.0 * C12;
        std::cout << "\n--- Cubic analysis ---\n";
        std::cout << "Cauchy pressure (C12-C44) = " << std::fixed << std::setprecision(2) << (C12 - C44) << " GPa"
                  << "  (" << ((C12 - C44) > 0 ? "ductile" : "brittle") << " tendency)\n";
        if (std::abs(denom) > 1e-10)
            std::cout << "Kleinman parameter zeta = " << std::fixed << std::setprecision(4) << (C11 + 8.0 * C12) / denom
                      << "  (0=no internal strain, 1=max)\n";
    }

    // ── Born stability check ─────────────────────────────────────────────────
    const StabilityResult stab = checkBornStability(C_GPa);
    std::cout << "\nBorn stability: " << (stab.stable ? "STABLE" : "** UNSTABLE **") << "\n";

    // ── Optional VRH averages ────────────────────────────────────────────────
    std::optional<VRHResult> vrh;
    if (averages) {
        vrh = computeVRH(C_GPa);
        if (!vrh)
            std::cerr << "Warning: C_ij matrix is singular; VRH averages cannot be computed.\n";
    }

    // ── Write output ─────────────────────────────────────────────────────────
    if (!terminal_only) {
        std::ofstream fout(outputFile);
        if (!fout) {
            std::cerr << "Error: cannot open output file: " << outputFile << "\n";
            return 1;
        }
        writeCij(fout, C_GPa, log, stab, vrh);
        std::cout << "Written: " << outputFile << "\n";
    }

    // Always echo the table to terminal
    writeCij(std::cout, C_GPa, log, stab, vrh);

    return 0;
}
