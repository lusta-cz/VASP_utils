#include <lapacke.h>

#include <array>
#include <CLI/CLI.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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
// gives -σ (compressive positive). We negate to recover the Cauchy stress (tensile positive)
// so that C_ij > 0 for a mechanically stable crystal.
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

// ── Linear algebra (LAPACKE) ──────────────────────────────────────────────────

// Invert a 6×6 matrix in-place using LU factorisation (LAPACKE_dgetrf/dgetri).
// Returns false if singular.
static std::optional<Cmat> invert6(Cmat M) {
    // Flatten row-major
    std::array<double, 36> flat;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            flat[i * 6 + j] = M[i][j];

    std::array<lapack_int, 6> ipiv{};
    if (LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 6, 6, flat.data(), 6, ipiv.data()) != 0)
        return std::nullopt;
    if (LAPACKE_dgetri(LAPACK_ROW_MAJOR, 6, flat.data(), 6, ipiv.data()) != 0)
        return std::nullopt;

    Cmat inv{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            inv[i][j] = flat[i * 6 + j];
    return inv;
}

// ── Stress-strain fitting ─────────────────────────────────────────────────────

// Fit C_ij (kBar) from stress-strain data via least squares (LAPACKE_dgels).
// For each data point k: strain row ε[j] = voigt[j]*amp, stress column σ[i].
// Solves min ||A*X - B|| where A is (n_data×6), B is (n_data×6).
// Solution X[j][i] = C[i][j]; symmetrized afterwards.
static std::optional<Cmat> fitStressStrain(const std::vector<ModeData>& mode_data) {
    // Build flat row-major A (n_data×6) and B (n_data×6)
    std::vector<double> A_flat, B_flat;
    int n_data = 0;
    for (const auto& md : mode_data) {
        for (size_t a = 0; a < md.amplitudes.size(); ++a) {
            for (int j = 0; j < 6; ++j)
                A_flat.push_back(md.mode.voigt[j] * md.amplitudes[a]);
            for (int i = 0; i < 6; ++i)
                B_flat.push_back(md.stresses[a][i]);
            ++n_data;
        }
    }
    if (n_data < 6)
        return std::nullopt;

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

// ── Energy-strain fitting ─────────────────────────────────────────────────────

// Fit polynomial ΔE(δ) and return the quadratic coefficient a₂ (eV).
// symmetric=true → fit a₂δ² + a₄δ⁴  (only even powers, matches the mode symmetry)
// symmetric=false → fit a₁δ + a₂δ² + a₃δ³ + a₄δ⁴
// Uses LAPACKE_dgels for the overdetermined least-squares problem.
static std::optional<double> fitPolynomial(const std::vector<double>& amps, const std::vector<double>& dE,
                                           bool symmetric) {
    const int m = static_cast<int>(amps.size());
    const int ncols = symmetric ? 2 : 4;  // number of basis functions
    const int a2col = symmetric ? 0 : 1;  // column index of the δ² term

    if (m < ncols)
        return std::nullopt;

    // Build design matrix A (m×ncols) row-major and rhs b (m×1)
    std::vector<double> A_flat(m * ncols), b(m);
    for (int k = 0; k < m; ++k) {
        const double d = amps[k];
        b[k] = dE[k];
        if (symmetric) {
            A_flat[k * ncols + 0] = d * d;          // δ²
            A_flat[k * ncols + 1] = d * d * d * d;  // δ⁴
        } else {
            A_flat[k * ncols + 0] = d;              // δ
            A_flat[k * ncols + 1] = d * d;          // δ²
            A_flat[k * ncols + 2] = d * d * d;      // δ³
            A_flat[k * ncols + 3] = d * d * d * d;  // δ⁴
        }
    }

    if (LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, ncols, 1, A_flat.data(), ncols, b.data(), 1) != 0)
        return std::nullopt;

    return b[a2col];  // solution is in first ncols rows of b
}

// Independent C_ij pairs (i≤j, 0-based) for each crystal system.
// These are exactly the free elements solved for by energyStrainModes().
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
            // unique axis b (standard setting)
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
            C[1][5] = -C[0][5];  // C26 = -C16
            break;
        case CrystalSystem::TrigonalII:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[5][5] = 0.5 * (C[0][0] - C[0][1]);
            C[1][3] = -C[0][3];  // C24 = -C14
            C[4][5] = C[0][3];   // C56 = C14
            break;
        case CrystalSystem::TrigonalI:
            C[1][1] = C[0][0];
            C[1][2] = C[0][2];
            C[4][4] = C[3][3];
            C[5][5] = 0.5 * (C[0][0] - C[0][1]);
            C[1][3] = -C[0][3];  // C24 = -C14
            C[4][5] = C[0][3];   // C56 = C14
            C[1][4] = -C[0][4];  // C25 = -C15
            C[3][5] = C[0][4];   // C46 = C15
            break;
        default:
            break;  // orthorhombic, monoclinic, triclinic: no extra constraints
    }
    // Enforce C_ij = C_ji
    for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j)
            C[j][i] = C[i][j];
}

// Fit C_ij (kBar) from energy-strain data.
//
// For each mode m: a₂_m = (V₀/2) * Σ_{ij} C_{ij} * v_m[i] * v_m[j]
// → build square n×n system (n = n_independent) and solve with LAPACKE_dgesv.
// Result is converted from eV/Å³ to kBar (× 1602.1766).
static std::optional<Cmat> fitEnergyStrain(const std::vector<ModeData>& mode_data, double volume_ang3,
                                           CrystalSystem cs) {
    const auto pairs = independentPairs(cs);
    const int n = static_cast<int>(pairs.size());

    if (static_cast<int>(mode_data.size()) != n)
        return std::nullopt;

    // Polynomial fit → a₂ for each mode
    std::vector<double> a2(n);
    for (int m = 0; m < n; ++m) {
        auto val = fitPolynomial(mode_data[m].amplitudes, mode_data[m].energies, mode_data[m].mode.symmetric);
        if (!val)
            return std::nullopt;
        a2[m] = *val;
    }

    // Coefficient matrix A[m][k] and rhs[m] = 2*a2[m]/V0  (eV/Å³)
    // A[m][k] = factor * v_m[i_k] * v_m[j_k],  factor=1 diagonal, 2 off-diagonal
    std::vector<double> A_flat(n * n), rhs(n);
    for (int m = 0; m < n; ++m) {
        const auto& v = mode_data[m].mode.voigt;
        for (int k = 0; k < n; ++k) {
            const auto [i, j] = pairs[k];
            A_flat[m * n + k] = (i == j ? 1.0 : 2.0) * v[i] * v[j];
        }
        rhs[m] = 2.0 * a2[m] / volume_ang3;
    }

    std::vector<lapack_int> ipiv(n);
    if (LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A_flat.data(), n, ipiv.data(), rhs.data(), 1) != 0)
        return std::nullopt;

    // Fill C from solution and apply crystal symmetry
    Cmat C{};
    for (int k = 0; k < n; ++k) {
        const auto [i, j] = pairs[k];
        C[i][j] = rhs[k];
    }
    applySymmetry(C, cs);

    // Convert eV/Å³ → kBar  (1 eV/Å³ = 1602.1766208 kBar)
    constexpr double eVA3_to_kBar = 1602.1766208;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            C[i][j] *= eVA3_to_kBar;

    return C;
}

// ── Voigt-Reuss-Hill polycrystalline averages ─────────────────────────────────

struct VRHResult {
    double K_V, G_V;   // Voigt bulk and shear moduli
    double K_R, G_R;   // Reuss bulk and shear moduli
    double K_H, G_H;   // Hill (arithmetic mean) bulk and shear moduli
    double E_H, nu_H;  // Hill Young's modulus and Poisson ratio
};

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

    return VRHResult{K_V, G_V, K_R, G_R, K_H, G_H, E_H, nu_H};
}

// ── Output ────────────────────────────────────────────────────────────────────

static void writeCij(std::ostream& out, const Cmat& C, const ElasticDeformLog& log,
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

    if (vrh) {
        out << "#\n# Voigt-Reuss-Hill polycrystalline averages (GPa)\n";
        out << std::fixed << std::setprecision(2);
        out << "#  K_Voigt = " << std::setw(8) << vrh->K_V << "   G_Voigt = " << std::setw(8) << vrh->G_V << "\n";
        out << "#  K_Reuss = " << std::setw(8) << vrh->K_R << "   G_Reuss = " << std::setw(8) << vrh->G_R << "\n";
        out << "#  K_Hill  = " << std::setw(8) << vrh->K_H << "   G_Hill  = " << std::setw(8) << vrh->G_H << "\n";
        out << "#  E_Hill  = " << std::setw(8) << vrh->E_H << "   nu_Hill = " << std::setw(8) << vrh->nu_H << "\n";
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
            std::cout << "Reference stress (kBar):" << "  xx=" << std::fixed << std::setprecision(3) << (*ref_stress)[0]
                      << "  yy=" << (*ref_stress)[1] << "  zz=" << (*ref_stress)[2] << "\n\n";
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
        writeCij(fout, C_GPa, log, vrh);
        std::cout << "Written: " << outputFile << "\n";
    }

    // Always echo the table to terminal
    writeCij(std::cout, C_GPa, log, vrh);

    return 0;
}
