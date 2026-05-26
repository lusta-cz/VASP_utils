#include <algorithm>
#include <array>
#include <CLI/CLI.hpp>
#define _USE_MATH_DEFINES
#include <lapacke.h>

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

#include "atomic_weights.h"
#include "elastic.h"
#include "error_parse.h"

// 6×6 elastic constant matrix, C[i][j] in Voigt notation (0-based).
using Cmat = std::array<std::array<double, 6>, 6>;

// ── VRH result ────────────────────────────────────────────────────────────────

struct VRHResult {
    double K_V, G_V;   // Voigt bulk and shear moduli
    double K_R, G_R;   // Reuss bulk and shear moduli
    double K_H, G_H;   // Hill (arithmetic mean) bulk and shear moduli
    double E_H, nu_H;  // Hill Young's modulus and Poisson ratio
    double A_U;        // Universal anisotropy: 5*(G_V/G_R) + (K_V/K_R) - 6
    double pugh;       // Pugh ratio K_H/G_H (>1.75 → ductile tendency)
};

struct SoundVelocitiesResult {
    // Isotropic method (Anderson 1963)
    double v_l_iso{0.0};
    double v_t_iso{0.0};
    double v_m_iso{0.0};

    // Anisotropic method (Christoffel integration)
    double v_l_aniso{0.0};
    double v_t_aniso{0.0};
    double v_m_aniso{0.0};
};

/**
 * @brief Helper to map 4th-rank tensor indices (0,1,2 for x,y,z) to 6x6 Voigt indices.
 */
inline int tensorToVoigt(int i, int j) {
    if (i == j)
        return i;  // 00->0 (C11), 11->1 (C22), 22->2 (C33)
    if ((i == 1 && j == 2) || (i == 2 && j == 1))
        return 3;  // 23, 32 -> 3 (C44)
    if ((i == 0 && j == 2) || (i == 2 && j == 0))
        return 4;  // 03, 30 -> 4 (C55)
    if ((i == 0 && j == 1) || (i == 1 && j == 0))
        return 5;  // 01, 10 -> 5 (C66)
    return 0;
}

// ── Thermal properties result structure ──────────────────────────────────────

struct ThermalPropertiesResult {
    double debye_iso{0.0};           // Debye temperature from Isotropic velocity (K)
    double debye_aniso{0.0};         // Debye temperature from Anisotropic velocity (K)
    double gruneisen{0.0};           // Dimensionless acoustic Grüneisen parameter
    double kappa_clarke{0.0};        // Clarke min thermal conductivity (W/m*K)
    double kappa_cahill_iso{0.0};    // Cahill min thermal conductivity (Isotropic velocities)
    double kappa_cahill_aniso{0.0};  // Cahill min thermal conductivity (Anisotropic velocities)
    double alpha_V{0.0};             // Volumetric thermal expansion coefficient (1/K, High-T limit)
};

/**
 * @brief Computes sound velocities using both Isotropic and Anisotropic approaches.
 * @param C Elastic stiffness matrix C_ij in GPa.
 * @param density Mass density in g/cm^3.
 * @param vrh Polycrystalline Hill averages used for the isotropic model.
 */
SoundVelocitiesResult compute_sound_velocities(const Cmat& C, double density, const VRHResult& vrh) {
    SoundVelocitiesResult res_vel;
    if (density <= 0.0)
        return res_vel;

    // ─────────────────────────────────────────────────────────────────────────
    // APPROACH 1: Isotropic baseline (Anderson 1963)
    // ─────────────────────────────────────────────────────────────────────────
    // Conversion factor from (GPa / (g/cm^3)) to m/s is exactly 1000.0
    double factor = 1000.0;
    res_vel.v_l_iso = std::sqrt((vrh.K_H + (4.0 / 3.0) * vrh.G_H) / density) * factor;
    res_vel.v_t_iso = std::sqrt(vrh.G_H / density) * factor;

    double inv_vm3_iso =
        (1.0 / 3.0) * ((1.0 / std::pow(res_vel.v_l_iso, 3.0)) + (2.0 / std::pow(res_vel.v_t_iso, 3.0)));
    res_vel.v_m_iso = std::pow(inv_vm3_iso, -1.0 / 3.0);

    // ─────────────────────────────────────────────────────────────────────────
    // APPROACH 2: Anisotropic Integration (Christoffel Equation)
    // ─────────────────────────────────────────────────────────────────────────
    double sum_inv_v3 = 0.0;
    double sum_vl = 0.0;
    double sum_vt = 0.0;
    double total_weight = 0.0;

    // High density grid for numerical integration over the unit sphere
    const int n_theta = 45;
    const int n_phi = 90;
    const double d_theta = M_PI / n_theta;
    const double d_phi = 2.0 * M_PI / n_phi;

    for (int i = 0; i < n_theta; ++i) {
        double theta = d_theta * (i + 0.5);
        double sin_theta = std::sin(theta);
        double weight = sin_theta * d_theta * d_phi;  // Area element dA

        for (int j = 0; j < n_phi; ++j) {
            double phi = d_phi * (j + 0.5);

            // Unit propagation vector n
            std::array<double, 3> n = {sin_theta * std::cos(phi), sin_theta * std::sin(phi), std::cos(theta)};

            // Build Christoffel Matrix M_ik = C_ijkl * n_j * n_l (Symmetric 3x3)
            double M[9] = {0.0};
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    double val = 0.0;
                    for (int j_idx = 0; j_idx < 3; ++j_idx) {
                        for (int l_idx = 0; l_idx < 3; ++l_idx) {
                            int alpha = tensorToVoigt(r, j_idx);
                            int beta = tensorToVoigt(c, l_idx);
                            val += C[alpha][beta] * n[j_idx] * n[l_idx];
                        }
                    }
                    M[r * 3 + c] = val;
                }
            }

            // Solve eigenvalues using LAPACKE (Symmetric Matrix Eigensolver)
            double eigenvalues[3] = {0.0};
            int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', 3, M, 3, eigenvalues);

            if (info == 0) {
                // Convert eigenvalues (GPa) to absolute velocities (m/s)
                double v0 = std::sqrt(std::max(0.0, eigenvalues[0]) / density) * 1000.0;  // Quasi-transverse 1
                double v1 = std::sqrt(std::max(0.0, eigenvalues[1]) / density) * 1000.0;  // Quasi-transverse 2
                double v2 = std::sqrt(std::max(0.0, eigenvalues[2]) / density) * 1000.0;  // Quasi-longitudinal

                if (v0 > 1e-3 && v1 > 1e-3 && v2 > 1e-3) {
                    sum_inv_v3 +=
                        weight * ((1.0 / std::pow(v0, 3.0)) + (1.0 / std::pow(v1, 3.0)) + (1.0 / std::pow(v2, 3.0)));
                    sum_vt += weight * (v0 + v1);
                    sum_vl += weight * v2;
                    total_weight += weight;
                }
            }
        }
    }

    if (total_weight > 0.0) {
        // Average the anisotropic branches globally over the integrated surface
        double mean_inv_v3 = sum_inv_v3 / (3.0 * total_weight);
        res_vel.v_m_aniso = std::pow(mean_inv_v3, -1.0 / 3.0);

        res_vel.v_l_aniso = sum_vl / total_weight;
        res_vel.v_t_aniso = sum_vt / (2.0 * total_weight);  // 2 transverse branches
    }

    return res_vel;
}

/**
 * @brief Computes thermal parameters using isotropic and anisotropic inputs.
 * @param vrh Polycrystalline moduli averages (K_H, G_H, E_H, nu_H in GPa).
 * @param vel Calculated sound velocities container.
 * @param total_atoms Total number of atoms inside the unit cell.
 * @param volume Cell volume in cubic Angstroms (Å³).
 * @param density Mass density in g/cm³.
 */
ThermalPropertiesResult compute_thermal_properties(const VRHResult& vrh, const SoundVelocitiesResult& vel,
                                                   int total_atoms, double volume, double density) {
    ThermalPropertiesResult res_thermal;
    if (volume <= 0.0 || total_atoms <= 0 || density <= 0.0)
        return res_thermal;

    // --- Fundamental Physical Constants ---
    constexpr double h_planck = 6.62607015e-34;   // J*s
    constexpr double k_boltzmann = 1.380649e-23;  // J/K

    // Number density: atoms per m³
    // (total_atoms / (volume * 1e-30 m³))
    double atom_density_m3 = static_cast<double>(total_atoms) / (volume * 1e-30);

    // ─────────────────────────────────────────────────────────────────────────
    // 1. DEBYE TEMPERATURES (θ_D)
    // ─────────────────────────────────────────────────────────────────────────
    double debye_factor = (h_planck / k_boltzmann) * std::pow((3.0 * atom_density_m3) / (4.0 * M_PI), 1.0 / 3.0);

    res_thermal.debye_iso = debye_factor * vel.v_m_iso;
    res_thermal.debye_aniso = debye_factor * vel.v_m_aniso;

    // ─────────────────────────────────────────────────────────────────────────
    // 2. GRÜNEISEN PARAMETER (γ)
    // ─────────────────────────────────────────────────────────────────────────
    // Handle potential denominator divergence if Poisson ratio approaches ~0.66
    double gruneisen_denom = 2.0 - 3.0 * vrh.nu_H;
    if (std::abs(gruneisen_denom) > 1e-5) {
        res_thermal.gruneisen = 1.5 * ((1.0 + vrh.nu_H) / gruneisen_denom);
    } else {
        res_thermal.gruneisen = 0.0;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // 3. MINIMUM THERMAL CONDUCTIVITY (κ_min)
    // ─────────────────────────────────────────────────────────────────────────
    double density_factor_23 = std::pow(atom_density_m3, 2.0 / 3.0);

    // A. Clarke Model
    // E_Pa = vrh.E_H * 1e9, density_kg_m3 = density * 1000.0
    double E_over_rho = (vrh.E_H * 1e9) / (density * 1000.0);
    res_thermal.kappa_clarke = 0.87 * k_boltzmann * density_factor_23 * std::sqrt(E_over_rho);

    // B. Cahill Model (Isotropic velocity profile branch)
    res_thermal.kappa_cahill_iso = 0.5 * k_boltzmann * density_factor_23 * (vel.v_l_iso + 2.0 * vel.v_t_iso);

    // C. Cahill Model (Anisotropic velocity profile branch)
    res_thermal.kappa_cahill_aniso = 0.5 * k_boltzmann * density_factor_23 * (vel.v_l_aniso + 2.0 * vel.v_t_aniso);

    // ─────────────────────────────────────────────────────────────────────────
    // 4. THERMAL EXPANSION COEFFICIENT (α_V) - High Temperature Limit
    // ─────────────────────────────────────────────────────────────────────────
    // Cv per unit volume in classical limit = 3 * atom_density_m3 * k_boltzmann
    double Cv_per_volume = 3.0 * atom_density_m3 * k_boltzmann;

    // Bulk modulus in Pa = vrh.K_H * 1e9
    double K_Pa = vrh.K_H * 1e9;

    if (K_Pa > 0.0) {
        res_thermal.alpha_V = (res_thermal.gruneisen * Cv_per_volume) / K_Pa;
    }

    return res_thermal;
}

//--Input
static bool fetchCij(std::ifstream& file, Cmat& C, double& volume, ParseError& error) {
    std::string line;
    int line_number = 0;
    bool volume_found = false;

    // 1. Skip lines until we find the start of the C_ij matrix, scanning for Volume along the way
    while (readRequiredLine(file, line, line_number, error, "matrix search")) {
        // Check if this line contains the volume property
        if (line.rfind("# Volume: ", 0) == 0) {
            std::stringstream ss_vol(line.substr(10));
            if (ss_vol >> volume) {
                volume_found = true;
            }
        }

        // Look for a line starting with "C_i1" (ignoring leading whitespace)
        std::stringstream ss(line);
        std::string first_token;
        if (ss >> first_token && first_token == "C_i1") {
            // Check if we found volume before the matrix started
            if (!volume_found) {
                return fail(error, ParseErrorKind::Semantic, line_number,
                            "could not find or parse '# Volume:' metadata prefix before matrix data");
            }

            // Found the first row! Process it immediately.
            for (int j = 0; j < 6; ++j) {
                if (!(ss >> C[0][j])) {
                    return fail(error, ParseErrorKind::Parse, line_number,
                                "failed to parse matrix element C_1" + std::to_string(j + 1));
                }
            }
            break;  // Break out of the skip loop to handle the remaining 5 rows
        }
    }

    // If we dropped out of the loop without finding "C_i1", readRequiredLine already set the error
    if (line_number == 0 || (C[0][0] == 0 && line.find("C_i1") == std::string::npos)) {
        return fail(error, ParseErrorKind::Parse, line_number, "Could not find start of C_ij matrix data");
    }

    // 2. Read the remaining 5 rows (rows index 1 to 5)
    for (int i = 1; i < 6; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "matrix row " + std::to_string(i + 1))) {
            return false;
        }

        std::stringstream ss(line);
        std::string row_label;

        if (!(ss >> row_label)) {
            return fail(error, ParseErrorKind::Parse, line_number, "missing row label");
        }

        // Validate that rows are sequential and not malformed
        if (row_label != "C_i" + std::to_string(i + 1)) {
            return fail(error, ParseErrorKind::Parse, line_number,
                        "expected row label 'C_i" + std::to_string(i + 1) + "', got '" + row_label + "'");
        }

        // Parse the 6 numeric values for this row
        for (int j = 0; j < 6; ++j) {
            if (!(ss >> C[i][j])) {
                return fail(error, ParseErrorKind::Parse, line_number,
                            "failed to parse matrix element C_" + std::to_string(i + 1) + std::to_string(j + 1));
            }
        }
    }

    return true;
}

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

// ── VRH polycrystalline averages ──────────────────────────────────────────────

// All inputs and outputs in the same units as C (GPa when called with C_GPa).
static bool computeVRH(const Cmat& C, VRHResult& res_vrh) {
    res_vrh.K_V = (1.0 / 9.0) * (C[0][0] + C[1][1] + C[2][2] + 2.0 * (C[0][1] + C[0][2] + C[1][2]));
    res_vrh.G_V = (1.0 / 15.0) *
                  (C[0][0] + C[1][1] + C[2][2] - (C[0][1] + C[0][2] + C[1][2]) + 3.0 * (C[3][3] + C[4][4] + C[5][5]));

    auto S_opt = invert6(C);
    if (!S_opt) {
        return false;  // Matrix inversion failed, return failure state
    }
    const Cmat& S = *S_opt;

    res_vrh.K_R = 1.0 / (S[0][0] + S[1][1] + S[2][2] + 2.0 * (S[0][1] + S[0][2] + S[1][2]));
    res_vrh.G_R = 15.0 / (4.0 * (S[0][0] + S[1][1] + S[2][2]) - 4.0 * (S[0][1] + S[0][2] + S[1][2]) +
                          3.0 * (S[3][3] + S[4][4] + S[5][5]));

    res_vrh.K_H = 0.5 * (res_vrh.K_V + res_vrh.K_R);
    res_vrh.G_H = 0.5 * (res_vrh.G_V + res_vrh.G_R);
    res_vrh.E_H = 9.0 * res_vrh.K_H * res_vrh.G_H / (3.0 * res_vrh.K_H + res_vrh.G_H);
    res_vrh.nu_H = (3.0 * res_vrh.K_H - 2.0 * res_vrh.G_H) / (2.0 * (3.0 * res_vrh.K_H + res_vrh.G_H));

    res_vrh.A_U = (res_vrh.G_R > 1e-10 && res_vrh.K_R > 1e-10)
                      ? 5.0 * (res_vrh.G_V / res_vrh.G_R) + (res_vrh.K_V / res_vrh.K_R) - 6.0
                      : std::numeric_limits<double>::quiet_NaN();
    res_vrh.pugh = (res_vrh.G_H > 1e-10) ? res_vrh.K_H / res_vrh.G_H : std::numeric_limits<double>::quiet_NaN();

    return true;
}

static void writeVRH(std::ostream& out, const std::optional<VRHResult>& vrh) {
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

static void writeVelocitiesOutput(std::ostream& out, const SoundVelocitiesResult& vel) {
    out << "#\n# Isotropic & Anisotropic Acoustic Sound Velocities (m/s)\n";
    out << std::fixed << std::setprecision(1);
    out << "#  Isotropic (Anderson 1963):\n";
    out << "#    v_longitudinal = " << std::setw(8) << vel.v_l_iso << " m/s\n";
    out << "#    v_transverse   = " << std::setw(8) << vel.v_t_iso << " m/s\n";
    out << "#    v_mean         = " << std::setw(8) << vel.v_m_iso << " m/s\n";
    out << "#  Anisotropic (Christoffel Integration):\n";
    out << "#    v_longitudinal = " << std::setw(8) << vel.v_l_aniso << " m/s\n";
    out << "#    v_transverse   = " << std::setw(8) << vel.v_t_aniso << " m/s\n";
    out << "#    v_mean         = " << std::setw(8) << vel.v_m_aniso << " m/s\n";
}

static void writeThermalOutput(std::ostream& out, const ThermalPropertiesResult& therm) {
    out << "#\n# Thermal & Quasi-Harmonic Thermodynamic Properties\n";
    out << std::fixed << std::setprecision(1);
    out << "#  Debye Temperature (Isotropic v_m)   = " << std::setw(6) << therm.debye_iso << " K\n";
    out << "#  Debye Temperature (Anisotropic v_m) = " << std::setw(6) << therm.debye_aniso << " K\n";
    out << std::setprecision(3);
    out << "#  Acoustic Gruneisen Parameter (γ)    = " << std::setw(6) << therm.gruneisen << "\n";
    out << std::setprecision(2);
    out << "#  Minimum Thermal Conductivity (W/m*K):\n";
    out << "#    Clarke Model                      = " << std::setw(6) << therm.kappa_clarke << " W/(m·K)\n";
    out << "#    Cahill Model (Isotropic Profile)  = " << std::setw(6) << therm.kappa_cahill_iso << " W/(m·K)\n";
    out << "#    Cahill Model (Anisotropic Profile)= " << std::setw(6) << therm.kappa_cahill_aniso << " W/(m·K)\n";
    out << std::scientific << std::setprecision(3);
    out << "#  Volumetric Thermal Expansion α_V    = " << therm.alpha_V << " K⁻¹  (High-T limit)\n";
}

double calculateDensityFromLog(const ElasticDeformLog& log) {
    if (log.volume <= 0.0) {
        std::cerr << "Error: Cell volume in log file is zero or negative.\n";
        return 0.0;
    }
    if (log.elems.size() != log.num_atms.size() || log.elems.empty()) {
        std::cerr << "Error: Mismatched or empty elements/atom-counts vectors in log.\n";
        return 0.0;
    }

    double total_mass = 0.0;
    for (size_t i = 0; i < log.elems.size(); ++i) {
        double weight = getAtomicWeight(log.elems[i]);
        if (weight == 0.0) {
            std::cerr << "Warning: Atomic weight for element '" << log.elems[i]
                      << "' is unknown! Calculated density will be low.\n";
        }
        total_mass += log.num_atms[i] * weight;
    }

    // Mathematical Conversion Factor:
    // (Total Molar Mass / Avogadro's Number) / (Volume * 10^-24 cm^3 per Angstrom^3)
    // Simplifies to: (mass / volume) * 1.660539067
    constexpr double AMU_TO_GRAMS_FACTOR = 1.660539067;
    return (total_mass / log.volume) * AMU_TO_GRAMS_FACTOR;
}

int main(int argc, char* argv[]) {
    CLI::App app{"Calculate mechanical, acoustic, and thermal properties from elastic constants"};

    std::string inputFile{"elastic_constants.dat"};
    std::string logFile{"elastic_deform.log"};
    std::string outputFile{"thermal_properties.txt"};

    // Density and Volume are often needed for acoustic/thermal properties.
    // We can either parse them from the comment/header if available or let users override them.
    double density{0.0};  // g/cm^3
    double volume{0.0};   // Angstrom^3
    int total_atoms{0};

    app.add_option("--input,-i", inputFile, "File containing elastic constants table C_ij")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--log,-l", logFile, "Path to elastic_deform.log manifest")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--output,-o", outputFile, "Output file for calculated properties")->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    // -- Reading C_ij from input file -----
    Cmat C_ij{};
    ParseError error;

    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "I/O error while reading " << inputFile << ": cannot open file\n";
        return 1;
    }

    if (!fetchCij(file, C_ij, volume, error)) {
        reportParseError(inputFile, error);
        return 1;
    }

    // ── Read manifest ────────────────────────────────────────────────────────
    auto log_opt = readElasticLog(logFile);
    if (!log_opt) {
        std::cerr << "Error: cannot read manifest: " << logFile << "\n";
        return 1;
    }
    const ElasticDeformLog& log = *log_opt;

    density = calculateDensityFromLog(log);
    volume = log.volume;
    total_atoms = log.total_atms;

    if (density <= 0.0) {
        std::cerr << "Fatal Error: Invalid structure properties parsed. Cannot compute density.\n";
        return 1;
    }

    std::cout << "Successfully parsed structural metadata:\n";
    std::cout << "  Chemical Formula : ";
    for (size_t i = 0; i < log.elems.size(); ++i) {
        std::cout << log.elems[i] << log.num_atms[i] << " ";
    }
    std::cout << "\n  Total Atoms      : " << total_atoms << "\n";
    std::cout << "  Cell Volume      : " << volume << " Ang^3\n";
    std::cout << "  Calculated Density: " << std::fixed << std::setprecision(4) << density << " g/cm^3\n\n";

    std::ofstream fout(outputFile);
    if (!fout) {
        std::cerr << "Error: cannot open output file: " << outputFile << "\n";
        return 1;
    }

    // ── Optional VRH averages ────────────────────────────────────────────────
    VRHResult vrh;
    SoundVelocitiesResult vel;
    ThermalPropertiesResult therm;
    if (!computeVRH(C_ij, vrh)) {
        std::cerr << "Warning: C_ij matrix is singular; VRH averages CANNOT be computed!\n"
                  << "Warning: sound velocities and thermal properties CANNOT be computed!\n";
    } else {
        writeVRH(fout, vrh);
        std::cout << "VRH averages written to: " << outputFile << "\n";

        vel = compute_sound_velocities(C_ij, density, vrh);
        writeVelocitiesOutput(fout, vel);
        std::cout << "Sound velocities written to: " << outputFile << "\n";

        therm = compute_thermal_properties(vrh, vel, total_atoms, volume, density);
        writeThermalOutput(fout, therm);
        std::cout << "Thermal properties written to: " << outputFile << "\n";
    }

    return 0;
}
