#include "elastic.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "poscar_file.h"

// ─── Crystal system classification ──────────────────────────────────────────

CrystalSystem crystalSystemFromSpaceGroup(int n) {
    if (n >= 1 && n <= 2)
        return CrystalSystem::Triclinic;
    if (n >= 3 && n <= 15)
        return CrystalSystem::Monoclinic;
    if (n >= 16 && n <= 74)
        return CrystalSystem::Orthorhombic;
    if (n >= 75 && n <= 88)
        return CrystalSystem::TetragonalI;
    if (n >= 89 && n <= 142)
        return CrystalSystem::TetragonalII;
    if (n >= 143 && n <= 148)
        return CrystalSystem::TrigonalI;
    if (n >= 149 && n <= 167)
        return CrystalSystem::TrigonalII;
    if (n >= 168 && n <= 194)
        return CrystalSystem::Hexagonal;
    if (n >= 195 && n <= 230)
        return CrystalSystem::Cubic;
    throw std::out_of_range("Space group number out of range 1–230");
}

std::string crystalSystemName(CrystalSystem cs) {
    switch (cs) {
        case CrystalSystem::Triclinic:
            return "triclinic";
        case CrystalSystem::Monoclinic:
            return "monoclinic";
        case CrystalSystem::Orthorhombic:
            return "orthorhombic";
        case CrystalSystem::TetragonalI:
            return "tetragonal-I (4, -4, 4/m)";
        case CrystalSystem::TetragonalII:
            return "tetragonal-II (422, 4mm, -42m, 4/mmm)";
        case CrystalSystem::TrigonalI:
            return "trigonal-I (3, -3)";
        case CrystalSystem::TrigonalII:
            return "trigonal-II (32, 3m, -3m)";
        case CrystalSystem::Hexagonal:
            return "hexagonal";
        case CrystalSystem::Cubic:
            return "cubic";
    }
    return "unknown";
}

int nIndependentConstants(CrystalSystem cs) {
    switch (cs) {
        case CrystalSystem::Triclinic:
            return 21;
        case CrystalSystem::Monoclinic:
            return 13;
        case CrystalSystem::Orthorhombic:
            return 9;
        case CrystalSystem::TetragonalI:
            return 7;
        case CrystalSystem::TetragonalII:
            return 6;
        case CrystalSystem::TrigonalI:
            return 7;
        case CrystalSystem::TrigonalII:
            return 6;
        case CrystalSystem::Hexagonal:
            return 5;
        case CrystalSystem::Cubic:
            return 3;
    }
    return 21;
}

// ─── Strain mode tables ──────────────────────────────────────────────────────

std::vector<ElasticStrainMode> energyStrainModes(CrystalSystem cs) {
    // Voigt vector layout: [ε1(xx), ε2(yy), ε3(zz), ε4(2yz), ε5(2xz), ε6(2xy)]
    // Strategy: pure modes give diagonal Cii; mixed modes give Cii+Cjj+2Cij.
    // Reference: Le Page & Saxe, PRB 65, 104104 (2002).

    switch (cs) {
        case CrystalSystem::Cubic:
            // C11, C12, C44  (C11=C22=C33, C12=C13=C23, C44=C55=C66)
            return {
                {"vol", {1, 1, 1, 0, 0, 0}, false},    // E ~ (C11+2C12)/2
                {"e1-e2", {1, -1, 0, 0, 0, 0}, true},  // E ~ (C11-C12)
                {"e4", {0, 0, 0, 1, 0, 0}, true},      // E ~ C44/2
            };

        case CrystalSystem::Hexagonal:
            // C11, C12, C13, C33, C44  (C66=(C11-C12)/2 by symmetry)
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},       {"e3", {0, 0, 1, 0, 0, 0}, false},
                {"e1-e2", {1, -1, 0, 0, 0, 0}, true},     // → C11-C12
                {"e1+e2+e3", {1, 1, 1, 0, 0, 0}, false},  // → C11+C12+C33+2C13
                {"e4", {0, 0, 0, 1, 0, 0}, true},
            };

        case CrystalSystem::TrigonalII:
            // C11, C12, C13, C33, C44, C14  (C66=(C11-C12)/2)
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e3", {0, 0, 1, 0, 0, 0}, false},
                {"e1+e2", {1, 1, 0, 0, 0, 0}, false}, {"e1+e2+e3", {1, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, false},    {"e1+e4", {1, 0, 0, 1, 0, 0}, false},  // → C11+C44+2C14
            };

        case CrystalSystem::TrigonalI:
            // C11, C12, C13, C33, C44, C14, C25  (C66=(C11-C12)/2)
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e3", {0, 0, 1, 0, 0, 0}, false},
                {"e1+e2", {1, 1, 0, 0, 0, 0}, false}, {"e1+e2+e3", {1, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, false},    {"e1+e4", {1, 0, 0, 1, 0, 0}, false},
                {"e2+e5", {0, 1, 0, 0, 1, 0}, false},  // → C22+C55+2C25
            };

        case CrystalSystem::TetragonalII:
            // C11, C12, C13, C33, C44, C66  (C11=C22, C44=C55)
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e3", {0, 0, 1, 0, 0, 0}, false},
                {"e1-e2", {1, -1, 0, 0, 0, 0}, true}, {"e1+e2+e3", {1, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, true},     {"e6", {0, 0, 0, 0, 0, 1}, true},
            };

        case CrystalSystem::TetragonalI:
            // C11, C12, C13, C33, C44, C66, C16
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e3", {0, 0, 1, 0, 0, 0}, false},
                {"e1-e2", {1, -1, 0, 0, 0, 0}, true}, {"e1+e2+e3", {1, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, true},     {"e6", {0, 0, 0, 0, 0, 1}, true},
                {"e1+e6", {1, 0, 0, 0, 0, 1}, false},  // → C11+C66+2C16
            };

        case CrystalSystem::Orthorhombic:
            // C11, C22, C33, C12, C13, C23, C44, C55, C66
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e2", {0, 1, 0, 0, 0, 0}, false},
                {"e3", {0, 0, 1, 0, 0, 0}, false},    {"e1+e2", {1, 1, 0, 0, 0, 0}, false},
                {"e1+e3", {1, 0, 1, 0, 0, 0}, false}, {"e2+e3", {0, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, true},     {"e5", {0, 0, 0, 0, 1, 0}, true},
                {"e6", {0, 0, 0, 0, 0, 1}, true},
            };

        case CrystalSystem::Monoclinic:
            // 13 constants (monoclinic-b convention):
            // C11,C22,C33,C44,C55,C66,C12,C13,C23,C15,C25,C35,C46
            return {
                {"e1", {1, 0, 0, 0, 0, 0}, false},    {"e2", {0, 1, 0, 0, 0, 0}, false},
                {"e3", {0, 0, 1, 0, 0, 0}, false},    {"e1+e2", {1, 1, 0, 0, 0, 0}, false},
                {"e1+e3", {1, 0, 1, 0, 0, 0}, false}, {"e2+e3", {0, 1, 1, 0, 0, 0}, false},
                {"e4", {0, 0, 0, 1, 0, 0}, false},    {"e5", {0, 0, 0, 0, 1, 0}, false},
                {"e6", {0, 0, 0, 0, 0, 1}, false},    {"e1+e5", {1, 0, 0, 0, 1, 0}, false},  // → C11+C55+2C15
                {"e2+e5", {0, 1, 0, 0, 1, 0}, false},                                        // → C22+C55+2C25
                {"e3+e5", {0, 0, 1, 0, 1, 0}, false},                                        // → C33+C55+2C35
                {"e4+e6", {0, 0, 0, 1, 0, 1}, false},                                        // → C44+C66+2C46
            };

        case CrystalSystem::Triclinic:
        default: {
            // 21 constants: 6 pure + 15 mixed (all unique pairs i<j)
            std::vector<ElasticStrainMode> modes;
            // 6 pure modes
            const char* pure_labels[6] = {"e1", "e2", "e3", "e4", "e5", "e6"};
            for (int i = 0; i < 6; ++i) {
                std::array<double, 6> v{};
                v[i] = 1.0;
                modes.push_back({pure_labels[i], v, false});
            }
            // 15 mixed modes
            for (int i = 0; i < 6; ++i) {
                for (int j = i + 1; j < 6; ++j) {
                    std::array<double, 6> v{};
                    v[i] = 1.0;
                    v[j] = 1.0;
                    modes.push_back({"e" + std::to_string(i + 1) + "+e" + std::to_string(j + 1), v, false});
                }
            }
            return modes;
        }
    }
}

std::vector<ElasticStrainMode> stressStrainModes() {
    return {
        {"e1", {1, 0, 0, 0, 0, 0}}, {"e2", {0, 1, 0, 0, 0, 0}}, {"e3", {0, 0, 1, 0, 0, 0}},
        {"e4", {0, 0, 0, 1, 0, 0}}, {"e5", {0, 0, 0, 0, 1, 0}}, {"e6", {0, 0, 0, 0, 0, 1}},
    };
}

// ─── Strain application ──────────────────────────────────────────────────────

POSCAR applyStrain(const POSCAR& poscar, const std::array<double, 6>& voigt, double amplitude) {
    // Build symmetric strain tensor ε[i][j] from Voigt vector scaled by amplitude.
    // Voigt: ε4 = 2*ε_yz, ε5 = 2*ε_xz, ε6 = 2*ε_xy → off-diagonal = voigt/2.
    double eps[3][3] = {};
    eps[0][0] = voigt[0] * amplitude;
    eps[1][1] = voigt[1] * amplitude;
    eps[2][2] = voigt[2] * amplitude;
    eps[1][2] = eps[2][1] = voigt[3] * amplitude / 2.0;
    eps[0][2] = eps[2][0] = voigt[4] * amplitude / 2.0;
    eps[0][1] = eps[1][0] = voigt[5] * amplitude / 2.0;

    // Deformation matrix F = I + ε (symmetric, no rotation).
    double F[3][3] = {};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j ? 1.0 : 0.0) + eps[i][j];
        }
    }

    // Transform lattice row vectors: new_L[i] = old_L[i] · F  (F symmetric → F = F^T).
    // new_lattice[i][j] = sum_k old_lattice[i][k] * F[k][j]
    POSCAR out = poscar;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out.lattice[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                out.lattice[i][j] += poscar.lattice[i][k] * F[k][j];
            }
        }
    }
    return out;
}

// ─── Log file I/O ────────────────────────────────────────────────────────────

namespace {

std::string formatAmplitude(double amp) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    if (amp >= 0.0)
        oss << "+";
    oss << amp;
    return oss.str();
}

std::string crystalSystemKey(CrystalSystem cs) {
    switch (cs) {
        case CrystalSystem::Triclinic:
            return "triclinic";
        case CrystalSystem::Monoclinic:
            return "monoclinic";
        case CrystalSystem::Orthorhombic:
            return "orthorhombic";
        case CrystalSystem::TetragonalI:
            return "tetragonal_I";
        case CrystalSystem::TetragonalII:
            return "tetragonal_II";
        case CrystalSystem::TrigonalI:
            return "trigonal_I";
        case CrystalSystem::TrigonalII:
            return "trigonal_II";
        case CrystalSystem::Hexagonal:
            return "hexagonal";
        case CrystalSystem::Cubic:
            return "cubic";
    }
    return "triclinic";
}

CrystalSystem crystalSystemFromKey(const std::string& key) {
    if (key == "monoclinic")
        return CrystalSystem::Monoclinic;
    if (key == "orthorhombic")
        return CrystalSystem::Orthorhombic;
    if (key == "tetragonal_I")
        return CrystalSystem::TetragonalI;
    if (key == "tetragonal_II")
        return CrystalSystem::TetragonalII;
    if (key == "trigonal_I")
        return CrystalSystem::TrigonalI;
    if (key == "trigonal_II")
        return CrystalSystem::TrigonalII;
    if (key == "hexagonal")
        return CrystalSystem::Hexagonal;
    if (key == "cubic")
        return CrystalSystem::Cubic;
    return CrystalSystem::Triclinic;
}

}  // namespace

bool writeElasticLog(const std::string& path, const ElasticDeformLog& log) {
    std::ofstream f(path);
    if (!f) {
        std::cerr << "Error: cannot open log file for writing: " << path << "\n";
        return false;
    }

    f << std::fixed << std::setprecision(6);
    f << "METHOD=" << log.method << "\n";
    f << "AMPLITUDE=" << log.amplitude << "\n";
    f << "NPOINTS=" << log.npoints << "\n";
    f << "SYMMETRY_MODE=" << log.symmetry_mode << "\n";
    f << "CRYSTAL_SYSTEM=" << crystalSystemKey(log.crystal_system) << "\n";
    f << "SPACE_GROUP_NUMBER=" << log.space_group_number << "\n";
    f << "SPACE_GROUP_SYMBOL=" << log.space_group_symbol << "\n";
    f << "POINT_GROUP=" << log.point_group << "\n";
    f << "N_INDEPENDENT=" << log.n_independent << "\n";
    f << "VOLUME=" << log.volume << "\n";
    f << "N_MODES=" << static_cast<int>(log.modes.size()) << "\n";
    f << "REFERENCE_DIR=" << log.reference_dir << "\n";

    for (size_t m = 0; m < log.modes.size(); ++m) {
        const auto& mode = log.modes[m];
        const auto& amps = log.amplitudes[m];
        const auto& dirs = log.dirs[m];

        f << "MODE_" << (m + 1) << "_LABEL=" << mode.label << "\n";

        f << "MODE_" << (m + 1) << "_VOIGT=";
        for (int k = 0; k < 6; ++k)
            f << mode.voigt[k] << (k < 5 ? " " : "");
        f << "\n";

        f << "MODE_" << (m + 1) << "_AMPLITUDES=";
        for (size_t a = 0; a < amps.size(); ++a)
            f << formatAmplitude(amps[a]) << (a + 1 < amps.size() ? " " : "");
        f << "\n";

        f << "MODE_" << (m + 1) << "_DIRS=";
        for (size_t a = 0; a < dirs.size(); ++a)
            f << dirs[a] << (a + 1 < dirs.size() ? " " : "");
        f << "\n";
    }

    return f.good();
}

std::optional<ElasticDeformLog> readElasticLog(const std::string& path) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Error: cannot open log file: " << path << "\n";
        return std::nullopt;
    }

    ElasticDeformLog log;
    int n_modes = 0;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos)
            continue;
        const std::string key = line.substr(0, eq);
        const std::string val = line.substr(eq + 1);

        if (key == "METHOD")
            log.method = val;
        else if (key == "AMPLITUDE")
            log.amplitude = std::stod(val);
        else if (key == "NPOINTS")
            log.npoints = std::stoi(val);
        else if (key == "SYMMETRY_MODE")
            log.symmetry_mode = val;
        else if (key == "CRYSTAL_SYSTEM")
            log.crystal_system = crystalSystemFromKey(val);
        else if (key == "SPACE_GROUP_NUMBER")
            log.space_group_number = std::stoi(val);
        else if (key == "SPACE_GROUP_SYMBOL")
            log.space_group_symbol = val;
        else if (key == "POINT_GROUP")
            log.point_group = val;
        else if (key == "N_INDEPENDENT")
            log.n_independent = std::stoi(val);
        else if (key == "VOLUME")
            log.volume = std::stod(val);
        else if (key == "REFERENCE_DIR")
            log.reference_dir = val;
        else if (key == "N_MODES") {
            n_modes = std::stoi(val);
            log.modes.resize(n_modes);
            log.amplitudes.resize(n_modes);
            log.dirs.resize(n_modes);
        } else {
            // MODE_N_LABEL, MODE_N_VOIGT, MODE_N_AMPLITUDES, MODE_N_DIRS
            for (int m = 0; m < n_modes; ++m) {
                const std::string prefix = "MODE_" + std::to_string(m + 1) + "_";
                if (key == prefix + "LABEL") {
                    log.modes[m].label = val;
                } else if (key == prefix + "VOIGT") {
                    std::istringstream ss(val);
                    for (int k = 0; k < 6; ++k)
                        ss >> log.modes[m].voigt[k];
                } else if (key == prefix + "AMPLITUDES") {
                    std::istringstream ss(val);
                    std::string tok;
                    while (ss >> tok)
                        log.amplitudes[m].push_back(std::stod(tok));
                } else if (key == prefix + "DIRS") {
                    std::istringstream ss(val);
                    std::string tok;
                    while (ss >> tok)
                        log.dirs[m].push_back(tok);
                }
            }
        }
    }

    if (log.method.empty() || n_modes == 0) {
        std::cerr << "Error: log file missing required fields\n";
        return std::nullopt;
    }
    return log;
}
