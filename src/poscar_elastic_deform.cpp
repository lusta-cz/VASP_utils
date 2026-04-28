#include <CLI/CLI.hpp>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "elastic.h"
#include "poscar_file.h"
#include "symmetry.h"

namespace fs = std::filesystem;

namespace {

// Format a signed amplitude for directory names: +0.010000 / -0.010000
std::string fmtAmp(double amp) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    if (amp >= 0.0)
        oss << "+";
    oss << amp;
    return oss.str();
}

// Build the list of non-zero deformation amplitudes (symmetric about zero).
// npoints must be even; step = amplitude / (npoints/2).
std::vector<double> buildAmplitudes(double amplitude, int npoints) {
    const int n_half = npoints / 2;
    const double step = amplitude / n_half;
    std::vector<double> amps;
    amps.reserve(npoints);
    for (int i = -n_half; i <= n_half; ++i) {
        if (i != 0)
            amps.push_back(i * step);
    }
    return amps;
}

// Zero-pad mode index for directory names: 01, 02, ..., 21
std::string modePrefix(int idx, int total) {
    const int width = total >= 10 ? 2 : 1;
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << idx;
    return oss.str();
}

bool writeDeformedPOSCAR(const POSCAR& poscar, const std::array<double, 6>& voigt, double amplitude,
                         const fs::path& dir) {
    fs::create_directories(dir);
    POSCAR deformed = applyStrain(poscar, voigt, amplitude);
    return deformed.writePOSCAR((dir / "POSCAR").string());
}

}  // namespace

int main(int argc, char* argv[]) {
    CLI::App app{"Generate defomred POSCAR files for elastic constants calculation"};

    std::string inputFile{"POSCAR"};
    std::string outputDir;
    std::string method;
    std::string symmetryMode{"auto"};
    double amplitude{0.01};
    int npoints{-1};
    bool yes_flag{false};

    app.add_option("--input,-i", inputFile, "Input equilibrium POSCAR")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--output,-o", outputDir, "Root output directory")->required();
    app.add_option("--method,-m", method, "Post-processing method: energy or stress")
        ->required()
        ->check(CLI::IsMember({"energy", "stress"}));
    app.add_option("--symmetry,-s", symmetryMode,
                   "Symmetry handling: auto (detect via spglib) or none (triclinic, 21 constants)")
        ->capture_default_str()
        ->check(CLI::IsMember({"auto", "none"}));
    app.add_option("--amplitude,-a", amplitude,
                   "Maximum strain amplitude (dimensionless). Recommended: 0.01. Limit: 0.05.")
        ->capture_default_str()
        ->check(CLI::Range(1e-6, 0.10));
    app.add_option("--npoints,-n", npoints,
                   "Number of non-zero deformed structures per mode (must be even). "
                   "Default: 6 for energy, 4 for stress.");
    app.add_flag("--yes,-y", yes_flag, "Skip confirmation prompt for large calculation sets");

    CLI11_PARSE(app, argc, argv);

    // ── Apply defaults and validate npoints ──────────────────────────────────
    if (npoints == -1)
        npoints = (method == "energy") ? 7 : 5;

    if (npoints < 3) {
        std::cerr << "Error: --npoints must be at least 3\n";
        return 1;
    }
    if (npoints % 2 == 0) {
        std::cerr << "Warning: --npoints must be odd for a symmetric distribution. " << "Using " << (npoints + 1)
                  << " instead of " << npoints << ".\n";
        npoints += 1;
    }

    if (amplitude > 0.05) {
        std::cerr << "Warning: amplitude " << amplitude
                  << " exceeds the recommended limit of 0.05. Results may be unreliable.\n";
    }

    POSCAR poscar;
    if (!poscar.readPOSCAR(inputFile)) {
        std::cerr << "Error reading POSCAR file: " << inputFile << "\n";
        return 1;
    }

    // Poscar is easier to handel in direct coordinates
    if (!poscar.is_direct) {
        poscar.toDirect();
    }

    // ── Determine crystal system ─────────────────────────────────────────────
    CrystalSystem cs = CrystalSystem::Triclinic;
    int spg_number = 1;
    std::string spg_symbol = "P1";
    std::string point_group = "1";

    if (symmetryMode == "auto") {
        constexpr double symprec = 1e-5;
        SpglibDatasetPtr dataset = analyzeSymmetry(poscar, symprec);
        if (!dataset || dataset->spacegroup_number <= 0) {
            std::cerr << "Warning: spglib symmetry analysis failed. Falling back to triclinic.\n";
        } else {
            spg_number = dataset->spacegroup_number;
            spg_symbol = dataset->international_symbol;
            point_group = dataset->pointgroup_symbol;
            cs = crystalSystemFromSpaceGroup(spg_number);
        }
    }

    const int n_independent = nIndependentConstants(cs);

    // ── Select strain modes ──────────────────────────────────────────────────
    const std::vector<ElasticStrainMode> modes = (method == "energy") ? energyStrainModes(cs) : stressStrainModes();

    const std::vector<double> amps = buildAmplitudes(amplitude, npoints);
    const int total_calculations = static_cast<int>(modes.size()) * static_cast<int>(amps.size()) + 1;

    // ── Print summary and warn if large ─────────────────────────────────────
    std::cout << "=== Elastic Deformation Generator ===\n";
    std::cout << "Input:               " << inputFile << "\n";
    std::cout << "Method:              " << method << "-strain\n";
    std::cout << "Symmetry:            "
              << (symmetryMode == "none"
                      ? "none (triclinic override)"
                      : crystalSystemName(cs) + " (space group " + spg_symbol + " #" + std::to_string(spg_number) + ")")
              << "\n";
    std::cout << "Point group:         " << (symmetryMode == "none" ? "1" : point_group) << "\n";
    std::cout << "Independent C_ij:    " << n_independent << "\n";
    std::cout << "Strain modes:        " << modes.size() << "\n";
    std::cout << "Structures per mode: " << amps.size() << " deformed  +  1 shared reference\n";
    std::cout << "Amplitude range:     [" << fmtAmp(-amplitude) << ", " << fmtAmp(amplitude) << "]\n";
    std::cout << "Total calculations:  " << total_calculations << "\n";

    if (total_calculations > 50) {
        std::cout << "\nWARNING: " << total_calculations << " VASP calculations required.\n";
        if (method == "energy" && symmetryMode == "none") {
            std::cout << "         Consider --method stress (" << 6 * npoints + 1
                      << " calculations) or --symmetry auto to reduce the count.\n";
        }
        if (!yes_flag) {
            std::cout << "Proceed? [y/N]: ";
            std::string answer;
            std::getline(std::cin, answer);
            if (answer != "y" && answer != "Y") {
                std::cout << "Aborted.\n";
                return 0;
            }
        }
    }

    // ── Create reference directory ───────────────────────────────────────────
    const fs::path root(outputDir);
    const fs::path ref_dir = root / "reference";
    fs::create_directories(ref_dir);
    if (!poscar.writePOSCAR((ref_dir / "POSCAR").string())) {
        std::cerr << "Error: cannot write reference POSCAR\n";
        return 1;
    }

    // ── Generate deformed structures ─────────────────────────────────────────
    ElasticDeformLog log;
    log.method = method;
    log.amplitude = amplitude;
    log.npoints = npoints;
    log.symmetry_mode = symmetryMode;
    log.crystal_system = cs;
    log.space_group_number = spg_number;
    log.space_group_symbol = spg_symbol;
    log.point_group = point_group;
    log.n_independent = n_independent;
    log.reference_dir = "reference";
    log.modes = modes;
    log.amplitudes.resize(modes.size());
    log.dirs.resize(modes.size());

    const int n_modes = static_cast<int>(modes.size());
    for (int m = 0; m < n_modes; ++m) {
        const auto& mode = modes[m];
        const std::string mode_dir_name = "mode_" + modePrefix(m + 1, n_modes) + "_" + mode.label;
        const fs::path mode_dir = root / mode_dir_name;

        for (double amp : amps) {
            const std::string amp_dir_name = "amp_" + fmtAmp(amp);
            const fs::path amp_dir = mode_dir / amp_dir_name;

            if (!writeDeformedPOSCAR(poscar, mode.voigt, amp, amp_dir)) {
                std::cerr << "Error: cannot write POSCAR to " << amp_dir << "\n";
                return 1;
            }

            log.amplitudes[m].push_back(amp);
            log.dirs[m].push_back((fs::path(mode_dir_name) / amp_dir_name).string());
        }
    }

    // ── Write manifest ───────────────────────────────────────────────────────
    const std::string log_path = (root / "elastic_deform.log").string();
    if (!writeElasticLog(log_path, log)) {
        std::cerr << "Error: cannot write manifest: " << log_path << "\n";
        return 1;
    }

    std::cout << "\nOutput directory: " << fs::absolute(root).string() << "\n";
    std::cout << "Manifest:         " << log_path << "\n";
    std::cout << "Done. Run VASP in each subdirectory, then use vasp_elastic to extract C_ij.\n";

    return 0;
}
