#include <CLI/CLI.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error_parse.h"

namespace fs = std::filesystem;

// Structure to hold optimization flags parsed from INCAR
struct IncarStatus {
    int ibrion{-1};  // Default in VASP when not specified is -1 (Static)
    int nsw{0};      // Default in VASP when not specified is 0
};

// Helper function to trim whitespaces from string views/strings
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// Scans INCAR to determine if ionic relaxation geometry optimization is on
bool checkGeometryOptimization(const std::string& incar_path, IncarStatus& status) {
    std::ifstream file(incar_path);
    if (!file)
        return false;

    std::string line;
    while (std::getline(file, line)) {
        // Strip out comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);
        comment_pos = line.find(';');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);

        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos)
            continue;

        std::string key = trim(line.substr(0, eq_pos));
        std::string val = trim(line.substr(eq_pos + 1));

        if (key == "IBRION") {
            try {
                status.ibrion = std::stoi(val);
            } catch (...) {
            }
        } else if (key == "NSW") {
            try {
                status.nsw = std::stoi(val);
            } catch (...) {
            }
        }
    }
    return true;
}

// Safely creates a modified INCAR injecting or changing specific target keys
bool writeModifiedIncar(const std::string& source, const std::string& dest, const std::string& key1,
                        const std::string& val1, const std::string& key2, const std::string& val2) {
    std::ifstream in(source);
    std::ofstream out(dest);
    if (!in || !out)
        return false;

    std::string line;
    bool found_key1 = false;
    bool found_key2 = false;

    while (std::getline(in, line)) {
        std::string clean_line = line;
        size_t comment_pos = clean_line.find_first_of("#;");
        if (comment_pos != std::string::npos)
            clean_line = clean_line.substr(0, comment_pos);

        size_t eq_pos = clean_line.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = trim(clean_line.substr(0, eq_pos));
            if (key == key1) {
                out << key1 << " = " << val1 << "\n";
                found_key1 = true;
                continue;
            }
            if (key == key2) {
                out << key2 << " = " << val2 << "\n";
                found_key2 = true;
                continue;
            }
        }
        out << line << "\n";
    }

    // If flags weren't originally in the template INCAR, append them safely
    if (!found_key1)
        out << key1 << " = " << val1 << "\n";
    if (!found_key2)
        out << key2 << " = " << val2 << "\n";

    return true;
}

// Generic utility to copy optional secondary runtime files if present
void copyOptionalVaspFiles(const fs::path& target_dir) {
    const std::vector<std::string> optional_files = {"POSCAR", "POTCAR"};
    for (const auto& f : optional_files) {
        if (fs::exists(f)) {
            fs::copy_file(f, target_dir / f, fs::copy_options::overwrite_existing);
        }
    }
}

int main(int argc, char* argv[]) {
    CLI::App app{"Pre-processing utility for generating VASP convergence test directories"};

    // Setup Parameters with specified boundaries and defaults
    double encut_min{250.0}, encut_max{550.0}, encut_step{50.0};
    double kspace_min{0.20}, kspace_max{0.05}, kspace_step{-0.025};
    bool force_relax{false};

    app.add_option("--encut-min", encut_min, "Minimum ENCUT value for sweep")->capture_default_str();
    app.add_option("--encut-max", encut_max, "Maximum ENCUT value for sweep")->capture_default_str();
    app.add_option("--encut-step", encut_step, "Step increase increment for ENCUT")->capture_default_str();

    app.add_option("--kspacing-min", kspace_min, "Minimum numerical KSPACING (Coarse start)")->capture_default_str();
    app.add_option("--kspacing-max", kspace_max, "Maximum numerical KSPACING (Fine finish)")->capture_default_str();
    app.add_option("--kspacing-step", kspace_step, "Step decrement change size for KSPACING")->capture_default_str();

    app.add_flag("--force-relax", force_relax, "Override geometry optimization block safety warning");

    CLI11_PARSE(app, argc, argv);

    // 1. Strict input validation parameters checks
    if (encut_min < 250.0 || encut_max > 900.0 || encut_min > encut_max) {
        std::cerr << "Semantic Error: ENCUT parameters out of strict bounds [250 to 900 eV].\n";
        return 1;
    }
    // Note: Lower KSPACING value means denser grid, threshold cannot be lower than 0.009
    if (kspace_max < 0.009 || kspace_min > 0.5 || kspace_max > kspace_min) {
        std::cerr << "Semantic Error: KSPACING parameters out of strict bounds [0.009 to 0.5 A^-1].\n";
        return 1;
    }

    if (!fs::exists("INCAR")) {
        std::cerr << "I/O error: Baseline 'INCAR' template file not found in current directory.\n";
        return 1;
    }

    // 2. Optimization Status Checking Logic
    IncarStatus status;
    checkGeometryOptimization("INCAR", status);
    if (status.ibrion != -1 && status.nsw > 0) {
        if (!force_relax) {
            std::cerr << "=========================================================================\n"
                      << "WARNING: Active geometry relaxation detected (IBRION = " << status.ibrion
                      << ", NSW = " << status.nsw << ").\n"
                      << "Running convergence sweeps with relaxation enabled can be very expensive\n"
                      << "if the initial structure is far from minima.\n"
                      << "Please disable optimization inside INCAR or pass '--force-relax' to proceed.\n"
                      << "=========================================================================\n";
            return 1;
        }
    }

    std::cout << "Starting initialization of convergence environment testing tree matrices...\n";

    // 3. Generate ENCUT Test Subdirectories Matrix Block
    fs::path encut_parent("convergence_encut");
    fs::create_directories(encut_parent);

    // Generate exactly 7 steps or follow user loops
    int encut_calculations_count = 0;
    for (double e = encut_min; e <= encut_max + 1e-3; e += encut_step) {
        std::stringstream ss;
        ss << "encut_" << static_cast<int>(e);
        fs::path sub = encut_parent / ss.str();
        fs::create_directories(sub);

        // Lock companion KSPACING to a sound fallback parameter (0.15) during ENCUT testing
        if (!writeModifiedIncar("INCAR", (sub / "INCAR").string(), "ENCUT", std::to_string(static_cast<int>(e)),
                                "KSPACING", "0.15")) {
            std::cerr << "Error writing modifications to " << sub / "INCAR" << "\n";
            return 1;
        }
        copyOptionalVaspFiles(sub);
        encut_calculations_count++;
    }
    std::cout << " -> Generated " << encut_calculations_count << " test branches inside 'convergence_encut/'.\n";

    // 4. Generate KSPACING Test Subdirectories Matrix Block
    fs::path kspace_parent("convergence_kspacing");
    fs::create_directories(kspace_parent);

    int kspace_calculations_count = 0;
    // Step loop goes down from coarse numerical min value to tight max accuracy boundary
    for (double k = kspace_min; k >= kspace_max - 1e-5; k += kspace_step) {
        std::stringstream folder_name, val_stream;
        folder_name << "kspacing_" << std::fixed << std::setprecision(4) << k;
        val_stream << std::fixed << std::setprecision(4) << k;

        fs::path sub = kspace_parent / folder_name.str();
        fs::create_directories(sub);

        // Lock companion ENCUT to a safe sound fallback baseline parameter (500 eV) during KSPACING testing
        if (!writeModifiedIncar("INCAR", (sub / "INCAR").string(), "KSPACING", val_stream.str(), "ENCUT", "500")) {
            std::cerr << "Error writing modifications to " << sub / "INCAR" << "\n";
            return 1;
        }
        copyOptionalVaspFiles(sub);
        kspace_calculations_count++;
    }
    std::cout << " -> Generated " << kspace_calculations_count << " test branches inside 'convergence_kspacing/'.\n";

    std::cout << "\nSetup complete! You can now launch calculations across directory matrices.\n";
    return 0;
}
