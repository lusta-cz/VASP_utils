#include <algorithm>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error_parse.h"

namespace fs = std::filesystem;

// Structure to pair a parsed test parameter with its extracted final energy value and grid
struct DataPoint {
    double parameter{0.0};
    double energy{0.0};
    std::string kgrid{"N/A"};  // Stores the extracted grid dimensions (e.g., "8x8x8")
};

// Helper function to trim whitespace padding from strings
std::string trimString(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

// Scrapes the final converged total free energy and the automatically generated k-grid from an OUTCAR file
bool extractOutcarData(const std::string& outcar_path, double& energy, std::string& kgrid, ParseError& error) {
    std::ifstream file(outcar_path);
    if (!file) {
        return fail(error, ParseErrorKind::Io, 0, "Cannot open file " + outcar_path);
    }

    std::string line;
    bool found_energy = false;
    double last_energy = 0.0;
    int line_num = 0;
    kgrid = "N/A";

    while (std::getline(file, line)) {
        ++line_num;

        // VASP OUTCAR k-grid configuration line identifier: "Successfully generated k-point grid:"
        if (line.find("Successfully generated k-point grid:") != std::string::npos) {
            size_t pos = line.find(":");
            if (pos != std::string::npos) {
                std::string grid_raw = line.substr(pos + 1);
                std::stringstream ss(grid_raw);
                int k1, k2, k3;
                if (ss >> k1 >> k2 >> k3) {
                    kgrid = std::to_string(k1) + "x" + std::to_string(k2) + "x" + std::to_string(k3);
                }
            }
        }

        // VASP OUTCAR energy benchmark line identifier: "energy  without entropy="
        if (line.find("energy  without entropy=") != std::string::npos) {
            size_t energy_pos = line.find_last_of(" \t");
            if (energy_pos != std::string::npos) {
                try {
                    last_energy = std::stod(trimString(line.substr(energy_pos)));
                    found_energy = true;
                } catch (...) {
                    continue;  // Fail gracefully if a partial line format error occurs during active runs
                }
            }
        }
    }

    if (!found_energy) {
        return fail(error, ParseErrorKind::Parse, line_num,
                    "Could not locate completed free energy blocks in " + outcar_path);
    }

    energy = last_energy;
    return true;
}

// Scans a parent directory to process nested test calculations
void processConvergenceDirectory(const std::string& parent_dir, const std::string& prefix,
                                 const std::string& output_filename) {
    if (!fs::exists(parent_dir)) {
        std::cout << " -> Notice: Parent matrix track '" << parent_dir
                  << "' not found. Skipping calculation collection.\n";
        return;
    }

    std::vector<DataPoint> dataset;
    ParseError error;

    std::cout << "Parsing simulation paths inside '" << parent_dir << "'...\n";

    for (const auto& entry : fs::directory_iterator(parent_dir)) {
        if (!entry.is_directory())
            continue;

        std::string folder_name = entry.path().filename().string();
        if (folder_name.rfind(prefix, 0) == 0) {
            std::string val_str = folder_name.substr(prefix.length());
            double param_val = 0.0;
            try {
                param_val = std::stod(val_str);
            } catch (...) {
                std::cerr << " -> Warning: Could not parse numerical parameter from path: " << folder_name << "\n";
                continue;
            }

            fs::path outcar_file = entry.path() / "OUTCAR";
            if (!fs::exists(outcar_file)) {
                std::cerr << " -> Warning: Missing OUTCAR in calculation folder: " << folder_name
                          << " (Job may still be running)\n";
                continue;
            }

            double final_energy = 0.0;
            std::string generated_kgrid = "N/A";
            if (extractOutcarData(outcar_file.string(), final_energy, generated_kgrid, error)) {
                dataset.push_back({param_val, final_energy, generated_kgrid});
            } else {
                reportParseError(outcar_file.string(), error);
            }
        }
    }

    if (dataset.empty()) {
        std::cerr << " -> Error: No valid completed energy sets extracted from '" << parent_dir << "'.\n";
        return;
    }

    // Sort dataset rows cleanly
    if (prefix == "encut_") {
        std::sort(dataset.begin(), dataset.end(),
                  [](const DataPoint& a, const DataPoint& b) { return a.parameter < b.parameter; });
    } else {
        std::sort(dataset.begin(), dataset.end(),
                  [](const DataPoint& a, const DataPoint& b) { return a.parameter > b.parameter; });
    }

    // Write matrix values cleanly out to a space/tab-separated data file
    std::ofstream out(output_filename);
    if (!out) {
        std::cerr << "I/O Error: Cannot write results matrix to target file: " << output_filename << "\n";
        return;
    }

    if (prefix == "encut_") {
        out << "# ENCUT (eV)\tTotal_Energy (eV)\tGenerated_KGrid\n";
        for (const auto& pt : dataset) {
            out << std::setw(8) << static_cast<int>(pt.parameter) << "\t" << std::fixed << std::setprecision(6)
                << pt.energy << "\t" << pt.kgrid << "\n";
        }
    } else {
        out << "# KSPACING (1/A)\tTotal_Energy (eV)\tGenerated_KGrid\n";
        for (const auto& pt : dataset) {
            out << std::fixed << std::setprecision(4) << pt.parameter << "\t\t" << std::fixed << std::setprecision(6)
                << pt.energy << "\t" << pt.kgrid << "\n";
        }
    }

    std::cout << "Successfully exported collected profiles to: " << output_filename << "\n";
}

int main(int argc, char* argv[]) {
    CLI::App app{"Post-processing utility to extract energy profiles from convergence test subdirectories"};

    std::string encut_out{"encut_convergence.dat"};
    std::string kspace_out{"kspacing_convergence.dat"};

    app.add_option("--output-encut,-e", encut_out, "Output filename for ENCUT data profile matrix")
        ->capture_default_str();
    app.add_option("--output-kspacing,-k", kspace_out, "Output filename for KSPACING data profile matrix")
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    std::cout << "Starting data matrix extraction from completed directories...\n\n";

    // Parse the ENCUT calculations directory tree
    processConvergenceDirectory("convergence_encut", "encut_", encut_out);
    std::cout << "\n-----------------------------------------------------------------------\n\n";
    // Parse the KSPACING calculations directory tree
    processConvergenceDirectory("convergence_kspacing", "kspacing_", kspace_out);

    std::cout << "\nData processing run completed.\n";
    return 0;
}
