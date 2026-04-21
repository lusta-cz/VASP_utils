#include <CLI/CLI.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "eigenval2bands.h"

int main(int argc, char** argv) {
    CLI::App app{"VASP EIGENVAL to Band Structure Utility"};

    // Configuration and default values
    std::string eigen_file{"EIGENVAL"};
    std::string kpoints_file{"KPOINTS"};
    std::string doscar_file{"DOSCAR"};
    std::string outcar_file{"OUTCAR"};
    double manual_fermi{0.0};
    std::string fermiSource{"doscar"};

    // CLI Options
    app.add_option("-i,--input", eigen_file, "Input EIGENVAL file")->default_val("EIGENVAL")->check(CLI::ExistingFile);
    app.add_option("-k,--kpoints", kpoints_file, "KPOINTS file used for the band structure calculation")
        ->default_val("KPOINTS")
        ->check(CLI::ExistingFile);
    app.add_option("--fermi-source", fermiSource, "Fermi level source: doscar (default) | outcar | manual | none")
        ->capture_default_str()
        ->check(CLI::IsMember({"doscar", "outcar", "manual", "none"}));
    app.add_option("--doscar", doscar_file, "DOSCAR file path (used with --fermi-source=doscar)")
        ->capture_default_str();
    app.add_option("--outcar", outcar_file, "OUTCAR file path (used with --fermi-source=outcar)")
        ->capture_default_str();
    auto* fermi_opt = app.add_option("-e,--fermi", manual_fermi, "Fermi level in eV (used with --fermi-source=manual)")
                          ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    double e_fermi{0.0};
    bool success{true};

    // Logic to determine the active Fermi source
    if (fermiSource == "doscar") {
        if (!parseFromDoscar(doscar_file, e_fermi)) {
            success = false;
        } else {
            std::cout << "Using Fermi level from " << doscar_file << ": " << e_fermi << " eV\n";
        }
    } else if (fermiSource == "outcar") {
        if (!parseFromOutcar(outcar_file, e_fermi)) {
            success = false;
        } else {
            std::cout << "Using Fermi level from " << outcar_file << ": " << e_fermi << " eV\n";
        }
    } else if (fermiSource == "manual") {
        if (fermi_opt->count() == 0) {
            std::cerr << "Error: --fermi-source=manual requires -e/--fermi VALUE\n";
            return 1;
        }
        e_fermi = manual_fermi;
        std::cout << "Info: Fermi energy is set to: " << e_fermi << " eV\n";
    } else {
        std::cout << "Warning: Energy shift is disabled. 0.0 eV is NOT the Fermi level!\n";
    }

    if (!success) {
        std::cout << "Warning: Reading of the Fermi level failed! Continuing without energy shift!\n";
    }

    int kpoints_between{0};
    // --- Parsing from KPOINTS file
    if (!parseKpoints(kpoints_file, kpoints_between)) {
        return 1;
    }

    EigenvalData data;
    // --- Parsing from EIGENVAL file
    if (!parseEigenval(eigen_file, data)) {
        return 1;
    }

    // --- Path Calculation and Output Generation ---
    if (data.total_kpoints == 0 || kpoints_between == 0) {
        std::cerr << "Error: Something went wrong while reading the number of k points!\n";
        return 1;
    }

    int num_paths = data.total_kpoints / kpoints_between;
    double cumulative_dist = 0.0;

    for (int p = 0; p < num_paths; ++p) {
        std::string out_name = "band_path_" + std::to_string(p + 1) + ".dat";
        std::ofstream out(out_name);
        out << std::fixed << std::setprecision(10);

        int s_idx = p * kpoints_between;
        int e_idx = s_idx + kpoints_between - 1;

        double dx = data.kpoints[e_idx].x - data.kpoints[s_idx].x;
        double dy = data.kpoints[e_idx].y - data.kpoints[s_idx].y;
        double dz = data.kpoints[e_idx].z - data.kpoints[s_idx].z;
        double segment_dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        double step = (kpoints_between > 1) ? (segment_dist / (kpoints_between - 1)) : 0.0;

        for (int b = 0; b < data.nbands; ++b) {
            for (int k = 0; k < kpoints_between; ++k) {
                double dist = cumulative_dist + k * step;
                const auto& kp = data.kpoints[s_idx + k];
                if (data.nspin == 2)
                    out << dist << " " << kp.energies_up[b] - e_fermi << " " << kp.energies_dn[b] - e_fermi << "\n";
                else
                    out << dist << " " << kp.energies_up[b] - e_fermi << "\n";
            }
            out << "\n";
        }

        cumulative_dist += segment_dist;
    }

    std::cout << "Success: Generated " << num_paths << " data files for plotting.\n";

    return 0;
}
