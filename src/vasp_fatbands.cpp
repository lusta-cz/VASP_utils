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
#include "poscar_file.h"

int main(int argc, char** argv) {
    CLI::App app{"VASP EIGENVAL to Band Structure Utility"};

    // Configuration and default values
    std::string eigen_file{"EIGENVAL"};
    std::string kpoints_file{"KPOINTS"};
    std::string doscar_file{"DOSCAR"};
    std::string outcar_file{"OUTCAR"};
    std::string procar_file{"PROCAR"};
    std::string poscar_file{"CONTCAR"};
    double manual_fermi{0.0};
    std::string fermiSource{"doscar"};
    std::string element{};

    // CLI Options
    app.add_option("-i,--input", eigen_file, "Input EIGENVAL file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_option("-k,--kpoints", kpoints_file, "KPOINTS file used for the band structure calculation")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--fermi-source", fermiSource, "Fermi level source: doscar (default) | outcar | manual | none")
        ->capture_default_str()
        ->check(CLI::IsMember({"doscar", "outcar", "manual", "none"}));
    app.add_option("--doscar", doscar_file, "DOSCAR file path (used with --fermi-source=doscar)")
        ->capture_default_str();
    app.add_option("--outcar", outcar_file, "OUTCAR file path (required even without --fermi-source=outcar)")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--procar", procar_file, "PROCAR file from the band structure calculation")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("--poscar", poscar_file, "POSCAR/CONTCAR file used for the band structure calculation")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
    app.add_option("-e,--element", element, "Target element label (e.g., Ni, Fe, O)")->required();
    auto* fermi_opt = app.add_option("-e,--fermi", manual_fermi, "Fermi level in eV (used with --fermi-source=manual)")
                          ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    double e_fermi{0.0};
    bool success{true};

    bool is_ncl = false;
    if (!checkNonCollinear(outcar_file, is_ncl)) {
        return 1;
    }

    if (is_ncl) {
        std::cerr << "Error: LNONCOLLINEAR tag was .TRUE.. That is NOT supported!\n";
        return 1;
    }

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

    POSCAR structure;
    if (!structure.readPOSCAR(poscar_file)) {
        return 1;
    }

    ProcarData procar_data;
    // --- Parsing from PROCAR file
    if (!parseProcar(procar_file, procar_data)) {
        return 1;
    }

    int global_atom_offset{0};
    int num_element_wanted{0};
    for (int i = 0; i < static_cast<int>(structure.elements.size()); i++) {
        if (element == structure.elements[i]) {
            num_element_wanted = structure.num_atoms[i];
            break;
        }
        global_atom_offset += structure.num_atoms[i];
    }

    if (num_element_wanted == 0) {
        std::cerr << "Error: Required element '" << element << "' not found in the structure file '" << poscar_file
                  << "'!\n";
        return 1;
    }

    std::cout << "Info: Identified " << num_element_wanted << " atom(s) matching element '" << element << "'.\n";

    // Deduce path segments
    int num_paths = data.total_kpoints / kpoints_between;

    for (int i_num_files = 0; i_num_files < num_element_wanted; i_num_files++) {
        // This calculates the exact global row index for this specific site inside the PROCAR data matrix
        int target_ion_idx = global_atom_offset + i_num_files;

        std::string out_name = "fatband_structure_" + element + "_" + std::to_string(i_num_files + 1) + ".dat";
        std::ofstream out(out_name);
        out << std::fixed << std::setprecision(10);

        for (int b = 0; b < data.nbands; ++b) {
            double cumulative_dist = 0.0;

            for (int p = 0; p < num_paths; ++p) {
                int s_idx = p * kpoints_between;
                int e_idx = s_idx + kpoints_between - 1;

                double dx = data.kpoints[e_idx].x - data.kpoints[s_idx].x;
                double dy = data.kpoints[e_idx].y - data.kpoints[s_idx].y;
                double dz = data.kpoints[e_idx].z - data.kpoints[s_idx].z;
                double segment_dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                double step = (kpoints_between > 1) ? (segment_dist / (kpoints_between - 1)) : 0.0;

                for (int k = 0; k < kpoints_between; ++k) {
                    double dist = cumulative_dist + k * step;
                    int global_k_idx = s_idx + k;
                    const auto& kp = data.kpoints[global_k_idx];

                    // Extract the orbital projection weight for this specific atom site.
                    // For example, we grab the 'total' projection column (which is the last index).
                    int total_column_idx = procar_data.orbital_labels.size() - 1;
                    double weight =
                        procar_data.kpoints[global_k_idx].bands[b].atom_weights[target_ion_idx][total_column_idx];

                    if (data.nspin == 2) {
                        // Include weight in output files for both spin channels
                        out << dist << " " << kp.energies_up[b] - e_fermi << " " << weight << " "
                            << kp.energies_dn[b] - e_fermi << " " << weight << "\n";
                    } else {
                        out << dist << " " << kp.energies_up[b] - e_fermi << " " << weight << "\n";
                    }
                }
                cumulative_dist += segment_dist;
            }
            out << "\n";
        }

        std::cout << "Success: Generated file '" << out_name << "' for element '" << element << "' number "
                  << (i_num_files + 1) << ".\n";
    }

    return 0;
}
