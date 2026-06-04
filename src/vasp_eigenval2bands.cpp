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
    std::string fermiSource{"outcar"};

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
    app.add_option("--outcar", outcar_file, "OUTCAR file path (required even without --fermi-source=outcar)")
        ->capture_default_str()
        ->check(CLI::ExistingFile);
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

    BZPath bz_path;
    // --- Parsing from KPOINTS file
    if (!parseKpoints(kpoints_file, bz_path)) {
        return 1;
    }

    EigenvalData data;
    // --- Parsing from EIGENVAL file
    if (!parseEigenval(eigen_file, data)) {
        return 1;
    }

    int kpoints_between = bz_path.kpts_per_seg;
    // --- Path Calculation and Output Generation ---
    if (data.total_kpoints == 0 || kpoints_between == 0) {
        std::cerr << "Error: Something went wrong while reading the number of k points!\n";
        return 1;
    }

    int num_paths = data.total_kpoints / kpoints_between;
    // Integrity Verification
    if (num_paths != bz_path.num_segments) {
        std::cerr << "Warning: EIGENVAL tracks count (" << num_paths << ") does not match KPOINTS structural count ("
                  << bz_path.num_segments << ")!\n";
    }

    // =========================================================================
    // CALL SHARED REFACTORED LOG FUNCTION
    // =========================================================================
    double log_cumulative_dist = 0.0;
    if (!writeKpointsLog("kpoints.log", bz_path, data, log_cumulative_dist)) {
        return 1;
    }
    std::cout << "Success: Generated automated high-symmetry coordinate mapping -> 'kpoints.log'.\n";

    // =========================================================================
    // GENERATING INDEPENDENT SPLIT PATH FILES (With Discontinuity tracking)
    // =========================================================================
    double split_cumulative_dist = 0.0;
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

        // ─── HERE IS THE CONTINUITY CHECK FOR SEPARATE FILES ───
        if (p > 0) {
            std::string prev_end = bz_path.segments[p - 1].end_label;
            std::string curr_start = bz_path.segments[p].start_label;
            if (prev_end != curr_start) {
                split_cumulative_dist += 0.25;  // Shifts the starting coordinate of the new file forward!
            }
        }

        for (int b = 0; b < data.nbands; ++b) {
            for (int k = 0; k < kpoints_between; ++k) {
                double dist = split_cumulative_dist + k * step;
                const auto& kp = data.kpoints[s_idx + k];
                if (data.nspin == 2)
                    out << dist << " " << kp.energies_up[b] - e_fermi << " " << kp.energies_dn[b] - e_fermi << "\n";
                else
                    out << dist << " " << kp.energies_up[b] - e_fermi << "\n";
            }
            out << "\n";
        }
        split_cumulative_dist += segment_dist;
    }

    // =========================================================================
    // GENERATING COMBINED TRACE PATH FILE (bands_all.dat with continuity check)
    // =========================================================================
    std::string combined_name = "bands_all.dat";
    std::ofstream comb_out(combined_name);
    comb_out << std::fixed << std::setprecision(10);
    comb_out << "# X_Distance   Energy_Up   [Energy_Dn]\n";

    double comb_cumulative_dist = 0.0;

    for (int b = 0; b < data.nbands; ++b) {
        // Reset distance tracker for each discrete band layout sequence
        comb_cumulative_dist = 0.0;

        for (int p = 0; p < num_paths; ++p) {
            int s_idx = p * kpoints_between;
            int e_idx = s_idx + kpoints_between - 1;

            double dx = data.kpoints[e_idx].x - data.kpoints[s_idx].x;
            double dy = data.kpoints[e_idx].y - data.kpoints[s_idx].y;
            double dz = data.kpoints[e_idx].z - data.kpoints[s_idx].z;
            double segment_dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            double step = (kpoints_between > 1) ? (segment_dist / (kpoints_between - 1)) : 0.0;

            if (p > 0) {
                std::string prev_end = bz_path.segments[p - 1].end_label;
                std::string curr_start = bz_path.segments[p].start_label;
                if (prev_end != curr_start) {
                    comb_cumulative_dist += 0.25;
                    // TRACE BREAK: Double blank rows inform Gnuplot to drop drawing connectivity lines
                    // across the empty gap length
                    comb_out << "\n\n";
                }
            }

            for (int k = 0; k < kpoints_between; ++k) {
                double dist = comb_cumulative_dist + k * step;
                const auto& kp = data.kpoints[s_idx + k];

                if (data.nspin == 2) {
                    comb_out << dist << " " << kp.energies_up[b] - e_fermi << " " << kp.energies_dn[b] - e_fermi
                             << "\n";
                } else {
                    comb_out << dist << " " << kp.energies_up[b] - e_fermi << "\n";
                }
            }
            comb_cumulative_dist += segment_dist;
        }
        comb_out << "\n\n";  // Isolate complete band sequence trails cleanly from each other
    }

    // =========================================================================
    // SANITY CHECK
    // =========================================================================
    if (std::abs(split_cumulative_dist - log_cumulative_dist) > 1e-5 ||
        std::abs(comb_cumulative_dist - log_cumulative_dist) > 1e-5) {
        std::cerr << "Warning: Something went wrong and the 'kpoints.log' (" << log_cumulative_dist
                  << "), split files (" << split_cumulative_dist << "), or 'bands_all.dat' (" << comb_cumulative_dist
                  << ") do not have matching final x-coordinates!\n";
    }

    std::cout << "Success: Generated " << num_paths << " data files for plotting.\n";
    return 0;
}
