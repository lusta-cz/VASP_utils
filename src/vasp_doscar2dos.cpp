#include <CLI/CLI.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "doscar.h"

namespace fs = std::filesystem;

static void writeHeader(std::ofstream& out, const DoscarData& data, const std::string& title, bool shifted) {
    out << "# " << title << "\n"
        << "# NIONS = " << data.nions << "   NKPTS = " << data.nkpts << "   NEDOS = " << data.nedos
        << "   nspin = " << data.nspin << "\n";

    out << std::fixed << std::setprecision(6);
    if (shifted) {
        out << "# E_Fermi = " << data.e_fermi << " eV (shifted to 0)\n"
            << "# E_range = [" << (data.emin - data.e_fermi) << ", " << (data.emax - data.e_fermi) << "] eV\n";
    } else {
        out << "# E_Fermi = " << data.e_fermi << " eV (NOT shifted — energies are absolute)\n"
            << "# E_range = [" << data.emin << ", " << data.emax << "] eV\n";
    }
}

static void writeColNames(std::ofstream& out, const std::vector<std::string>& cols) {
    out << "#";
    for (const auto& c : cols)
        out << "  " << c;
    out << "\n";
}

static bool writeRows(std::ofstream& out, const std::vector<std::vector<double>>& rows, double energy_shift) {
    out << std::fixed << std::setprecision(8);
    for (const auto& row : rows) {
        out << (row[0] - energy_shift);
        for (std::size_t j = 1; j < row.size(); ++j)
            out << "  " << row[j];
        out << "\n";
    }
    return out.good();
}

int main(int argc, char** argv) {
    CLI::App app{
        "VASP DOSCAR to DOS data files\n"
        "Writes tdos.dat and (if partial DOS present) pdos_atom_N.dat.\n"
        "Energy axis is shifted so E_Fermi = 0 by default."};

    std::string doscar_file{"DOSCAR"};
    std::string output_dir{"."};
    std::string prefix;
    bool no_shift{false};
    bool total_only{false};

    app.add_option("--input,-i", doscar_file, "Input DOSCAR file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_flag("--no-shift", no_shift, "Disable E_Fermi shift (energies written as absolute values)");
    app.add_option("--output-dir,-d", output_dir, "Directory for output files")->capture_default_str();
    app.add_option("--prefix,-p", prefix, "Filename prefix, e.g. 'run1_' produces run1_tdos.dat");
    app.add_flag("--total-only", total_only, "Write only tdos.dat, skip partial DOS files");

    CLI11_PARSE(app, argc, argv);

    // --- Parse DOSCAR ---
    DoscarData data;
    if (!parseDoscar(doscar_file, data))
        return 1;

    const bool shifted = !no_shift;
    const double energy_shift = shifted ? data.e_fermi : 0.0;

    if (shifted)
        std::cout << "E_Fermi from " << doscar_file << ": " << data.e_fermi << " eV\n";
    else
        std::cout << "E_Fermi shift disabled. E_Fermi = " << data.e_fermi << " eV (from " << doscar_file
                  << ", not applied).\n";

    // --- Validate output directory ---
    fs::path outdir{output_dir};
    if (!fs::is_directory(outdir)) {
        std::cerr << "Error: output directory does not exist: " << output_dir << "\n";
        return 1;
    }

    auto outpath = [&](const std::string& name) { return (outdir / (prefix + name)).string(); };

    int files_written = 0;

    // --- Total DOS ---
    {
        auto path = outpath("tdos.dat");
        std::ofstream out(path);
        if (!out) {
            std::cerr << "Error: cannot write " << path << "\n";
            return 1;
        }
        writeHeader(out, data, "VASP DOSCAR: Total Density of States", shifted);
        writeColNames(out, data.tdos_col_names);
        writeRows(out, data.tdos, energy_shift);
        std::cout << "Written: " << path << "\n";
        ++files_written;
    }

    // --- Partial DOS ---
    if (!total_only && data.has_partial) {
        for (int ia = 0; ia < static_cast<int>(data.pdos.size()); ++ia) {
            auto path = outpath("pdos_atom_" + std::to_string(ia + 1) + ".dat");
            std::ofstream out(path);
            if (!out) {
                std::cerr << "Error: cannot write " << path << "\n";
                return 1;
            }
            writeHeader(out, data, "VASP DOSCAR: Partial Density of States — Atom " + std::to_string(ia + 1), shifted);
            writeColNames(out, data.pdos_col_names);
            writeRows(out, data.pdos[ia], energy_shift);
            std::cout << "Written: " << path << "\n";
            ++files_written;
        }
    } else if (!total_only && !data.has_partial) {
        std::cout << "Info: no partial DOS in " << doscar_file << "\n";
    }

    std::cout << "Done: " << files_written << " file(s) written.\n";
    return 0;
}
