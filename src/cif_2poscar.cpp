#include <CLI/CLI.hpp>
#include <iostream>
#include <string>

#include "cif_file.h"
#include "poscar_file.h"

int main(int argc, char** argv) {
    CLI::App app{
        "Convert a CIF file to VASP POSCAR format\n"
        "Symmetry operations in the CIF are applied to produce the full unit cell."};

    std::string input_file{"input.cif"};
    std::string output_file{"POSCAR"};
    double occ_threshold{0.5};

    app.add_option("--input,-i", input_file, "Input CIF file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_option("--output,-o", output_file, "Output POSCAR file")->capture_default_str();
    app.add_option("--occ-threshold", occ_threshold, "Minimum occupancy to include a site (default: 0.5)")
        ->capture_default_str()
        ->check(CLI::Range(0.0, 1.0));

    CLI11_PARSE(app, argc, argv);

    POSCAR poscar;
    if (!readCif(input_file, poscar, occ_threshold))
        return 1;

    if (!poscar.writePOSCAR(output_file))
        return 1;

    std::cout << "Written: " << output_file << "  (" << poscar.total_atoms << " atoms)\n";
    return 0;
}
