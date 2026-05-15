#include <CLI/CLI.hpp>
#include <iostream>
#include <string>

#include "cif_file.h"
#include "poscar_file.h"

int main(int argc, char** argv) {
    CLI::App app{"Convert a VASP POSCAR file to CIF format (P1, all atoms explicit)"};

    std::string input_file{"POSCAR"};
    std::string output_file{"output.cif"};

    app.add_option("--input,-i", input_file, "Input POSCAR file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_option("--output,-o", output_file, "Output CIF file")->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    POSCAR poscar;
    if (!poscar.readPOSCAR(input_file))
        return 1;

    if (!writeCif(output_file, poscar))
        return 1;

    std::cout << "Written: " << output_file << "\n";
    return 0;
}
