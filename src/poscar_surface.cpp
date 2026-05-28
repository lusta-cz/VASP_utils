#include <spglib.h>

#include <CLI/CLI.hpp>
#include <iostream>
#include <string>

#include "poscar_file.h"
#include "surface.h"
#include "symmetry.h"

int main(int argc, char* argv[]) {
    CLI::App app{"Build a surface slab POSCAR from a bulk structure"};

    std::string inputFile{"POSCAR"};
    std::string outputFile{""};
    std::vector<int> miller(3);
    int n_layers{1};
    int n_frozen{1};
    double vacuum{15.0};
    bool use_primitive{false};
    double symprec{1e-5};
    bool center_slab{false};

    app.add_option("--input,-i", inputFile, "Input bulk POSCAR file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_option("--output,-o", outputFile, "Output POSCAR file (default: POSCAR_surface_hkl)");
    app.add_option("--miller,-m", miller, "Miller indices h k l (three integers)")->expected(3)->required();
    app.add_option("--layers,-n", n_layers, "Number of atomic layers")->required()->check(CLI::PositiveNumber);
    app.add_option("--frozen,-f", n_frozen, "Number of bottom layers frozen (Selective Dynamics F F F)")
        ->capture_default_str()
        ->check(CLI::NonNegativeNumber);
    app.add_option("--vacuum,-v", vacuum, "Vacuum thickness in Angstroms")
        ->capture_default_str()
        ->check(CLI::NonNegativeNumber);
    app.add_flag("--primitive,-p", use_primitive, "Convert input to primitive cell before building the slab");
    app.add_option("--symprec", symprec, "Symmetry tolerance for primitive-cell conversion")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    app.add_flag("--center,-c", center_slab, "Center the slab with equal vacuum layers on both sides");

    CLI11_PARSE(app, argc, argv);

    if (n_frozen > n_layers) {
        std::cerr << "Error: --frozen (" << n_frozen << ") must not exceed --layers (" << n_layers << ").\n";
        return 1;
    }

    if (outputFile.empty()) {
        outputFile =
            "POSCAR_surface_" + std::to_string(miller[0]) + std::to_string(miller[1]) + std::to_string(miller[2]);
    }

    POSCAR poscar;
    if (!poscar.readPOSCAR(inputFile)) {
        std::cerr << "Error: failed to read " << inputFile << "\n";
        return 1;
    }

    if (use_primitive) {
        if (symprec > 1e-3)
            std::cerr << "Warning: symprec " << symprec << " is high — consider using default 1e-5.\n";

        auto prim = standardizeCell(poscar, symprec, /*primitive=*/1, /*idealize=*/0);
        if (!prim) {
            std::cerr << "Error: failed to convert to primitive cell.\n";
            return 1;
        }
        poscar = *prim;
        std::cout << "Info: using primitive cell (" << poscar.total_atoms << " atoms).\n";
    }

    SlabOptions opts;
    opts.h = miller[0];
    opts.k = miller[1];
    opts.l = miller[2];
    opts.n_layers = n_layers;
    opts.n_frozen = n_frozen;
    opts.vacuum = vacuum;
    opts.symprec = symprec;
    opts.center = center_slab;

    POSCAR slab = buildSlab(poscar, opts);
    if (slab.total_atoms == 0) {
        std::cerr << "Error: slab construction failed.\n";
        return 1;
    }

    std::cout << "Slab: " << slab.total_atoms << " atoms, " << "vacuum " << vacuum << " Å, "
              << "frozen layers: " << n_frozen << "\n";

    if (!slab.writePOSCAR(outputFile)) {
        std::cerr << "Error: failed to write " << outputFile << "\n";
        return 1;
    }

    std::cout << "Written to " << outputFile << "\n";
    return 0;
}
