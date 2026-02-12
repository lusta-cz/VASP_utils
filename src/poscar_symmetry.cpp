#include "poscar_symmetry.h"

#include <spglib.h>

#include <fstream>
#include <iostream>
#include <optional>
#include <string>

#include "poscar_file.h"
#include "symmetry.h"

bool readInput(int argc, char* argv[], std::string& inputFile, std::string& outputFile, double& symprec,
               bool& primitive, bool& wyckoff, bool& symoperation) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            printHelp();
            return false;
        } else if (arg == "--input") {
            if (i + 1 >= argc)
                return false;
            inputFile = argv[++i];
        } else if (arg == "--primitive") {
            primitive = true;
        } else if (arg == "--wyckoff") {
            wyckoff = true;
        } else if (arg == "--symoper") {
            symoperation = true;
        } else if (arg == "--output") {
            if (i + 1 >= argc)
                return false;
            outputFile = argv[++i];
        } else if (arg == "--symprec") {
            if (i + 1 >= argc)
                return false;
            try {
                symprec = std::stod(argv[++i]);
            } catch (...) {
                return false;
            }
        } else {
            std::cerr << "Warning: unknown argument! Ignoring it!\n";
            printHelp();
        }
    }
    return true;
}

bool validateInput(const std::string& inputFile, const std::string& outputFile, double& symprec, bool& primitive) {
    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Error: cannot open file " << inputFile << "\n";
        return false;
    }
    std::ifstream file_2(outputFile);
    if (file_2 && primitive) {
        std::cerr << "Warning: output file " << outputFile << " already exists!!! Will be overwritten!!!\n";
    }
    if (symprec < 0) {
        std::cerr << "Error: symprec cannot be negative!!!\n";
        return false;
    }
    if (symprec == 0) {
        std::cerr << "Error: symprec cannot be 0!!!\n";
        return false;
    }
    if (symprec > 1e-3) {
        std::cerr << "Warning: symprec is too high! Consider using default value 1e-5.\n";
    }

    return true;
}

void printHelp() {
    std::cerr << "Usage:\n"
                 "  poscar_symmetry [options]\n\n"
                 "Options:\n"
                 "  --input      input POSCAR file name (default: POSCAR)\n"
                 "  --symprec    symmetry tolerance (spglib symprec) (default: 1e-5)\n"
                 "  --wyckoff    print the Wyckoff positions\n"
                 "  --symoper    print the symmetry operations\n"
                 "  --primitive  generate primitive POSCAR file if the input cell is not primitive\n"
                 "  --output     output POSCAR file name (used with --primitive) (default: POSCAR_primitive)\n"
                 "  --help       show this help message\n\n"
                 "Example:\n"
                 "  poscar_symmetry --input POSCARin --symprec 1e-5 --primitive\n";
}

int main(int argc, char* argv[]) {
    std::string inputFile{"POSCAR"};
    std::string outputFile{"POSCAR_primitive"};
    bool primitive{false};
    double symprec{1e-5};
    bool wyckoff{false};
    bool symoperation{false};

    if (!readInput(argc, argv, inputFile, outputFile, symprec, primitive, wyckoff, symoperation)) {
        return 1;
    }

    if (!validateInput(inputFile, outputFile, symprec, primitive)) {
        return 1;
    }

    POSCAR poscar;
    poscar.readPOSCAR(inputFile);

    auto dataset = analyzeSymmetry(poscar, symprec);

    if (dataset) {
        printSymmetryInfo(*dataset, wyckoff, symoperation);
    } else {
        std::cerr << "Error: failed to analyze symmetry.\n";
        return 1;
    }

    if (primitive) {
        std::cout << "Creating the primitive cell file.\n";
        auto maybePrimitive = makePrimitiveCell(poscar, symprec);

        if (!maybePrimitive) {
            std::cerr << "Error: failed to create primitive POSCAR.\n";
            return 1;
        }

        if (!maybePrimitive->writePOSCAR(outputFile)) {
            std::cerr << "Error: failed to write primitive POSCAR file to " << outputFile << "\n";
            return 1;
        }
    }

    return 0;
}
