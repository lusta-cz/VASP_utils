#include "poscar_d2c.h"

#include <iostream>
#include <string>

#include "poscar_file.h"

bool readInput(int argc, char* argv[], std::string& inputFile, std::string& outputFile) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            printHelp();
            return false;
        } else if (arg == "--input") {
            if (i + 1 >= argc)
                return false;
            inputFile = argv[++i];
        } else if (arg == "--output") {
            if (i + 1 >= argc)
                return false;
            outputFile = argv[++i];
        } else {
            std::cerr << "Warning: unknown argument! Ignoring!\n";
            printHelp();
        }
    }
    return true;
}

void printHelp() {
    std::cerr << "Usage:\n"
                 "  poscar_d2c [options]\n\n"
                 "Options:\n"
                 "  --input   input POSCAR file name\n"
                 "  --output  output POSCAR file name\n"
                 "  --help    Show this help message\n\n"
                 "Example:\n"
                 "  poscar_d2c --input POSCARin --output POSCARout\n";
}

int main(int argc, char* argv[]) {
    std::string inputFile{"POSCAR"};
    std::string outputFile{"POSCAR_cartesian"};

    if (!readInput(argc, argv, inputFile, outputFile)) {
        if (argc > 1 && std::string(argv[1]) != "--help") {
            std::cerr << "Error parsing input arguments\n";
        }
        return 1;
    }

    POSCAR poscar;
    if (!poscar.readPOSCAR(inputFile)) {
        std::cerr << "Error reading POSCAR file: " << inputFile << "\n";
        return 1;
    }

    if (poscar.is_direct) {
        poscar.toCartesian();
        std::cout << "Converted Direct -> Cartesian\n";
    } else {
        std::cout << "POSCAR is already in Cartesian coordinates\n";
    }

    if (!poscar.writePOSCAR(outputFile)) {
        std::cerr << "Error writing output POSCAR file: " << outputFile << "\n";
        return 1;
    }

    std::cout << "Output written to: " << outputFile << "\n";

    return 0;
}
