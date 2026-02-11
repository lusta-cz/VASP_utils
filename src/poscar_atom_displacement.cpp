#include "poscar_atom_displacement.h"
#include "poscar_file.h"


#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>


bool testReadPOSCAR(const std::string& filename)
{
    POSCAR original;

    if (!original.readPOSCAR(filename)) {
        std::cerr << "Error: reading POSCAR file from " << filename << "\n";
        return false;
    }
    if (!original.writePOSCAR("POSCAR_out_1")) {
        std::cerr << "Error: reading POSCAR file from " << filename << "\n";
        return false;
    }

    /*// --- Print the POSCAR contents ---
    std::cout << "Comment: " << original.comment << "\n";
    std::cout << "Scale: " << original.scale << "\n";

    std::cout << "Lattice vectors:\n";
    for (int i = 0; i < 3; ++i)
        std::cout << std::fixed << std::setprecision(9)
                  << original.lattice[i][0] << " "
                  << original.lattice[i][1] << " "
                  << original.lattice[i][2] << "\n";

    std::cout << "Elements: ";
    for (auto &e : original.elements) std::cout << e << " ";
    std::cout << "\n";

    std::cout << "Number of atoms: ";
    for (auto &n : original.num_atoms) std::cout << n << " ";
    std::cout << "\n";

    std::cout << "Selective dynamics: "
              << (original.selective_dynamics ? "yes" : "no") << "\n";
    std::cout << "Coordinate type: "
              << (original.is_direct ? "Direct" : "Cartesian") << "\n";

    std::cout << "Coordinates:\n";
    for (size_t i = 0; i < original.coordinates.size(); ++i)
        std::cout << std::fixed << std::setprecision(9)
                  << original.coordinates[i].x << " "
                  << original.coordinates[i].y << " "
                  << original.coordinates[i].z << "\n";*/
    return true;
}

bool readInput(int argc,
               char* argv[],
               std::string& filename,
               int& n_files,
               int& n_atoms,
               double& amplitude, bool& allAtoms)
{
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            printHelp();
            return false;
        }
        else if (arg == "--nfiles") {
            if (i + 1 >= argc) return false;
            try {
                n_files = std::stoi(argv[++i]);
            }
            catch (...) {
                return false;
            }
        }
        else if (arg == "--input") {
            if (i + 1 >= argc) return false;
            filename = argv[++i];
        }
        else if (arg == "--allatoms") {
            allAtoms = true;
        }
        else if (arg == "--natoms") {
            if (i + 1 >= argc) return false;
            try {
                n_atoms = std::stoi(argv[++i]);
            }
            catch (...) {
                return false;
            }
        }
        else if (arg == "--amp") {
            if (i + 1 >= argc) return false;
            try {
                amplitude = std::stod(argv[++i]);
            }
            catch (...) {
                return false;
            }
        }
        else {
            std::cerr << "Warning: unknown argument! Ignoring it!\n";
            printHelp();
        }
    }
    return true;
}


bool validateInput(const std::string& filename,
                   int n_files,
                   int n_atoms,
                   double amplitude)
{
    //Testing input for valid values.
    if(n_atoms < 0)
    {
        std::cerr << "Error: number of atoms to displace is negative!\n";
        return false;
    }

    if(n_atoms == 0)
    {
        std::cerr << "Error: number of atoms to displace is 0!\n";
        return false;
    }

    if(n_files < 0)
    {
        std::cerr << "Error: number displaced structure files is negative!\n";
        return false;
    }

    if(n_files == 0)
    {
        std::cerr << "Error: displaced structure files is 0!\n";
        return false;
    }

    if(n_files > 1000)
    {
        std::cerr << "Error: number displaced structure files is too high!\n";
        return false;
    }

    if(amplitude < 0)
    {
        std::cerr << "Error: amplitude of displacement is negative!\n";
        return false;
    }
    if(amplitude > 0.1)
    {
        std::cerr << "Warning: amplitude is above 10%!!! That can cause errors running VASP!\n";
    }

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }
    return true;
}

void printHelp()
{
    std::cerr <<
        "Usage:\n"
        "  poscar_atom_displace [options]\n\n"
        "Options:\n"
        "  --input      POSCAR file name\n"
        "  --nfiles     number of displaced structure files to create\n"
        "  --natoms     number of atoms to displace\n"
        "  --allatoms   displace all atoms in the input file\n"
        "  --amp        maximal norm of the displacement vector in Angstroms\n"
        "  --help       Show this help message\n\n"
        "Example:\n"
        "  poscar_atom_displace --input POSCAR --nfiles 10 --natoms 1 --amp 0.1\n";
}

int main(int argc, char* argv[])
{
//Iniciaization of Input file and parameters.
    std::string filename {"POSCAR"};
    int n_atoms {1};
    int n_files {1};
    double amplitude {0.01};
    bool allAtoms {false};

//Read input arguments
    if (!readInput(argc, argv, filename, n_files, n_atoms, amplitude, allAtoms))
        return 1;




    if (!validateInput(filename, n_files, n_atoms, amplitude))
        return 1;



    //if (!testReadPOSCAR(filename)) return 1;

    POSCAR original;

    if (!original.readPOSCAR(filename)) {
        std::cerr << "Error: reading POSCAR file " << filename << "\n";
        return 1;
    }

    if(allAtoms){
        n_atoms = original.total_atoms;
        std::cout << "Displacing all atoms in input file.\n";
        if (allAtoms && n_atoms != 1) {
            std::cerr << "Note: --allatoms overrides --natoms\n";
        }
    }

    else if(n_atoms > original.total_atoms)
    {
        std::cerr << "Warning: number of atoms in POSCAR file is smaller than number of atoms to be displaced that you required!\n"
                    "Displacing all avaiable atoms!\n";
        n_atoms = original.total_atoms;
        std::cout << "Number of atoms to displace: " << n_atoms << "\n";
    }

    for(int j = 0; j < n_files; j++)
    {
    POSCAR output(original);
    std::string filenameOut = "POSCAR_modified" + std::to_string(j + 1);

    output.displaceAtoms(n_atoms, amplitude);
    output.writePOSCAR(filenameOut);
    }


    /*
    std::cout << "Input file: " << filename << "\n";
    std::cout << "Number of atoms to displace: " << n_atoms << "\n";
    std::cout << "Displacement amplitude: " << amplitude << "\n";
    */


    return 0;
}
