#ifndef POSCAR_FILE_H_INCLUDED
#define POSCAR_FILE_H_INCLUDED


#include <string>
#include <vector>

struct Atom {
    double x, y, z;  // fractional or Cartesian
};

struct POSCAR {
    std::string comment {"System"};          // first line
    double scale {1.000000};                 // scaling factor
    double lattice[3][3];         // 3x3 lattice vectors
    std::vector<std::string> elements; // element symbols
    std::vector<int> num_atoms;        // number of atoms per element
    bool selective_dynamics = false;   // optional
    bool is_direct = true;             // true = Direct, false = Cartesian
    std::vector<Atom> coordinates;     // Nx3 coordinates
    int total_atoms {0};


    bool readPOSCAR(const std::string& filename);
    bool writePOSCAR(const std::string& filenameOut);
    void displaceAtoms(int n_atoms, double amplitude);
    void toDirect();
    void toCartesian();

    private:
    bool readPOSCARHeader(const std::string& filename);
    bool readPOSCAROptional(const std::string& filename);
    bool skipLines(std::ifstream& file, int n);
    bool readPOSCARCoordinates(const std::string& filename);
    void displaceAtom(size_t atom_index, double amplitude);

};
#endif // POSCAR_FILE_H_INCLUDED
