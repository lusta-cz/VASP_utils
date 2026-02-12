#include "poscar_file.h"

#include "random_utility.h"

// Linear algebra

#include <cblas.h>
#include <lapacke.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

bool POSCAR::skipLines(std::ifstream& file, int n) {
    std::string tmp;
    for (int i = 0; i < n; ++i) {
        if (!std::getline(file, tmp))
            return false;
    }
    return true;
}

bool POSCAR::readPOSCARHeader(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }

    std::string line;

    // Reading Line 1: comment
    if (!std::getline(file, comment))
        return false;

    // Reading Line 2: scale
    if (!std::getline(file, line))
        return false;

    scale = std::stod(line);

    // Reading Lines 3-5: lattice vectors
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(file, line))
            return false;

        std::istringstream iss(line);
        iss >> lattice[i][0] >> lattice[i][1] >> lattice[i][2];
    }

    // Reading Line 6: element symbols
    if (!std::getline(file, line))
        return false;

    std::istringstream iss_elements(line);
    std::string elem;
    elements.clear();
    while (iss_elements >> elem) {
        elements.push_back(elem);
    }

    // Reading Line 7: number of atoms per element
    if (!std::getline(file, line))
        return false;

    std::istringstream iss_counts(line);
    num_atoms.clear();
    int n;
    while (iss_counts >> n) {
        num_atoms.push_back(n);
    }

    return true;
}

bool POSCAR::readPOSCAROptional(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }

    std::string line;

    // Skip first 7 header lines (already read)
    if (!skipLines(file, 7))
        return false;

    // Reading Optional: Selective dynamics or Direct/Cartesian
    if (!std::getline(file, line))
        return false;
    if (line[0] == 'S' || line[0] == 's') {
        selective_dynamics = true;

        // For now selective dynamics is NOT supported!!! Remove if implemented with this keyword!!!
        std::cerr << "Error: selective dynamics is NOT supported yet!\n";
        return false;

        // Reading next line (Direct/Cartesian) if Selective dynamic is present
        if (!std::getline(file, line))
            return false;
    } else {
        selective_dynamics = false;
    }

    if (line[0] == 'D' || line[0] == 'd') {
        is_direct = true;
    } else
        is_direct = false;

    return true;
}

bool POSCAR::readPOSCARCoordinates(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }

    std::string line;

    // Skip first 7 header lines (already read)
    if (!skipLines(file, 7))
        return false;

    // Skip optional line for Selective dynamics if present
    if (selective_dynamics) {
        if (!std::getline(file, line))
            return false;
    }

    // Skip Direct/Cartesian line (already validated)
    if (!std::getline(file, line))
        return false;

    // Now read all coordinates
    for (size_t i = 0; i < coordinates.size(); ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: not enough coordinate lines in POSCAR\n";
            return false;
        }

        std::istringstream iss(line);
        if (!(iss >> coordinates[i].x >> coordinates[i].y >> coordinates[i].z)) {
            std::cerr << "Error: failed to parse coordinates for atom " << i << "\n";
            return false;
        }
    }

    return true;
}

bool POSCAR::readPOSCAR(const std::string& filename) {
    if (!readPOSCARHeader(filename)) {
        std::cerr << "Error: reading POSCAR header from " << filename << "\n";
        return false;
    }

    // Compute total number of atoms
    for (int count : num_atoms)
        total_atoms += count;

    coordinates.resize(total_atoms);

    if (!readPOSCAROptional(filename)) {
        std::cerr << "Error reading POSCAR structural keywords from " << filename << "\n";
        return false;
    }

    if (!readPOSCARCoordinates(filename)) {
        std::cerr << "Error reading POSCAR coordinates from " << filename << "\n";
        return false;
    }

    return true;
}

bool POSCAR::writePOSCAR(const std::string& filenameOut) {
    std::ifstream fileTest(filenameOut);
    if (fileTest.good()) {
        std::cerr << "Warning: file \"" << filenameOut << "\" already exists and will be overwritten.\n";
    }

    std::ofstream file(filenameOut);
    if (!file) {
        std::cerr << "Error: cannot create file " << filenameOut << "\n";
        return false;
    }

    // Writing Line 1: comment
    file << comment << "\n";

    // Writing Line 2: scale
    file << std::fixed << std::setprecision(10) << scale << "\n";

    // Writing Lines 3-5: lattice vectors
    for (int i = 0; i < 3; ++i)
        file << std::fixed << std::setprecision(10) << lattice[i][0] << " " << lattice[i][1] << " " << lattice[i][2]
             << "\n";

    // Writing Line 6: element symbols
    for (size_t i = 0; i < elements.size(); ++i)
        file << elements[i] << " ";
    file << "\n";

    // Writing Line 7: number of atoms
    for (size_t i = 0; i < num_atoms.size(); ++i)
        file << num_atoms[i] << " ";
    file << "\n";

    // Writing Optional: selective dynamics
    if (selective_dynamics)
        file << "Selective Dynamics\n";

    // Writing Direct/Cartesian
    file << (is_direct ? "Direct" : "Cartesian") << "\n";

    // Writing Atomic coordinates
    for (size_t i = 0; i < coordinates.size(); ++i)
        file << std::fixed << std::setprecision(10) << coordinates[i].x << " " << coordinates[i].y << " "
             << coordinates[i].z << "\n";

    // TO DO for Selective dynamics: include selective dynamics flags if needed

    file.flush();
    if (!file) {
        std::cerr << "Error: failed writing to " << filenameOut << "\n";
        return false;
    }

    file.close();

    return true;
}

void POSCAR::displaceAtom(size_t atom_index, double amplitude) {
    if (atom_index >= coordinates.size())
        return;

    // Generate random vector and normalize it
    double ex = randomDouble(-1.0, 1.0);
    double ey = randomDouble(-1.0, 1.0);
    double ez = randomDouble(-1.0, 1.0);

    double norm = std::sqrt(ex * ex + ey * ey + ez * ez);

    if (norm < 1e-12)
        return;

    ex = (ex / norm);
    ey = (ey / norm);
    ez = (ez / norm);

    // Generate random norm of the random vector, cubicroot needed for uniform representation of volume due to expanding
    // sphere with r^3

    double r = amplitude * std::cbrt(randomDouble(0.0, 1.0));

    // Apply to the atom
    coordinates[atom_index].x += ex * r;
    coordinates[atom_index].y += ey * r;
    coordinates[atom_index].z += ez * r;
}

void POSCAR::displaceAtoms(int n_atoms, double amplitude) {
    // Create vector with numbers from 0 to total_atoms-1 and then doing random permutation
    std::vector<size_t> indices(total_atoms);
    for (int i = 0; i < total_atoms; ++i)
        indices[i] = i;

    std::shuffle(indices.begin(), indices.end(), getGenerator());

    // Convert to Cartesian if needed
    bool was_direct = is_direct;
    if (is_direct)
        toCartesian();

    // Displace selected atoms
    for (int i = 0; i < n_atoms; ++i) {
        displaceAtom(indices[i], amplitude);  // already Cartesian
    }

    // Convert back to Direct if needed
    if (was_direct)
        toDirect();
}

/*
OLD version
void POSCAR::displaceAtoms(int n_atoms, AmpMode amp_mode, double amplitude)
{
    // Generate random atom indices to displace
    std::vector<size_t> indices(total_atoms);
    for (int i = 0; i < total_atoms; ++i) indices[i] = i;

    // Shuffle the indices randomly
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.end(), gen);

    // Take the first n_atoms indices and displace those atoms
    for (int i = 0; i < n_atoms; ++i){
        if (amp_mode == AmpMode::Direct) {
            displaceAtomDirect(indices[i], amplitude);
        }
        else if (amp_mode == AmpMode::Cartesian) {
            displaceAtomCartesian(indices[i], amplitude);
        }
    }

}*/

void POSCAR::toDirect() {
    if (is_direct)
        return;

    double A[9];

    // Copy lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A[i * 3 + j] = lattice[i][j];

    lapack_int ipiv[3];

    // LU factorization
    if (LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, A, 3, ipiv) != 0) {
        std::cerr << "Error: LU decomposition failed.\n";
        return;
    }

    // Compute inverse
    if (LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, A, 3, ipiv) != 0) {
        std::cerr << "Error: Matrix inversion failed.\n";
        return;
    }

    // Transform coordinates
    for (auto& atom : coordinates) {
        double x[3] = {atom.x, atom.y, atom.z};
        double y[3];

        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 0.0, y, 1);

        atom.x = y[0];
        atom.y = y[1];
        atom.z = y[2];
    }

    is_direct = true;
}

void POSCAR::toCartesian() {
    if (!is_direct)
        return;

    double A[9];

    // Flatten lattice (row-major)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A[i * 3 + j] = lattice[i][j];

    for (auto& atom : coordinates) {
        double x[3] = {atom.x, atom.y, atom.z};
        double y[3];

        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 0.0, y, 1);

        atom.x = y[0];
        atom.y = y[1];
        atom.z = y[2];
    }

    is_direct = false;
}
