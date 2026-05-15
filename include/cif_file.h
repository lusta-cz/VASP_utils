#ifndef CIF_FILE_H_INCLUDED
#define CIF_FILE_H_INCLUDED

#include <string>

#include "poscar_file.h"

// Read a CIF file into a POSCAR structure.
// Symmetry operations from the CIF are applied to generate the full unit cell.
// Sites with occupancy below occ_threshold are discarded.
bool readCif(const std::string& filename, POSCAR& poscar, double occ_threshold = 0.5);

// Write a POSCAR structure as a P1 CIF file (all atoms listed explicitly).
bool writeCif(const std::string& filename, const POSCAR& poscar);

#endif  // CIF_FILE_H_INCLUDED
