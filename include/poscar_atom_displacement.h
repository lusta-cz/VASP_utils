#ifndef POSCAR_ATOM_DISPLACEMENT_H_INCLUDED
#define POSCAR_ATOM_DISPLACEMENT_H_INCLUDED

#include <string>

bool readInput(int argc, char* argv[], std::string& filename, int& n_files, int& n_atoms, double& amplitude);

bool validateInput(const std::string& filename, int n_files, int n_atoms, double amplitude);

void printHelp();

#endif  // POSCAR_ATOM_DISPLACEMENT_H_INCLUDED
