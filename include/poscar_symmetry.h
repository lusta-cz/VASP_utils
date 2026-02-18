#ifndef POSCAR_SYMMETRY_H_INCLUDED
#define POSCAR_SYMMETRY_H_INCLUDED

#include <string>

bool readInput(int argc, char* argv[], std::string& inputFile, double& symprec, bool& wyckoff, bool& symoperation);
bool validateInput(const std::string& inputFile, double& symprec);
void printHelp();

#endif  // POSCAR_SYMMETRY_H_INCLUDED
