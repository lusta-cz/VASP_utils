#ifndef POSCAR_SYMMETRY_H_INCLUDED
#define POSCAR_SYMMETRY_H_INCLUDED

#include <string>

bool readInput(int argc, char* argv[], std::string& inputFile, std::string& outputFile, double& symprec,
               bool& primitive, bool& wyckoff, bool& symoperation);
bool validateInput(const std::string& inputFile, const std::string& outputFile, double& symprec, bool& primitive);
void printHelp();

#endif  // POSCAR_SYMMETRY_H_INCLUDED
