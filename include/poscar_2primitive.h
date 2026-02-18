#ifndef POSCAR_2PRIMITIVE_H_INCLUDED
#define POSCAR_2PRIMITIVE_H_INCLUDED

#include <string>

bool readInput(int argc, char* argv[], std::string& inputFile, std::string& outputFile, double& symprec);
bool validateInput(const std::string& inputFile, const std::string& outputFile, double& symprec);
void printHelp();

#endif  // POSCAR_2PRIMITIVE_H_INCLUDED
