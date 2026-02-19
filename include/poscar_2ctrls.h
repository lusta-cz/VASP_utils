#ifndef POSCAR_2CTRLS_INCLUDED
#define POSCAR_2CTRLS_INCLUDED

#include <string>

bool readInput(int argc, char* argv[], std::string& inputFile, std::string& outputFile);
void printHelp();

#endif  // POSCAR_2CTRLS_INCLUDED
