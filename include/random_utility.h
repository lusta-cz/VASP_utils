#ifndef RANDOM_UTILITY_H_INCLUDED
#define RANDOM_UTILITY_H_INCLUDED

#include <random>

double randomDouble(double min, double max);
void seedRandom(unsigned int seed);
extern std::mt19937& getGenerator();


#endif // RANDOM_UTILITY_H_INCLUDED
