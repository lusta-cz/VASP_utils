#include "random_utility.h"

// Single RNG engine for the whole program
static std::random_device rd;
static std::mt19937 gen(rd());

double randomDouble(double min, double max)
{
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

void seedRandom(unsigned int seed)
{
    gen.seed(seed);
}

std::mt19937& getGenerator()
{
    return gen;
}
