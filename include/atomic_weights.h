#ifndef ATOMIC_WEIGHTS_H_INCLUDED
#define ATOMIC_WEIGHTS_H_INCLUDED

#include <string>
#include <unordered_map>

/**
 * @brief Retrieves the standard atomic weight (g/mol) of a chemical element.
 * @param element The chemical symbol (e.g., "Fe", "O").
 * @return Molar mass in g/mol, or 0.0 if not found.
 */
inline double getAtomicWeight(const std::string& element) {
    static const std::unordered_map<std::string, double> weights = {
        {"H", 1.008},   {"He", 4.0026}, {"Li", 6.94},   {"Be", 9.0122}, {"B", 10.81},   {"C", 12.011},  {"N", 14.007},
        {"O", 15.999},  {"F", 18.998},  {"Ne", 20.180}, {"Na", 22.990}, {"Mg", 24.305}, {"Al", 26.982}, {"Si", 28.085},
        {"P", 30.974},  {"S", 32.06},   {"Cl", 35.45},  {"Ar", 39.948}, {"K", 39.098},  {"Ca", 40.078}, {"Sc", 44.956},
        {"Ti", 47.867}, {"V", 50.942},  {"Cr", 51.996}, {"Mn", 54.938}, {"Fe", 55.845}, {"Co", 58.933}, {"Ni", 58.693},
        {"Cu", 63.546}, {"Zn", 65.38},  {"Ga", 69.723}, {"Ge", 72.63},  {"As", 74.922}, {"Se", 78.971}, {"Br", 79.904},
        {"Kr", 83.798}, {"Rb", 85.468}, {"Sr", 87.62},  {"Y", 88.906},  {"Zr", 91.224}, {"Nb", 92.906}, {"Mo", 95.95},
        {"Tc", 98.0},   {"Ru", 101.07}, {"Rh", 102.91}, {"Pd", 106.42}, {"Ag", 107.87}, {"Cd", 112.41}, {"In", 114.82},
        {"Sn", 118.71}, {"Sb", 121.76}, {"Te", 127.60}, {"I", 126.90},  {"Xe", 131.29}, {"Cs", 132.91}, {"Ba", 137.33},
        {"Hf", 178.49}, {"Ta", 180.95}, {"W", 183.84},  {"Re", 186.21}, {"Os", 190.23}, {"Ir", 192.22}, {"Pt", 195.08},
        {"Au", 196.97}, {"Hg", 200.59}, {"Tl", 204.38}, {"Pb", 207.2},  {"Bi", 208.98}};

    auto it = weights.find(element);
    if (it != weights.end()) {
        return it->second;
    }
    return 0.0;
}

#endif  // ATOMIC_WEIGHTS_H_INCLUDED
