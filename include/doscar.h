#ifndef DOSCAR_H_INCLUDED
#define DOSCAR_H_INCLUDED

#include <string>
#include <vector>

struct DoscarData {
    int nions{};
    int nkpts{};
    int nedos{};
    int nspin{1};
    double emax{};
    double emin{};
    double e_fermi{};
    bool has_partial{false};

    // [nedos][ncols]: first col is energy (raw, unshifted)
    std::vector<std::vector<double>> tdos;
    std::vector<std::string> tdos_col_names;

    // [nions][nedos][ncols]: first col is energy (raw, unshifted)
    std::vector<std::vector<std::vector<double>>> pdos;
    std::vector<std::string> pdos_col_names;
};

bool parseDoscar(const std::string& filename, DoscarData& data);

#endif  // DOSCAR_H_INCLUDED
