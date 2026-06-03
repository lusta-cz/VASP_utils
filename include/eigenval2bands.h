#ifndef EIGENVAL2BANDS_H_INCLUDED
#define EIGENVAL2BANDS_H_INCLUDED

#include <string>
#include <vector>

#include "procar_file.h"

struct KPoint {
    double x, y, z;
    std::vector<double> energies_up;
    std::vector<double> energies_dn;  // empty if nspin == 1
};

struct EigenvalData {
    int total_kpoints{};
    int nbands{};
    int nspin{1};
    std::vector<KPoint> kpoints;
};

bool parseFromDoscar(const std::string& filename, double& e_fermi);
bool parseFromOutcar(const std::string& filename, double& e_fermi);
bool parseKpoints(const std::string& filename, int& kpts_per_seg);
bool parseEigenval(const std::string& filename, EigenvalData& data);
bool checkNonCollinear(std::string& filename, bool& is_ncl);
bool parseProcar(const std::string& filename, ProcarData& data);

#endif  // EIGENVAL2BANDS_H_INCLUDED
