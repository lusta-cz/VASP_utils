#ifndef SYMMETRY_H_INCLUDED
#define SYMMETRY_H_INCLUDED

#include <memory>
#include <optional>

#include <spglib.h>

struct POSCAR;



using SpglibDatasetPtr = std::unique_ptr<SpglibDataset, void(*)(SpglibDataset*)>;

SpglibDatasetPtr analyzeSymmetry(const POSCAR& poscar, const double& symprec);
std::optional<POSCAR> makePrimitiveCell(const POSCAR& poscar, const double& symprec);
void printSymmetryInfo(const SpglibDataset& dataset, const bool& wyckoff, const bool& symoperation);
void printSymmetryOperations(const SpglibDataset& dataset);

#endif // SYMMETRY_H_INCLUDED
