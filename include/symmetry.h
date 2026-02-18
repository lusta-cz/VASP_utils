#ifndef SYMMETRY_H_INCLUDED
#define SYMMETRY_H_INCLUDED

#include <spglib.h>

#include <memory>
#include <optional>

struct POSCAR;

using SpglibDatasetPtr = std::unique_ptr<SpglibDataset, void (*)(SpglibDataset*)>;

SpglibDatasetPtr analyzeSymmetry(const POSCAR& poscar, const double& symprec);
void printSymmetryInfo(const SpglibDataset& dataset, const bool& wyckoff, const bool& symoperation);
void printSymmetryOperations(const SpglibDataset& dataset);
std::optional<POSCAR> makePrimitiveCell(const POSCAR& poscar, const double& symprec);
std::optional<POSCAR> makeConventionalCell(const POSCAR& poscar, const double& symprec);

#endif  // SYMMETRY_H_INCLUDED
