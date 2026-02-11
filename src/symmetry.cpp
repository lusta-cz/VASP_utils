#include "symmetry.h"
#include "poscar_file.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <spglib.h>
#include <map>
#include <optional>
#include <set>

SpglibDatasetPtr analyzeSymmetry(const POSCAR& poscar,  const double& symprec)
{
    // Make a copy in fractional coordinates
    POSCAR poscarDirect = poscar;
    if (!poscarDirect.is_direct) {
        poscarDirect.toDirect();
    }

    int num_atoms = static_cast<int>(poscarDirect.total_atoms);

    // Prepare lattice as flat 3x3 vector
    double lattice[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[i][j] = poscarDirect.lattice[i][j];

    // Prepare atomic positions as flat vector (num_atoms*3)
    double positions[num_atoms][3];
    for (int i = 0; i < num_atoms; ++i) {
        positions[i][0] = poscarDirect.coordinates[i].x;
        positions[i][1] = poscarDirect.coordinates[i].y;
        positions[i][2] = poscarDirect.coordinates[i].z;
    }

    // Prepare atomic types
    int types[num_atoms];
    int type_counter = 1;
    std::map<std::string,int> element_map;

    for (size_t i = 0, idx = 0; i < poscarDirect.elements.size(); ++i) {
        const std::string& el = poscarDirect.elements[i];
        int n = poscarDirect.num_atoms[i];
        if (element_map.find(el) == element_map.end())
            element_map[el] = type_counter++;
        for (int j = 0; j < n; ++j)
            types[idx++] = element_map[el];
    }

    SpglibDataset* dataset = spg_get_dataset(lattice, positions, types, num_atoms, symprec);

    return SpglibDatasetPtr(dataset, &spg_free_dataset);
}

std::optional<POSCAR> makePrimitiveCell(const POSCAR& poscar, const double& symprec)
{
    // 1) Copy and ensure fractional coordinates
    POSCAR poscarDirect = poscar;
    if (!poscarDirect.is_direct) {
        poscarDirect.toDirect();
    }

    int num_atoms = static_cast<int>(poscarDirect.total_atoms);

    // 2) Prepare C‑arrays for spg_find_primitive
    double lattice[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[i][j] = poscarDirect.lattice[i][j];

    double positions[num_atoms][3];
    for (int i = 0; i < num_atoms; ++i) {
        positions[i][0] = poscarDirect.coordinates[i].x;
        positions[i][1] = poscarDirect.coordinates[i].y;
        positions[i][2] = poscarDirect.coordinates[i].z;
    }

    int types[num_atoms];
    int type_counter = 1;
    std::map<std::string, int> element_map;
    for (size_t i = 0, idx = 0; i < poscarDirect.elements.size(); ++i) {
        const std::string& el = poscarDirect.elements[i];
        int count = poscarDirect.num_atoms[i];
        if (element_map.find(el) == element_map.end())
            element_map[el] = type_counter++;
        for (int j = 0; j < count; ++j)
            types[idx++] = element_map[el];
    }

    // 3) Call spglib primitive finder
    int num_prim = spg_find_primitive(
        lattice,
        positions,
        types,
        num_atoms,
        symprec
    );

    if (num_prim <= 0) {
        return std::nullopt; // failed
    }

    // 4) Build POSCAR from arrays now holding primitive cell
    POSCAR prim;
    prim.comment = poscar.comment + " primitive cell";
    prim.is_direct = true;
    prim.total_atoms = num_prim;

    // Copy new lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            prim.lattice[i][j] = lattice[i][j];

    // Copy positions
    prim.coordinates.resize(num_prim);
    for (int i = 0; i < num_prim; ++i) {
        prim.coordinates[i].x = positions[i][0];
        prim.coordinates[i].y = positions[i][1];
        prim.coordinates[i].z = positions[i][2];
    }

    // Map types back to element symbols
    std::map<int, std::string> reverse_map;
    for (const auto& [el, t] : element_map)
        reverse_map[t] = el;

    std::map<std::string, int> elem_counts;
    for (int i = 0; i < num_prim; ++i)
        elem_counts[reverse_map[ types[i] ]]++;

    prim.elements.clear();
    prim.num_atoms.clear();
    for (auto& [el, cnt] : elem_counts) {
        prim.elements.push_back(el);
        prim.num_atoms.push_back(cnt);
    }

    return prim;
}



void printSymmetryInfo(const SpglibDataset& dataset, const bool& wyckoff, const bool& symoperation)
{
    std::cout << "=== Symmetry Information ===\n";

    // Space group
    std::cout << "Space group number: " << dataset.spacegroup_number << "\n";
    std::cout << "International symbol: " << dataset.international_symbol << "\n";

    // Hall symbol
    std::cout << "Hall symbol: " << dataset.hall_symbol << "\n";

    // Point group
    std::cout << "Point group: " << dataset.pointgroup_symbol << "\n";

    // Symmetry operations
    std::cout << "\nNumber of symmetry operations: " << dataset.n_operations << "\n";

    if(symoperation)
        printSymmetryOperations(dataset);

        // Compute number of irreducible atoms
    std::set<int> irreducible_atoms;
    for (int i = 0; i < dataset.n_atoms; ++i)
        irreducible_atoms.insert(dataset.equivalent_atoms[i]);

    std::cout << "\nNumber of Wyckoff positions (irreducible atoms): "
              << irreducible_atoms.size() << "\n";

    // Print Wyckoff letters for each atom
    if(wyckoff){
        std::cout << "Wyckoff letters: ";
        for (int i = 0; i < dataset.n_atoms; ++i) {
            char wyckoff_letter = 'a' + dataset.wyckoffs[i]; // convert 0->'a', 1->'b', ...
            std::cout << wyckoff_letter << " ";
        }
    }
    std::cout << "\n";
}

void printSymmetryOperations(const SpglibDataset& dataset){
    for (int i = 0; i < dataset.n_operations; ++i)
    {
        std::cout << "Operation " << i+1 << ":\n";
        std::cout << "  Rotation matrix:\n";
        for (int j = 0; j < 3; ++j)
            std::cout << "   "
                    << std::setw(2) << dataset.rotations[i][j][0] << " "
                    << std::setw(2) << dataset.rotations[i][j][1] << " "
                    << std::setw(2) << dataset.rotations[i][j][2] << "\n";
        std::cout << "  Translation vector: "
                << dataset.translations[i][0] << " "
                << dataset.translations[i][1] << " "
                << dataset.translations[i][2] << "\n";
    }

}

