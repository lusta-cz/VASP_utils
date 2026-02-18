#include "symmetry.h"

#include <spglib.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <vector>

#include "poscar_file.h"

SpglibDatasetPtr analyzeSymmetry(const POSCAR& poscar, const double& symprec) {
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
    std::map<std::string, int> element_map;

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

std::optional<POSCAR> makePrimitiveCell(const POSCAR& poscar, const double& symprec) {
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

    // Checking if empty spheres are present in input
    for (const auto& el : poscar.elements) {
        if (el == "X" || el == "E" || el == "V" || el == "Vac") {
            std::cerr << "Warning: Empty sphere detected.\n"
                      << "SPGLIB will treat them as real atoms and symmetry may change.\n";
        }
    }
    // 3) Call spglib primitive finder
    int num_prim = spg_find_primitive(lattice, positions, types, num_atoms, symprec);

    if (num_prim <= 0) {
        return std::nullopt;  // failed
    }

    // Determine maximum type index
    int max_type = 0;
    for (int i = 0; i < num_prim; ++i)
        if (types[i] > max_type)
            max_type = types[i];

    // Create lookup vector (type -> element)
    std::vector<std::string> type_to_element(max_type + 1);

    // Fill from original element_map
    for (const auto& [el, t] : element_map)
        type_to_element[t] = el;

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

    prim.elements = poscar.elements;
    prim.num_atoms.assign(prim.elements.size(), 0);

    for (int i = 0; i < num_prim; ++i) {
        std::string el = type_to_element[types[i]];

        // Find element index in original ordering
        for (size_t j = 0; j < prim.elements.size(); ++j) {
            if (prim.elements[j] == el) {
                prim.num_atoms[j]++;
                break;
            }
        }
    }
    return prim;
}

std::optional<POSCAR> makeConventionalCell(const POSCAR& poscar, const double& symprec) {
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

    // Checking if empty spheres are present in input
    for (const auto& el : poscar.elements) {
        if (el == "X" || el == "E" || el == "V" || el == "Vac") {
            std::cerr << "Warning: Empty sphere detected.\n"
                      << "SPGLIB will treat them as real atoms and symmetry may change.\n";
        }
    }

    // 3) Call spglib standardization
    int num_std = spg_standardize_cell(lattice, positions, types, num_atoms, 0, 1, symprec);

    if (num_std <= 0) {
        return std::nullopt;
    }

    // Determine maximum type index
    int max_type = 0;
    for (int i = 0; i < num_std; ++i)
        if (types[i] > max_type)
            max_type = types[i];

    // Create lookup vector (type -> element)
    std::vector<std::string> type_to_element(max_type + 1);

    // Fill from original element_map
    for (const auto& [el, t] : element_map)
        type_to_element[t] = el;

    // 4) Build POSCAR from arrays now holding primitive cell
    POSCAR std_poscar;
    std_poscar.comment = poscar.comment + " standardized cell";
    std_poscar.is_direct = true;
    std_poscar.total_atoms = num_std;

    // Copy new lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            std_poscar.lattice[i][j] = lattice[i][j];

    // Copy positions
    std_poscar.coordinates.resize(num_std);
    for (int i = 0; i < num_std; ++i) {
        std_poscar.coordinates[i].x = positions[i][0];
        std_poscar.coordinates[i].y = positions[i][1];
        std_poscar.coordinates[i].z = positions[i][2];
    }

    std_poscar.elements = poscar.elements;
    std_poscar.num_atoms.assign(std_poscar.elements.size(), 0);

    for (int i = 0; i < num_std; ++i) {
        std::string el = type_to_element[types[i]];

        // Find element index in original ordering
        for (size_t j = 0; j < std_poscar.elements.size(); ++j) {
            if (std_poscar.elements[j] == el) {
                std_poscar.num_atoms[j]++;
                break;
            }
        }
    }
    return std_poscar;
}

void printSymmetryInfo(const SpglibDataset& dataset, const bool& wyckoff, const bool& symoperation) {
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

    if (symoperation)
        printSymmetryOperations(dataset);

    // Compute number of irreducible atoms
    std::set<int> irreducible_atoms;
    for (int i = 0; i < dataset.n_atoms; ++i)
        irreducible_atoms.insert(dataset.equivalent_atoms[i]);

    std::cout << "\nNumber of Wyckoff positions (irreducible atoms): " << irreducible_atoms.size() << "\n";

    // Print Wyckoff letters for each atom
    if (wyckoff) {
        std::cout << "Wyckoff letters: ";
        for (int i = 0; i < dataset.n_atoms; ++i) {
            char wyckoff_letter = 'a' + dataset.wyckoffs[i];  // convert 0->'a', 1->'b', ...
            std::cout << wyckoff_letter << " ";
        }
    }
    std::cout << "\n";
}

void printSymmetryOperations(const SpglibDataset& dataset) {
    for (int i = 0; i < dataset.n_operations; ++i) {
        std::cout << "Operation " << i + 1 << ":\n";
        std::cout << "  Rotation matrix:\n";
        for (int j = 0; j < 3; ++j)
            std::cout << "   " << std::setw(2) << dataset.rotations[i][j][0] << " " << std::setw(2)
                      << dataset.rotations[i][j][1] << " " << std::setw(2) << dataset.rotations[i][j][2] << "\n";
        std::cout << "  Translation vector: " << dataset.translations[i][0] << " " << dataset.translations[i][1] << " "
                  << dataset.translations[i][2] << "\n";
    }
}
