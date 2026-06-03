#ifndef PROCAR_FILE_H_INCLUDED
#define PROCAR_FILE_H_INCLUDED

#include <string>
#include <vector>

// Holds orbital projection weights for a single band at a specific k-point
struct BandProjection {
    // Outer vector: over atoms (0 to NIONS-1)
    // Inner vector: over orbitals (e.g., s, py, pz, px, dxy...) plus the last element which is the total
    std::vector<std::vector<double>> atom_weights;
};

// Holds projections for all bands at a specific k-point
struct KPointProjection {
    std::vector<BandProjection> bands;  // Size matches data.nbands
};

struct ProcarData {
    int num_kpoints{0};
    int num_bands{0};
    int num_ions{0};

    // List of orbital labels parsed from the header line (e.g., "s", "py", "pz", "px"...)
    std::vector<std::string> orbital_labels;

    // Master matrix: [kpoint_idx]
    std::vector<KPointProjection> kpoints;
};

#endif  // PROCAR_FILE_H_INCLUDED
