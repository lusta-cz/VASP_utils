#ifndef SURFACE_H_INCLUDED
#define SURFACE_H_INCLUDED

#include "poscar_file.h"

struct SlabOptions {
    int h{0}, k{0}, l{0};  ///< Miller indices defining the surface normal.
    int n_layers{0};       ///< Number of atomic layers to include.
    int n_frozen{0};       ///< Number of bottom layers frozen (F F F).
    double vacuum{15.0};   ///< Vacuum thickness to add above the slab (Å).
    double symprec{1e-5};  ///< Symmetry tolerance for primitive-cell conversion.
};

/**
 * @brief Build a surface slab POSCAR from a bulk structure.
 *
 * The slab is oriented so that the (hkl) plane is perpendicular to the
 * c-axis.  Vacuum is added on top; the bottom @p n_frozen layers are
 * frozen with Selective Dynamics F F F.  If the in-plane surface unit
 * cell is larger than the projected bulk cell, an informational message
 * is printed to stdout.
 *
 * @param bulk    Bulk POSCAR (Direct or Cartesian — converted internally).
 * @param opts    Slab parameters.
 * @return        New POSCAR representing the slab, or an empty POSCAR on error.
 */
POSCAR buildSlab(const POSCAR& bulk, const SlabOptions& opts);

#endif  // SURFACE_H_INCLUDED
