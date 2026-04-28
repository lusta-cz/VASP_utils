#ifndef ELASTIC_H_INCLUDED
#define ELASTIC_H_INCLUDED

#include <array>
#include <optional>
#include <string>
#include <vector>

struct POSCAR;

/**
 * @brief Crystal system classification used to select minimal strain modes.
 *
 * TetragonalI  = point groups 4, -4, 4/m             (7 independent constants, has C16)
 * TetragonalII = point groups 422, 4mm, -42m, 4/mmm  (6 independent constants)
 * TrigonalI    = point groups 3, -3                  (7 independent constants, has C14, C25)
 * TrigonalII   = point groups 32, 3m, -3m            (6 independent constants, has C14)
 */
enum class CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    TetragonalI,
    TetragonalII,
    TrigonalI,
    TrigonalII,
    Hexagonal,
    Cubic
};

/**
 * @brief One independent strain pattern for elastic constant extraction.
 *
 * @c voigt holds the unit-amplitude Voigt strain vector [ε1,ε2,ε3,ε4,ε5,ε6].
 * The actual deformation amplitude is multiplied in at generation time.
 * Voigt convention: ε4 = 2ε_yz, ε5 = 2ε_xz, ε6 = 2ε_xy.
 */
struct ElasticStrainMode {
    std::string label;            ///< Human-readable label, e.g. "C11-C12" or "e1+e2".
    std::array<double, 6> voigt;  ///< Unit Voigt strain vector.
};

/**
 * @brief Manifest written by poscar_elastic_deform and consumed by vasp_elastic.
 *
 * Relative paths in @c dirs are relative to the directory that contains the log file.
 */
struct ElasticDeformLog {
    std::string method;  ///< "energy" or "stress"
    double amplitude{0.01};
    int npoints{6};             ///< Number of non-zero deformed structures per mode.
    std::string source_hint;    ///< "vasprun", "outcar", or "oszicar"
    std::string symmetry_mode;  ///< "auto" or "none"
    CrystalSystem crystal_system{CrystalSystem::Triclinic};
    int space_group_number{1};
    std::string space_group_symbol;
    std::string point_group;
    int n_independent{21};
    std::string reference_dir;  ///< Relative path to the reference (zero-strain) directory.
    std::vector<ElasticStrainMode> modes;
    std::vector<std::vector<double>> amplitudes;  ///< amplitudes[mode][amp_idx]
    std::vector<std::vector<std::string>> dirs;   ///< dirs[mode][amp_idx], relative paths
};

/** @brief Map spglib space group number (1–230) to crystal system. */
CrystalSystem crystalSystemFromSpaceGroup(int spg_number);

/** @brief Human-readable name for a crystal system, e.g. "hexagonal". */
std::string crystalSystemName(CrystalSystem cs);

/** @brief Number of independent elastic constants for a crystal system. */
int nIndependentConstants(CrystalSystem cs);

/**
 * @brief Minimal set of strain modes for the energy–strain method.
 *
 * Returns the minimal independent strain patterns for the given crystal system
 * following Le Page & Saxe, PRB 65, 104104 (2002).
 */
std::vector<ElasticStrainMode> energyStrainModes(CrystalSystem cs);

/**
 * @brief Six Voigt basis strain modes for the stress–strain method.
 *
 * The stress–strain method always uses the same six modes regardless of
 * crystal symmetry; symmetry constraints are imposed on the solution.
 */
std::vector<ElasticStrainMode> stressStrainModes();

/**
 * @brief Apply a symmetric strain to a POSCAR's lattice vectors.
 *
 * Builds the symmetric deformation matrix F = I + ε from the Voigt vector
 * scaled by @p amplitude and transforms each lattice row vector.  Fractional
 * atomic coordinates are unchanged (they co-deform with the lattice).
 *
 * Using the fully symmetric form (off-diagonal elements split equally between
 * upper and lower triangle) avoids any rigid-body rotation, which is essential
 * for noncollinear magnetic calculations where spin directions are fixed in
 * Cartesian space.
 *
 * @param poscar     Input structure (not modified).
 * @param voigt      Unit Voigt strain vector [ε1..ε6].
 * @param amplitude  Scaling factor applied to @p voigt.
 * @return           New POSCAR with deformed lattice.
 */
POSCAR applyStrain(const POSCAR& poscar, const std::array<double, 6>& voigt, double amplitude);

/**
 * @brief Write an elastic deformation manifest to disk.
 * @return false on I/O error.
 */
bool writeElasticLog(const std::string& path, const ElasticDeformLog& log);

/**
 * @brief Read an elastic deformation manifest from disk.
 * @return nullopt on parse or I/O error.
 */
std::optional<ElasticDeformLog> readElasticLog(const std::string& path);

#endif  // ELASTIC_H_INCLUDED
