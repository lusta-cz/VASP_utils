#include "eigenval2bands.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "error_parse.h"
#include "procar_file.h"

namespace {

// ---------------------------------------------------------------------------
// DOSCAR
// ---------------------------------------------------------------------------

bool parseDoscarFermi(std::ifstream& file, double& e_fermi, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Skip lines 1-5
    for (int i = 0; i < 5; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "DOSCAR header"))
            return false;
    }

    // Line 6: Emax  Emin  NEDOS  Efermi  weight
    if (!readRequiredLine(file, line, line_number, error, "DOSCAR line 6 (Emax Emin NEDOS Efermi weight)"))
        return false;

    std::istringstream ss(line);
    double emax, emin, nedos_d;
    if (!(ss >> emax >> emin >> nedos_d >> e_fermi))
        return fail(error, ParseErrorKind::Parse, line_number,
                    "expected 4 values (Emax Emin NEDOS Efermi), got: " + line);

    return true;
}

// ---------------------------------------------------------------------------
// OUTCAR
// ---------------------------------------------------------------------------

bool checkNonCollinear(std::ifstream& file, bool& is_ncl, ParseError& error) {
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("LNONCOLLINEAR") != std::string::npos) {
            if (line.find("T") != std::string::npos) {
                is_ncl = true;
            }
            break;
        }
        if (line.find("ENTER initialization of structure-lists") != std::string::npos) {
            break;
        }
    }
    return true;
}

bool parseOutcarFermi(std::ifstream& file, double& e_fermi, ParseError& error) {
    std::string line;
    int line_number = 0;
    bool found = false;

    // Take the last occurrence — VASP writes E-fermi after every SCF step
    while (std::getline(file, line)) {
        ++line_number;
        if (line.find("E-fermi") != std::string::npos) {
            std::istringstream ss(line);
            std::string t1, t2;
            double val;
            if (ss >> t1 >> t2 >> val) {
                e_fermi = val;
                found = true;
            }
        }
    }

    if (!found)
        return fail(error, ParseErrorKind::Parse, 0, "'E-fermi' entry not found");

    return true;
}

// ---------------------------------------------------------------------------
// KPOINTS
// ---------------------------------------------------------------------------

bool parseKpointsFile(std::ifstream& file, int& kpts_per_seg, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Line 1: comment — skip
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS comment line"))
        return false;

    // Line 2: k-points per segment
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS k-points per segment"))
        return false;

    std::istringstream ss(line);
    if (!(ss >> kpts_per_seg))
        return fail(error, ParseErrorKind::Parse, line_number, "expected integer k-points per segment, got: " + line);

    if (kpts_per_seg <= 1)
        return fail(error, ParseErrorKind::Semantic, line_number, "k-points per segment must be greater than 1");

    // Line 3: must start with 'L'/'l' to confirm line mode
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS mode line"))
        return false;

    if (line.empty() || std::toupper(static_cast<unsigned char>(line[0])) != 'L')
        return fail(error, ParseErrorKind::Semantic, line_number,
                    "KPOINTS file is not in line mode (line 3 must start with 'L'), got: " + line);

    return true;
}

// ---------------------------------------------------------------------------
// EIGENVAL
// ---------------------------------------------------------------------------

bool parseEigenvalFile(std::ifstream& file, EigenvalData& data, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Skip lines 1-5
    for (int i = 0; i < 5; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "EIGENVAL header"))
            return false;
    }

    // Line 6: NIONS  NKPTS  NBANDS
    if (!readRequiredLine(file, line, line_number, error, "EIGENVAL line 6 (NIONS NKPTS NBANDS)"))
        return false;

    {
        std::istringstream ss(line);
        int nions;
        if (!(ss >> nions >> data.total_kpoints >> data.nbands))
            return fail(error, ParseErrorKind::Parse, line_number, "expected NIONS NKPTS NBANDS, got: " + line);
    }

    if (data.total_kpoints <= 0)
        return fail(error, ParseErrorKind::Semantic, line_number, "total k-points must be positive");
    if (data.nbands <= 0)
        return fail(error, ParseErrorKind::Semantic, line_number, "number of bands must be positive");

    data.kpoints.reserve(data.total_kpoints);
    bool nspin_detected = false;

    for (int ikpt = 0; ikpt < data.total_kpoints; ++ikpt) {
        // Blank separator
        if (!readRequiredLine(file, line, line_number, error, "k-point separator"))
            return false;

        // k-point: kx ky kz weight
        if (!readRequiredLine(file, line, line_number, error, "k-point coordinates"))
            return false;

        KPoint kp;
        {
            std::istringstream ss(line);
            double weight;
            if (!(ss >> kp.x >> kp.y >> kp.z >> weight))
                return fail(error, ParseErrorKind::Parse, line_number,
                            "expected kx ky kz weight at k-point " + std::to_string(ikpt + 1));
        }

        kp.energies_up.reserve(data.nbands);
        kp.energies_dn.reserve(data.nbands);

        for (int ib = 0; ib < data.nbands; ++ib) {
            if (!readRequiredLine(file, line, line_number, error, "band eigenvalue"))
                return false;

            std::istringstream ss(line);
            int bidx;
            double e_up;
            if (!(ss >> bidx >> e_up))
                return fail(error, ParseErrorKind::Parse, line_number,
                            "expected band index and eigenvalue at k-point " + std::to_string(ikpt + 1) + ", band " +
                                std::to_string(ib + 1));

            kp.energies_up.push_back(e_up);

            // Detect nspin from first band line of first k-point
            if (!nspin_detected && ikpt == 0 && ib == 0) {
                double e_dn, occ;
                data.nspin = (ss >> e_dn >> occ) ? 2 : 1;
                nspin_detected = true;
            }

            if (data.nspin == 2) {
                std::istringstream ss2(line);
                int tmp;
                double eu, ed;
                if (!(ss2 >> tmp >> eu >> ed))
                    return fail(error, ParseErrorKind::Parse, line_number,
                                "expected spin-down eigenvalue at k-point " + std::to_string(ikpt + 1) + ", band " +
                                    std::to_string(ib + 1));
                kp.energies_dn.push_back(ed);
            }
        }

        data.kpoints.push_back(std::move(kp));
    }

    return true;
}

// ---------------------------------------------------------------------------
// PROCAR
// ---------------------------------------------------------------------------

bool parseProcarFile(std::ifstream& file, ProcarData& data, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Line 1: Header/Description text - Skip
    if (!readRequiredLine(file, line, line_number, error, "PROCAR header"))
        return false;

    // Line 2: # of k-points: X # of bands: Y # of ions: Z
    if (!readRequiredLine(file, line, line_number, error, "PROCAR dimensions line"))
        return false;

    {
        // Safely extract dimensions using text token filtering
        size_t k_pos = line.find("k-points:");
        size_t b_pos = line.find("bands:");
        size_t i_pos = line.find("ions:");

        if (k_pos == std::string::npos || b_pos == std::string::npos || i_pos == std::string::npos) {
            return fail(error, ParseErrorKind::Parse, line_number,
                        "Invalid PROCAR format. Dimensions line missing keywords: " + line);
        }

        std::istringstream ss_k(line.substr(k_pos + 9));
        ss_k >> data.num_kpoints;
        std::istringstream ss_b(line.substr(b_pos + 6));
        ss_b >> data.num_bands;
        std::istringstream ss_i(line.substr(i_pos + 5));
        ss_i >> data.num_ions;
    }

    if (data.num_kpoints <= 0 || data.num_bands <= 0 || data.num_ions <= 0) {
        return fail(error, ParseErrorKind::Semantic, line_number, "Dimensions must be positive values.");
    }

    data.kpoints.resize(data.num_kpoints);
    bool labels_captured = false;

    // Loop through each K-point block
    for (int ikpt = 0; ikpt < data.num_kpoints; ++ikpt) {
        data.kpoints[ikpt].bands.resize(data.num_bands);

        // Advance through file until we locate the exact string flag for the K-point header
        while (true) {
            if (!readRequiredLine(file, line, line_number, error, "K-point search scope"))
                return false;
            if (line.find("k-point") != std::string::npos)
                break;
        }

        // Loop through each Band block inside this K-point
        for (int ib = 0; ib < data.num_bands; ++ib) {
            // Find band line flag
            while (true) {
                if (!readRequiredLine(file, line, line_number, error, "Band search scope"))
                    return false;
                if (line.find("band") != std::string::npos)
                    break;
            }

            // Next line is the orbital character column header labels (ion  s  py  pz  px...)
            if (!readRequiredLine(file, line, line_number, error, "Orbital labels line"))
                return false;

            // Extract the headers on the fly the first time we see them
            if (!labels_captured) {
                std::istringstream ss_labels(line);
                std::string token;
                ss_labels >> token;  // Skip "ion" tag
                while (ss_labels >> token) {
                    data.orbital_labels.push_back(token);
                }
                labels_captured = true;
            }

            int num_orbitals = data.orbital_labels.size();
            auto& band_proj = data.kpoints[ikpt].bands[ib];
            band_proj.atom_weights.resize(data.num_ions);

            // Loop to parse weights for every individual ion line
            for (int ion = 0; ion < data.num_ions; ++ion) {
                if (!readRequiredLine(file, line, line_number, error, "Ion matrix data row"))
                    return false;

                std::istringstream ss_vals(line);
                int ion_idx_check;
                if (!(ss_vals >> ion_idx_check)) {
                    return fail(error, ParseErrorKind::Parse, line_number, "Expected ion row index: " + line);
                }

                band_proj.atom_weights[ion].resize(num_orbitals);
                for (int iorb = 0; iorb < num_orbitals; ++iorb) {
                    if (!(ss_vals >> band_proj.atom_weights[ion][iorb])) {
                        return fail(
                            error, ParseErrorKind::Parse, line_number,
                            "Failed to parse orbital projection index " + std::to_string(iorb) + " on line: " + line);
                    }
                }
            }

            // Skip the trailing total sum row summary line that VASP prints after the individual ions
            if (!readRequiredLine(file, line, line_number, error, "Skip total projection summary row"))
                return false;
        }
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

bool checkNonCollinear(std::string& filename, bool& is_ncl) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: cannot open file '" << filename << "' for LNONCOLLINEAR check \n";
        return false;
    }
    ParseError error;
    if (!checkNonCollinear(file, is_ncl, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}

bool parseFromDoscar(const std::string& filename, double& e_fermi) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseDoscarFermi(file, e_fermi, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}

bool parseFromOutcar(const std::string& filename, double& e_fermi) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseOutcarFermi(file, e_fermi, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}

bool parseKpoints(const std::string& filename, int& kpts_per_seg) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseKpointsFile(file, kpts_per_seg, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}

bool parseEigenval(const std::string& filename, EigenvalData& data) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseEigenvalFile(file, data, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}

bool parseProcar(const std::string& filename, ProcarData& data) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseProcarFile(file, data, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}
