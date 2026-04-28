#include "eigenval2bands.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "error_parse.h"

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

}  // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

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
