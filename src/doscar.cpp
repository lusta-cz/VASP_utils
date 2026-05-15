#include "doscar.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error_parse.h"

namespace {

static std::vector<double> parseLine(const std::string& line) {
    std::istringstream ss(line);
    std::vector<double> row;
    double v;
    while (ss >> v)
        row.push_back(v);
    return row;
}

static std::vector<std::string> makeTdosColNames(int nspin) {
    if (nspin == 1)
        return {"energy[eV]", "dos[states/eV]", "idos[states]"};
    return {"energy[eV]", "dos_up[states/eV]", "dos_dn[states/eV]", "idos_up[states]", "idos_dn[states]"};
}

// Build column name list for partial DOS from the detected column count.
// ncols includes the energy column.
static std::vector<std::string> makePdosColNames(int ncols, int nspin) {
    int data_cols = ncols - 1;
    std::vector<std::string> names{"energy[eV]"};

    // Orbital sets for known LORBIT decomposition levels
    static const std::vector<std::string> spd = {"s", "p", "d"};
    static const std::vector<std::string> spdf = {"s", "p", "d", "f"};
    static const std::vector<std::string> m9 = {"s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2"};
    static const std::vector<std::string> m16 = {"s",     "py",    "pz",   "px",   "dxy", "dyz",  "dz2",  "dxz",
                                                 "x2-y2", "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"};

    const std::vector<std::string>* orbs = nullptr;
    int n = (nspin == 1) ? data_cols : data_cols / 2;

    if (n == 3)
        orbs = &spd;
    else if (n == 4)
        orbs = &spdf;
    else if (n == 9)
        orbs = &m9;
    else if (n == 16)
        orbs = &m16;

    if (!orbs) {
        for (int i = 1; i <= data_cols; ++i)
            names.push_back("col_" + std::to_string(i));
        return names;
    }

    if (nspin == 1) {
        names.insert(names.end(), orbs->begin(), orbs->end());
    } else {
        for (const auto& o : *orbs) {
            names.push_back(o + "_up");
            names.push_back(o + "_dn");
        }
    }
    return names;
}

bool parseDoscarImpl(std::ifstream& file, DoscarData& data, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Line 1: NIONS NKPTS ... (remaining fields ignored)
    if (!readRequiredLine(file, line, line_number, error, "DOSCAR line 1"))
        return false;
    {
        std::istringstream ss(line);
        if (!(ss >> data.nions >> data.nkpts))
            return fail(error, ParseErrorKind::Parse, line_number, "expected NIONS NKPTS on line 1, got: " + line);
    }

    // Lines 2–5: skip
    for (int i = 1; i < 5; ++i)
        if (!readRequiredLine(file, line, line_number, error, "DOSCAR header"))
            return false;

    // Line 6: Emax Emin NEDOS Efermi weight
    if (!readRequiredLine(file, line, line_number, error, "DOSCAR line 6 (Emax Emin NEDOS Efermi weight)"))
        return false;
    {
        std::istringstream ss(line);
        double weight;
        if (!(ss >> data.emax >> data.emin >> data.nedos >> data.e_fermi >> weight))
            return fail(error, ParseErrorKind::Parse, line_number,
                        "expected Emax Emin NEDOS Efermi weight, got: " + line);
    }

    if (data.nedos <= 0)
        return fail(error, ParseErrorKind::Semantic, line_number,
                    "NEDOS must be positive, got: " + std::to_string(data.nedos));

    // ── Total DOS block ───────────────────────────────────────────────────────
    data.tdos.reserve(data.nedos);
    for (int i = 0; i < data.nedos; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "total DOS row"))
            return false;

        auto row = parseLine(line);
        if (row.empty())
            return fail(error, ParseErrorKind::Parse, line_number, "empty total DOS line");

        if (i == 0) {
            int nc = static_cast<int>(row.size());
            if (nc == 3)
                data.nspin = 1;
            else if (nc == 5)
                data.nspin = 2;
            else
                return fail(error, ParseErrorKind::Parse, line_number,
                            "unexpected total DOS column count (" + std::to_string(nc) +
                                "), expected 3 (nspin=1) or 5 (nspin=2)");
            data.tdos_col_names = makeTdosColNames(data.nspin);
        }
        data.tdos.push_back(std::move(row));
    }

    // ── Partial DOS blocks (one per ion) ─────────────────────────────────────
    for (int ia = 0; ia < data.nions; ++ia) {
        // Read the next non-blank line (block header or EOF)
        bool got_line = false;
        while (std::getline(file, line)) {
            ++line_number;
            if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                got_line = true;
                break;
            }
        }
        if (!got_line)
            break;  // clean EOF after tDOS — no partial DOS

        // Validate as block header (requires at least 5 numbers)
        {
            std::istringstream hss(line);
            double v1, v2, v3, v4, v5;
            if (!(hss >> v1 >> v2 >> v3 >> v4 >> v5)) {
                if (ia == 0)
                    break;  // first line after tDOS is not a header — no pDOS
                return fail(error, ParseErrorKind::Parse, line_number,
                            "expected partial DOS block header for atom " + std::to_string(ia + 1) + ", got: " + line);
            }
        }

        std::vector<std::vector<double>> atom_block;
        atom_block.reserve(data.nedos);
        for (int i = 0; i < data.nedos; ++i) {
            if (!readRequiredLine(file, line, line_number, error, "partial DOS atom " + std::to_string(ia + 1)))
                return false;

            auto row = parseLine(line);
            if (row.empty())
                return fail(error, ParseErrorKind::Parse, line_number,
                            "empty partial DOS line for atom " + std::to_string(ia + 1));

            if (ia == 0 && i == 0)
                data.pdos_col_names = makePdosColNames(static_cast<int>(row.size()), data.nspin);

            atom_block.push_back(std::move(row));
        }
        data.pdos.push_back(std::move(atom_block));
    }

    data.has_partial = !data.pdos.empty();

    if (data.has_partial && static_cast<int>(data.pdos.size()) != data.nions)
        std::cerr << "Warning: expected " << data.nions << " partial DOS blocks, found " << data.pdos.size() << "\n";

    return true;
}

}  // namespace

bool parseDoscar(const std::string& filename, DoscarData& data) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;
    if (!parseDoscarImpl(file, data, error)) {
        reportParseError(filename, error);
        return false;
    }
    return true;
}
