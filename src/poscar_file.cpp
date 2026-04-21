#include "poscar_file.h"

#include "random_utility.h"

// Linear algebra

#include <cblas.h>
#include <lapacke.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

enum class PoscarErrorKind { Io, Parse, Semantic };

struct PoscarError {
    PoscarErrorKind kind{PoscarErrorKind::Parse};
    int line_number{0};
    std::string message;
};

bool fail(PoscarError& error, PoscarErrorKind kind, int line_number, std::string message) {
    error.kind = kind;
    error.line_number = line_number;
    error.message = std::move(message);
    return false;
}

void reportReadError(const std::string& filename, const PoscarError& error) {
    const char* category = "Parse error";
    if (error.kind == PoscarErrorKind::Io) {
        category = "I/O error";
    } else if (error.kind == PoscarErrorKind::Semantic) {
        category = "Semantic error";
    }

    std::cerr << category << " while reading " << filename;
    if (error.line_number > 0) {
        std::cerr << " at line " << error.line_number;
    }
    std::cerr << ": " << error.message << "\n";
}

bool readRequiredLine(std::ifstream& file, std::string& line, int& line_number, PoscarError& error,
                      const std::string& context) {
    if (!std::getline(file, line)) {
        return fail(error, file.bad() ? PoscarErrorKind::Io : PoscarErrorKind::Parse, line_number + 1,
                    "unexpected end of file while reading " + context);
    }
    ++line_number;
    return true;
}

bool parseSingleDouble(const std::string& line, double& value) {
    std::istringstream iss(line);
    iss >> value;
    if (!iss)
        return false;

    iss >> std::ws;
    return iss.eof();
}

template <typename T>
bool parseVectorLine(const std::string& line, std::vector<T>& out) {
    std::istringstream iss(line);
    T value{};
    while (iss >> value) {
        out.push_back(value);
    }

    if (iss.fail() && !iss.eof()) {
        return false;
    }

    return true;
}

bool parseLatticeRow(const std::string& line, double row[3]) {
    std::istringstream iss(line);
    if (!(iss >> row[0] >> row[1] >> row[2])) {
        return false;
    }

    iss >> std::ws;
    return iss.eof();
}

bool parseCoordinateRow(const std::string& line, bool selective_dynamics, Atom& atom) {
    std::istringstream iss(line);
    if (!(iss >> atom.x >> atom.y >> atom.z)) {
        return false;
    }

    atom.selective_flags = {true, true, true};
    if (selective_dynamics) {
        std::string fx, fy, fz;
        if (!(iss >> fx >> fy >> fz)) {
            return false;
        }

        auto parseFlag = [](const std::string& token, bool& value) {
            if (token.empty()) {
                return false;
            }
            const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
            if (c == 'T') {
                value = true;
                return true;
            }
            if (c == 'F') {
                value = false;
                return true;
            }
            return false;
        };

        if (!parseFlag(fx, atom.selective_flags[0]) || !parseFlag(fy, atom.selective_flags[1]) ||
            !parseFlag(fz, atom.selective_flags[2])) {
            return false;
        }
    }

    return true;
}

bool parsePoscarFile(std::ifstream& file, POSCAR& parsed, PoscarError& error) {
    std::string line;
    int line_number = 0;

    if (!readRequiredLine(file, parsed.comment, line_number, error, "comment line")) {
        return false;
    }

    if (!readRequiredLine(file, line, line_number, error, "scaling factor")) {
        return false;
    }
    if (!parseSingleDouble(line, parsed.scale)) {
        return fail(error, PoscarErrorKind::Parse, line_number, "invalid scaling factor");
    }
    if (parsed.scale < 0.0) {
        return fail(error, PoscarErrorKind::Semantic, line_number, "scaling factor must not be negative");
    }
    if (parsed.scale == 0.0) {
        return fail(error, PoscarErrorKind::Semantic, line_number, "scaling factor must not be zero");
    }

    for (int i = 0; i < 3; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "lattice vector")) {
            return false;
        }
        if (!parseLatticeRow(line, parsed.lattice[i])) {
            return fail(error, PoscarErrorKind::Parse, line_number, "invalid lattice vector row");
        }
    }

    if (!readRequiredLine(file, line, line_number, error, "element symbols")) {
        return false;
    }
    if (!parseVectorLine(line, parsed.elements) || parsed.elements.empty()) {
        return fail(error, PoscarErrorKind::Parse, line_number, "invalid or empty element symbols line");
    }

    if (!readRequiredLine(file, line, line_number, error, "atom counts")) {
        return false;
    }
    if (!parseVectorLine(line, parsed.num_atoms) || parsed.num_atoms.empty()) {
        return fail(error, PoscarErrorKind::Parse, line_number, "invalid or empty atom counts line");
    }
    if (parsed.elements.size() != parsed.num_atoms.size()) {
        return fail(error, PoscarErrorKind::Semantic, line_number,
                    "number of element symbols does not match number of atom counts");
    }

    parsed.total_atoms = 0;
    for (int count : parsed.num_atoms) {
        if (count < 0) {
            return fail(error, PoscarErrorKind::Semantic, line_number, "atom counts must not be negative");
        }
        if (parsed.total_atoms > std::numeric_limits<int>::max() - count) {
            return fail(error, PoscarErrorKind::Semantic, line_number, "total atom count overflows int");
        }
        parsed.total_atoms += count;
    }

    if (!readRequiredLine(file, line, line_number, error, "coordinate mode")) {
        return false;
    }
    if (line.empty()) {
        return fail(error, PoscarErrorKind::Parse, line_number, "coordinate mode line must not be empty");
    }

    parsed.selective_dynamics = false;
    if (line[0] == 'S' || line[0] == 's') {
        parsed.selective_dynamics = true;
        if (!readRequiredLine(file, line, line_number, error, "coordinate mode")) {
            return false;
        }
        if (line.empty()) {
            return fail(error, PoscarErrorKind::Parse, line_number, "coordinate mode line must not be empty");
        }
    }

    if (line[0] == 'D' || line[0] == 'd') {
        parsed.is_direct = true;
    } else if (line[0] == 'C' || line[0] == 'c' || line[0] == 'K' || line[0] == 'k') {
        parsed.is_direct = false;
    } else {
        return fail(error, PoscarErrorKind::Semantic, line_number, "coordinate mode must be Direct or Cartesian");
    }

    parsed.coordinates.resize(parsed.total_atoms);
    for (int i = 0; i < parsed.total_atoms; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "atomic coordinates")) {
            return false;
        }
        if (!parseCoordinateRow(line, parsed.selective_dynamics, parsed.coordinates[static_cast<size_t>(i)])) {
            return fail(error, PoscarErrorKind::Parse, line_number,
                        "failed to parse coordinates for atom " + std::to_string(i));
        }
    }

    return true;
}

}  // namespace

void POSCAR::setScaleTo1() {
    if (std::abs(scale - 1.0) < 1e-12)
        return;

    // Scale lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[i][j] *= scale;

    // Scale Cartesian positions only
    if (!is_direct) {
        for (auto& atom : coordinates) {
            atom.x *= scale;
            atom.y *= scale;
            atom.z *= scale;
        }
    }

    scale = 1.0;
}

bool POSCAR::readPOSCAR(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error while reading " << filename << ": cannot open file\n";
        return false;
    }

    POSCAR parsed{};
    PoscarError error;
    if (!parsePoscarFile(file, parsed, error)) {
        reportReadError(filename, error);
        return false;
    }

    *this = std::move(parsed);
    setScaleTo1();
    return true;
}

bool POSCAR::writeCtrlsFile(const std::string& filenameOut) {
    if (is_direct) {
        std::cerr << "Error: Only POSCAR in cartesian coordinates can be writen to ctrls format!\n";
        return false;
    }

    std::ifstream fileTest(filenameOut);
    if (fileTest.good()) {
        std::cerr << "Warning: file \"" << filenameOut << "\" already exists and will be overwritten.\n";
    }

    fileTest.close();

    std::ofstream file(filenameOut);
    if (!file) {
        std::cerr << "Error: cannot create file " << filenameOut << "\n";
        return false;
    }

    // Transform units from angstroms to Bohr
    const double ang_to_bohr = 1.889726125;

    // Writing lattice information

    file << "STRUC    ALAT=" << scale * ang_to_bohr << "\n";

    file << "         PLAT=";
    for (int i = 0; i < 3; ++i) {
        if (i > 0) {
            file << "               ";
        } else {
            file << " ";  // Mezera pro první řádek hned za "PLAT="
        }
        file << std::fixed << std::setprecision(10) << lattice[i][0] << " " << lattice[i][1] << " " << lattice[i][2]
             << " \n";
    }

    // Writing atomic position and element
    file << "SITE \n";
    int atom_idx = 0;

    for (size_t i = 0; i < elements.size(); ++i) {
        for (int j = 0; j < num_atoms[i]; ++j) {
            if (atom_idx >= total_atoms)
                break;

            file << "      ATOM=" << std::left << std::setw(3) << elements[i] << " " << " POS=" << std::fixed
                 << std::setprecision(10) << coordinates[atom_idx].x * ang_to_bohr << " "
                 << coordinates[atom_idx].y * ang_to_bohr << " " << coordinates[atom_idx].z * ang_to_bohr << "\n";

            atom_idx++;
        }
    }

    file.close();

    return true;
}

bool POSCAR::writePOSCAR(const std::string& filenameOut) {
    std::ifstream fileTest(filenameOut);
    if (fileTest.good()) {
        std::cerr << "Warning: file \"" << filenameOut << "\" already exists and will be overwritten.\n";
    }

    std::ofstream file(filenameOut);
    if (!file) {
        std::cerr << "Error: cannot create file " << filenameOut << "\n";
        return false;
    }

    // Writing Line 1: comment
    file << comment << "\n";

    // Writing Line 2: scale
    file << std::fixed << std::setprecision(10) << scale << "\n";

    // Writing Lines 3-5: lattice vectors
    for (int i = 0; i < 3; ++i)
        file << std::fixed << std::setprecision(10) << lattice[i][0] << " " << lattice[i][1] << " " << lattice[i][2]
             << "\n";

    // Writing Line 6: element symbols
    for (size_t i = 0; i < elements.size(); ++i)
        file << elements[i] << " ";
    file << "\n";

    // Writing Line 7: number of atoms
    for (size_t i = 0; i < num_atoms.size(); ++i)
        file << num_atoms[i] << " ";
    file << "\n";

    // Writing Optional: selective dynamics
    if (selective_dynamics)
        file << "Selective Dynamics\n";

    // Writing Direct/Cartesian
    file << (is_direct ? "Direct" : "Cartesian") << "\n";

    // Writing Atomic coordinates
    for (size_t i = 0; i < coordinates.size(); ++i) {
        file << std::fixed << std::setprecision(10) << coordinates[i].x << " " << coordinates[i].y << " "
             << coordinates[i].z;
        if (selective_dynamics) {
            file << " " << (coordinates[i].selective_flags[0] ? "T" : "F") << " "
                 << (coordinates[i].selective_flags[1] ? "T" : "F") << " "
                 << (coordinates[i].selective_flags[2] ? "T" : "F");
        }
        file << "\n";
    }

    file.flush();
    if (!file) {
        std::cerr << "Error: failed writing to " << filenameOut << "\n";
        return false;
    }

    file.close();

    return true;
}

Atom POSCAR::displaceAtomWithVector(size_t atom_index, const DisplacementOptions& options) {
    Atom displacement{0.0, 0.0, 0.0};
    if (atom_index >= coordinates.size()) {
        return displacement;
    }

    displacement.x = randomDouble(-options.ampx, options.ampx);
    displacement.y = randomDouble(-options.ampy, options.ampy);
    displacement.z = randomDouble(-options.ampz, options.ampz);

    coordinates[atom_index].x += displacement.x;
    coordinates[atom_index].y += displacement.y;
    coordinates[atom_index].z += displacement.z;

    return displacement;
}

void POSCAR::displaceAtoms(int n_atoms, double amplitude) {
    DisplacementOptions options;
    options.ampx = amplitude;
    options.ampy = amplitude;
    options.ampz = amplitude;
    displaceAtoms(n_atoms, options);
}

void POSCAR::displaceAtoms(int n_atoms, const DisplacementOptions& options) {
    if (n_atoms <= 0 || total_atoms <= 0) {
        return;
    }

    std::vector<size_t> candidates;
    if (options.candidate_indices.empty()) {
        candidates.resize(static_cast<size_t>(total_atoms));
        for (int i = 0; i < total_atoms; ++i) {
            candidates[static_cast<size_t>(i)] = static_cast<size_t>(i);
        }
    } else {
        candidates = options.candidate_indices;
    }

    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(), [this](size_t idx) { return idx >= coordinates.size(); }),
        candidates.end());

    if (candidates.empty()) {
        return;
    }

    std::shuffle(candidates.begin(), candidates.end(), getGenerator());

    const size_t selected_count = std::min(static_cast<size_t>(n_atoms), candidates.size());

    bool was_direct = is_direct;
    if (was_direct) {
        toCartesian();
    }

    std::vector<size_t> selected_indices;
    selected_indices.reserve(selected_count);
    Atom net{0.0, 0.0, 0.0};
    for (size_t i = 0; i < selected_count; ++i) {
        const size_t idx = candidates[i];
        Atom d = displaceAtomWithVector(idx, options);
        net.x += d.x;
        net.y += d.y;
        net.z += d.z;
        selected_indices.push_back(idx);
    }

    if (options.zero_net && !selected_indices.empty()) {
        const double inv = 1.0 / static_cast<double>(selected_indices.size());
        const double mx = net.x * inv;
        const double my = net.y * inv;
        const double mz = net.z * inv;
        for (size_t idx : selected_indices) {
            coordinates[idx].x -= mx;
            coordinates[idx].y -= my;
            coordinates[idx].z -= mz;
        }
    }

    if (was_direct) {
        toDirect();
    }

    if (options.wrap_direct) {
        const bool started_direct = is_direct;
        if (!is_direct) {
            toDirect();
        }

        auto wrap01 = [](double v) {
            double wrapped = std::fmod(v, 1.0);
            if (wrapped < 0.0) {
                wrapped += 1.0;
            }
            return wrapped;
        };

        for (auto& atom : coordinates) {
            atom.x = wrap01(atom.x);
            atom.y = wrap01(atom.y);
            atom.z = wrap01(atom.z);
        }

        if (!started_direct) {
            toCartesian();
        }
    }
}

/*
OLD version
void POSCAR::displaceAtoms(int n_atoms, AmpMode amp_mode, double amplitude)
{
    // Generate random atom indices to displace
    std::vector<size_t> indices(total_atoms);
    for (int i = 0; i < total_atoms; ++i) indices[i] = i;

    // Shuffle the indices randomly
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.end(), gen);

    // Take the first n_atoms indices and displace those atoms
    for (int i = 0; i < n_atoms; ++i){
        if (amp_mode == AmpMode::Direct) {
            displaceAtomDirect(indices[i], amplitude);
        }
        else if (amp_mode == AmpMode::Cartesian) {
            displaceAtomCartesian(indices[i], amplitude);
        }
    }

}*/

void POSCAR::toDirect() {
    if (is_direct)
        return;

    double A[9];

    // Copy lattice
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A[i * 3 + j] = lattice[i][j];

    lapack_int ipiv[3];

    // LU factorization
    if (LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, A, 3, ipiv) != 0) {
        std::cerr << "Error: LU decomposition failed.\n";
        return;
    }

    // Compute inverse
    if (LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, A, 3, ipiv) != 0) {
        std::cerr << "Error: Matrix inversion failed.\n";
        return;
    }

    // Transform coordinates
    for (auto& atom : coordinates) {
        double x[3] = {atom.x, atom.y, atom.z};
        double y[3];

        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 0.0, y, 1);

        atom.x = y[0];
        atom.y = y[1];
        atom.z = y[2];
    }

    is_direct = true;
}

void POSCAR::toCartesian() {
    if (!is_direct)
        return;

    double A[9];

    // Flatten lattice (row-major)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A[i * 3 + j] = lattice[i][j];

    for (auto& atom : coordinates) {
        double x[3] = {atom.x, atom.y, atom.z};
        double y[3];

        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 0.0, y, 1);

        atom.x = y[0];
        atom.y = y[1];
        atom.z = y[2];
    }

    is_direct = false;
}

std::vector<std::string> POSCAR::atomSpecies() const {
    std::vector<std::string> species;
    species.reserve(static_cast<size_t>(total_atoms));

    for (size_t i = 0; i < elements.size() && i < num_atoms.size(); ++i) {
        for (int j = 0; j < num_atoms[i]; ++j) {
            species.push_back(elements[i]);
        }
    }

    if (species.size() > coordinates.size()) {
        species.resize(coordinates.size());
    }

    return species;
}

POSCAR POSCAR::makeSupercell(const double S[9]) {
    if (!is_direct) {
        toDirect();
    }

    POSCAR poscar_supercell;
    poscar_supercell.comment = comment + "_supercell";
    poscar_supercell.scale = scale;
    poscar_supercell.is_direct = true;
    poscar_supercell.elements = elements;

    int nx = static_cast<int>(S[0]);
    int ny = static_cast<int>(S[4]);
    int nz = static_cast<int>(S[8]);

    // Expanded lattice
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            poscar_supercell.lattice[i][j] = lattice[i][j] * S[i * 4];  // Diagonální expanze
        }
    }

    // Atoms number
    int multiplier = nx * ny * nz;
    for (int count : num_atoms) {
        poscar_supercell.num_atoms.push_back(count * multiplier);
        poscar_supercell.total_atoms += (count * multiplier);
    }

    // Atomic positions generation
    poscar_supercell.coordinates.reserve(poscar_supercell.total_atoms);
    int current_offset = 0;

    for (size_t e = 0; e < elements.size(); ++e) {
        for (int a = 0; a < num_atoms[e]; ++a) {
            Atom orig = coordinates[current_offset + a];

            // Trojitý cyklus pro replikaci v prostoru
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        Atom copy;
                        copy.x = (orig.x + i) / nx;
                        copy.y = (orig.y + j) / ny;
                        copy.z = (orig.z + k) / nz;
                        poscar_supercell.coordinates.push_back(copy);
                    }
                }
            }
        }
        current_offset += num_atoms[e];
    }

    return poscar_supercell;
}
