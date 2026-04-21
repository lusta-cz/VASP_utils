#include <CLI/CLI.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

struct KPoint {
    double x, y, z;
    std::vector<double> energies_up;
    std::vector<double> energies_dn;  // empty if nspin == 1
};

struct EigenvalData {
    int total_kpoints{0};
    int nbands{0};
    int nspin{1};
    std::vector<KPoint> kpoints;
};

namespace {

enum class ParseErrorKind { Io, Parse, Semantic };

struct ParseError {
    ParseErrorKind kind{ParseErrorKind::Parse};
    int line_number{0};
    std::string message;
};

bool parseError(ParseError& error, ParseErrorKind kind, int line_number, std::string message) {
    error.kind = kind;
    error.line_number = line_number;
    error.message = std::move(message);
    return false;
}

void reportParseError(const std::string& filename, const ParseError& error) {
    const char* category = (error.kind == ParseErrorKind::Io) ? "I/O error" : "Parse error";
    std::cerr << category << " reading from " << filename;
    if (error.line_number > 0)
        std::cerr << " at line " << error.line_number;
    std::cerr << ": " << error.message << "\n";
}

bool readRequiredLine(std::ifstream& file, std::string& line, int& line_number, ParseError& error,
                      const std::string& context) {
    if (!std::getline(file, line)) {
        return parseError(error, file.bad() ? ParseErrorKind::Io : ParseErrorKind::Parse, line_number + 1,
                          "unexpected end of file while reading " + context);
    }
    ++line_number;
    return true;
}

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
        return parseError(error, ParseErrorKind::Parse, line_number,
                          "expected 4 values (Emax Emin NEDOS Efermi), got: " + line);

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
        return parseError(error, ParseErrorKind::Parse, 0, "'E-fermi' entry not found");

    return true;
}

bool parseKpointsFile(std::ifstream& file, int& kpts_per_seg, ParseError& error) {
    std::string line;
    int line_number = 0;

    // Line 1: comment/title — skip
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS comment line"))
        return false;

    // Line 2: number of k-points per segment (first field)
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS k-points per segment"))
        return false;

    std::istringstream ss(line);
    if (!(ss >> kpts_per_seg))
        return parseError(error, ParseErrorKind::Parse, line_number,
                          "expected integer k-points per segment, got: " + line);

    if (kpts_per_seg <= 1)
        return parseError(error, ParseErrorKind::Parse, line_number, "k-points per segment must be greater than 1");

    // Line 3: must start with 'L' or 'l' to confirm line mode
    if (!readRequiredLine(file, line, line_number, error, "KPOINTS mode line"))
        return false;

    if (line.empty() || std::toupper(static_cast<unsigned char>(line[0])) != 'L')
        return parseError(error, ParseErrorKind::Semantic, line_number,
                          "KPOINTS file is not in line mode (line 3 must start with 'L' or 'l'), got: " + line);

    return true;
}

bool parseEigenvalFile(std::ifstream& file, EigenvalData& data, ParseError& error) {
    std::string line;
    int line_number = 0;

    for (int i = 0; i < 5; ++i) {
        if (!readRequiredLine(file, line, line_number, error, "EIGENVAL header"))
            return false;
    }

    if (!readRequiredLine(file, line, line_number, error, "EIGENVAL line 6 (NIONS NKPTS NBANDS)"))
        return false;

    {
        std::istringstream ss(line);
        int nions;
        if (!(ss >> nions >> data.total_kpoints >> data.nbands))
            return parseError(error, ParseErrorKind::Parse, line_number, "expected NIONS NKPTS NBANDS, got: " + line);
    }

    if (data.total_kpoints <= 0)
        return parseError(error, ParseErrorKind::Semantic, line_number, "total k-points must be positive");
    if (data.nbands <= 0)
        return parseError(error, ParseErrorKind::Semantic, line_number, "number of bands must be positive");

    data.kpoints.reserve(data.total_kpoints);
    bool nspin_detected = false;

    for (int ikpt = 0; ikpt < data.total_kpoints; ++ikpt) {
        if (!readRequiredLine(file, line, line_number, error, "k-point separator"))
            return false;

        if (!readRequiredLine(file, line, line_number, error, "k-point coordinates"))
            return false;

        KPoint kp;
        {
            std::istringstream ss(line);
            double weight;
            if (!(ss >> kp.x >> kp.y >> kp.z >> weight))
                return parseError(error, ParseErrorKind::Parse, line_number,
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
                return parseError(error, ParseErrorKind::Parse, line_number,
                                  "expected band index and eigenvalue at k-point " + std::to_string(ikpt + 1) +
                                      ", band " + std::to_string(ib + 1));

            kp.energies_up.push_back(e_up);

            // Detect nspin from first band line of first k-point
            if (!nspin_detected && ikpt == 0 && ib == 0) {
                double e_dn, occ;
                if (ss >> e_dn >> occ)
                    data.nspin = 2;
                else
                    data.nspin = 1;
                nspin_detected = true;
            }

            if (data.nspin == 2) {
                std::istringstream ss2(line);
                int tmp;
                double eu, ed;
                if (!(ss2 >> tmp >> eu >> ed))
                    return parseError(error, ParseErrorKind::Parse, line_number,
                                      "expected spin-down eigenvalue at k-point " + std::to_string(ikpt + 1) +
                                          ", band " + std::to_string(ib + 1));
                kp.energies_dn.push_back(ed);
            }
        }

        data.kpoints.push_back(std::move(kp));
    }

    return true;
}

}  // namespace

bool parseFromDoscar(const std::string& filename, double& e_fermi) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading Fermi level from " << filename << ": cannot open file\n";
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
        std::cerr << "I/O error reading Fermi level from " << filename << ": cannot open file\n";
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

int main(int argc, char** argv) {
    CLI::App app{"VASP EIGENVAL to Band Structure Utility"};

    // Configuration and default values
    std::string eigen_file{"EIGENVAL"};
    std::string kpoints_file{"KPOINTS"};
    std::string doscar_file{"DOSCAR"};
    std::string outcar_file{"OUTCAR"};
    double manual_fermi{0.0};
    std::string fermiSource{"doscar"};

    // CLI Options
    app.add_option("-i,--input", eigen_file, "Input EIGENVAL file")->default_val("EIGENVAL")->check(CLI::ExistingFile);
    app.add_option("-k,--kpoints", kpoints_file, "KPOINTS file used for the band structure calculation")
        ->default_val("KPOINTS")
        ->check(CLI::ExistingFile);
    app.add_option("--fermi-source", fermiSource, "Fermi level source: doscar (default) | outcar | manual | none")
        ->capture_default_str()
        ->check(CLI::IsMember({"doscar", "outcar", "manual", "none"}));
    app.add_option("--doscar", doscar_file, "DOSCAR file path (used with --fermi-source=doscar)")
        ->capture_default_str();
    app.add_option("--outcar", outcar_file, "OUTCAR file path (used with --fermi-source=outcar)")
        ->capture_default_str();
    auto* fermi_opt = app.add_option("-e,--fermi", manual_fermi, "Fermi level in eV (used with --fermi-source=manual)")
                          ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    double e_fermi{0.0};
    bool success{true};

    // Logic to determine the active Fermi source
    if (fermiSource == "doscar") {
        if (!parseFromDoscar(doscar_file, e_fermi)) {
            success = false;
        } else {
            std::cout << "Using Fermi level from " << doscar_file << ": " << e_fermi << " eV\n";
        }
    } else if (fermiSource == "outcar") {
        if (!parseFromOutcar(outcar_file, e_fermi)) {
            success = false;
        } else {
            std::cout << "Using Fermi level from " << outcar_file << ": " << e_fermi << " eV\n";
        }
    } else if (fermiSource == "manual") {
        if (fermi_opt->count() == 0) {
            std::cerr << "Error: --fermi-source=manual requires -e/--fermi VALUE\n";
            return 1;
        }
        e_fermi = manual_fermi;
        std::cout << "Info: Fermi energy is set to: " << e_fermi << " eV\n";
    } else {
        std::cout << "Warning: Energy shift is disabled. 0.0 eV is NOT the Fermi level!\n";
    }

    if (!success) {
        std::cout << "Warning: Reading of the Fermi level failed! Continuing without energy shift!\n";
    }

    int kpoints_between{0};
    // --- Parsing from KPOINTS file
    if (!parseKpoints(kpoints_file, kpoints_between)) {
        return 1;
    }

    EigenvalData data;
    // --- Parsing from EIGENVAL file
    if (!parseEigenval(eigen_file, data)) {
        return 1;
    }

    // --- Path Calculation and Output Generation ---
    if (data.total_kpoints == 0 || kpoints_between == 0) {
        std::cerr << "Error: Something went wrong while reading the number of k points!\n";
        return 1;
    }

    int num_paths = data.total_kpoints / kpoints_between;
    double cumulative_dist = 0.0;

    for (int p = 0; p < num_paths; ++p) {
        std::string out_name = "band_path_" + std::to_string(p + 1) + ".dat";
        std::ofstream out(out_name);
        out << std::fixed << std::setprecision(10);

        int s_idx = p * kpoints_between;
        int e_idx = s_idx + kpoints_between - 1;

        double dx = data.kpoints[e_idx].x - data.kpoints[s_idx].x;
        double dy = data.kpoints[e_idx].y - data.kpoints[s_idx].y;
        double dz = data.kpoints[e_idx].z - data.kpoints[s_idx].z;
        double segment_dist = std::sqrt(dx * dx + dy * dy + dz * dz);
        double step = (kpoints_between > 1) ? (segment_dist / (kpoints_between - 1)) : 0.0;

        for (int b = 0; b < data.nbands; ++b) {
            for (int k = 0; k < kpoints_between; ++k) {
                double dist = cumulative_dist + k * step;
                const auto& kp = data.kpoints[s_idx + k];
                if (data.nspin == 2)
                    out << dist << " " << kp.energies_up[b] - e_fermi << " " << kp.energies_dn[b] - e_fermi << "\n";
                else
                    out << dist << " " << kp.energies_up[b] - e_fermi << "\n";
            }
            out << "\n";
        }

        cumulative_dist += segment_dist;
    }

    std::cout << "Success: Generated " << num_paths << " data files for plotting.\n";

    return 0;
}
