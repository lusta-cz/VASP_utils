#include <CLI/CLI.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error_parse.h"

struct MagComponents {
    double s{0.0}, p{0.0}, d{0.0}, total{0.0};
};

struct AtomMagnetic {
    int index{0};
    std::string element{""};
    std::vector<MagComponents> mag;  // mag[0] = x, mag[1] = y, mag[2] = z
};

bool checkNonCollinear(std::ifstream& file) {
    std::streampos original_pos = file.tellg();
    file.seekg(0, std::ios::beg);

    std::string line;
    bool is_ncl = false;
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

    file.clear();
    file.seekg(original_pos);
    return is_ncl;
}

// Scrapes VRHFIN lines and extracts counts from "ions per type ="
bool parseOutcarInosPerElement(std::ifstream& file, std::vector<int>& ions_per_element,
                               std::vector<std::string>& elements, ParseError& error) {
    file.seekg(0, std::ios::beg);
    std::string line;
    int line_number = 0;
    bool found_counts = false;

    while (std::getline(file, line)) {
        ++line_number;

        // Extract chemical symbol from POTCAR/VRHFIN info (e.g., "VRHFIN =Na: s1")
        if (line.find("VRHFIN") != std::string::npos) {
            size_t eq_pos = line.find("=");
            size_t col_pos = line.find(":");
            if (eq_pos != std::string::npos && col_pos != std::string::npos && col_pos > eq_pos) {
                std::string sym = line.substr(eq_pos + 1, col_pos - eq_pos - 1);
                // Clean whitespace padding
                sym.erase(0, sym.find_first_not_of(" \t"));
                sym.erase(sym.find_last_not_of(" \t") + 1);
                elements.push_back(sym);
            }
        }

        if (line.find("ions per type =") != std::string::npos) {
            size_t pos = line.find("=");
            std::stringstream ss(line.substr(pos + 1));
            int count;
            while (ss >> count) {
                ions_per_element.push_back(count);
            }
            found_counts = true;
            break;
        }
    }

    if (!found_counts) {
        return fail(error, ParseErrorKind::Parse, line_number, "Could not find 'ions per type =' line.");
    }

    if (elements.size() < ions_per_element.size()) {
        return fail(error, ParseErrorKind::Semantic, line_number,
                    "Mismatch: Found fewer VRHFIN element tags than types in 'ions per type'.");
    }

    // Shrink element strings to match exactly the number of atom types listed
    if (elements.size() > ions_per_element.size()) {
        elements.resize(ions_per_element.size());
    }

    return true;
}

bool parseOutcarMagneticData(std::ifstream& file, const std::vector<int>& ions_per_element,
                             const std::vector<std::string>& elements, std::vector<AtomMagnetic>& magnetic_data,
                             ParseError& error) {
    bool is_ncl = checkNonCollinear(file);
    int directions_to_read = is_ncl ? 3 : 1;
    std::vector<std::string> targets = {"magnetization (x)", "magnetization (y)", "magnetization (z)"};

    // Unroll ions_per_element and elements vectors into a sequential per-atom mapping helper
    // e.g., if elements = {"Na", "Cl"} and counts = {2, 3}, maps to: {"Na", "Na", "Cl", "Cl", "Cl"}
    std::vector<std::string> atom_to_element_map;
    for (size_t i = 0; i < ions_per_element.size(); ++i) {
        for (int c = 0; c < ions_per_element[i]; ++c) {
            if (i < elements.size()) {
                atom_to_element_map.push_back(elements[i]);
            }
        }
    }

    std::string line;

    for (int d = 0; d < directions_to_read; ++d) {
        std::streampos last_target_pos;
        bool found_target = false;
        int line_number = 0;
        int target_line_num = 0;

        file.clear();
        file.seekg(0, std::ios::beg);

        while (std::getline(file, line)) {
            ++line_number;
            if (line.find(targets[d]) != std::string::npos) {
                last_target_pos = file.tellg();
                target_line_num = line_number;
                found_target = true;
            }
        }

        if (!found_target) {
            if (d == 0) {
                return fail(error, ParseErrorKind::Semantic, line_number,
                            "No magnetization data found. Check if ISPIN=2 is enabled.");
            } else {
                return fail(error, ParseErrorKind::Parse, line_number, "Missing expected NCL block: " + targets[d]);
            }
        }

        file.clear();
        file.seekg(last_target_pos);
        line_number = target_line_num;
        for (int i = 0; i < 3; ++i) {
            if (!readRequiredLine(file, line, line_number, error, targets[d] + " header"))
                return false;
        }

        int current_atom_idx = 0;
        while (std::getline(file, line)) {
            ++line_number;

            // Check for separating horizontal line
            if (line.find("--------------------------------------------------") != std::string::npos) {
                // Read the row tracking final global totals
                if (std::getline(file, line)) {
                    ++line_number;
                    std::stringstream ss(line);
                    std::string label;
                    double s_val, p_val, d_val, tot_val;

                    if (ss >> label >> s_val >> p_val >> d_val >> tot_val) {
                        MagComponents total_comps{s_val, p_val, d_val, tot_val};

                        if (d == 0) {
                            AtomMagnetic total_row;
                            total_row.index = 0;
                            total_row.element = "total";
                            total_row.mag.push_back(total_comps);
                            magnetic_data.push_back(total_row);
                        } else {
                            if (!magnetic_data.empty()) {
                                magnetic_data.back().mag.push_back(total_comps);
                            }
                        }
                    }
                }
                break;
            }

            // Read regular lines containing atomic coordinates matrices
            std::stringstream ss(line);
            int ion;
            double s_val, p_val, d_val, tot_val;

            if (ss >> ion >> s_val >> p_val >> d_val >> tot_val) {
                MagComponents comps{s_val, p_val, d_val, tot_val};

                if (d == 0) {
                    AtomMagnetic atom;
                    atom.index = ion;

                    // Assign mapped element string using your dynamic configuration mapping
                    if (current_atom_idx < atom_to_element_map.size()) {
                        atom.element = atom_to_element_map[current_atom_idx];
                    }

                    atom.mag.push_back(comps);
                    magnetic_data.push_back(atom);
                } else {
                    if (current_atom_idx < magnetic_data.size()) {
                        magnetic_data[current_atom_idx].mag.push_back(comps);
                    }
                }
                current_atom_idx++;
            }
        }
    }
    return true;
}

bool parseFromOutcar(const std::string& filename, std::vector<int>& ions_per_element,
                     std::vector<std::string>& elements, std::vector<AtomMagnetic>& magnetic_data) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "I/O error reading " << filename << ": cannot open file\n";
        return false;
    }
    ParseError error;

    if (!parseOutcarInosPerElement(file, ions_per_element, elements, error)) {
        reportParseError(filename, error);
        return false;
    }
    if (!parseOutcarMagneticData(file, ions_per_element, elements, magnetic_data, error)) {
        reportParseError(filename, error);
        return false;
    }

    return true;
}

int main(int argc, char* argv[]) {
    CLI::App app{"Export of magnetization and structural information from OUTCAR"};

    std::string inputFile{"OUTCAR"};
    std::string outputFile{"magnetization_info.dat"};

    app.add_option("--input,-i", inputFile, "Input OUTCAR file")->capture_default_str()->check(CLI::ExistingFile);
    app.add_option("--output,-o", outputFile, "Output file for magnetization and structural data")
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    std::vector<int> ions_per_element;
    std::vector<std::string> elements;
    std::vector<AtomMagnetic> magnetic_data;

    if (!parseFromOutcar(inputFile, ions_per_element, elements, magnetic_data)) {
        return 1;
    }

    // Write matrix values out to disk
    std::ofstream out(outputFile);
    if (!out) {
        std::cerr << "Error writing to output file " << outputFile << "\n";
        return 1;
    }

    out << "# Ion\tElement\t";
    if (!magnetic_data.empty() && magnetic_data[0].mag.size() == 1) {
        out << "s_x\tp_x\td_x\ttot_x\n";
        for (const auto& atom : magnetic_data) {
            out << atom.index << "\t" << atom.element << "\t" << atom.mag[0].s << "\t" << atom.mag[0].p << "\t"
                << atom.mag[0].d << "\t" << atom.mag[0].total << "\n";
        }
    } else {
        out << "s_x\tp_x\td_x\ttot_x\ts_y\tp_y\td_y\ttot_y\ts_z\tp_z\td_z\ttot_z\n";
        for (const auto& atom : magnetic_data) {
            out << atom.index << "\t" << atom.element << "\t";
            for (size_t d = 0; d < atom.mag.size(); ++d) {
                out << atom.mag[d].s << "\t" << atom.mag[d].p << "\t" << atom.mag[d].d << "\t" << atom.mag[d].total
                    << "\t";
            }
            out << "\n";
        }
    }

    std::cout << "Successfully generated data file: " << outputFile << "\n";
    return 0;
}
