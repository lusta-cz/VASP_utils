#include "cif_file.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// gemmi: header-only CIF parsing + crystallography
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>
#include <gemmi/symmetry.hpp>

#include "poscar_file.h"

namespace {

// ── Lattice ↔ cell-parameter conversions ─────────────────────────────────────

static void lattice_to_cell(const double lat[3][3], double& a, double& b, double& c, double& alpha_deg,
                            double& beta_deg, double& gamma_deg) {
    auto dot3 = [](const double* u, const double* v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
    auto safe_acos = [](double x) { return std::acos(std::max(-1.0, std::min(1.0, x))) * (180.0 / M_PI); };
    a = std::sqrt(dot3(lat[0], lat[0]));
    b = std::sqrt(dot3(lat[1], lat[1]));
    c = std::sqrt(dot3(lat[2], lat[2]));
    alpha_deg = safe_acos(dot3(lat[1], lat[2]) / (b * c));  // b·c
    beta_deg = safe_acos(dot3(lat[0], lat[2]) / (a * c));   // a·c
    gamma_deg = safe_acos(dot3(lat[0], lat[1]) / (a * b));  // a·b
}

// Standard crystallographic convention: a ∥ x, b in xy-plane, c general.
static void cell_to_lattice(double a, double b, double c, double alpha_d, double beta_d, double gamma_d,
                            double lat[3][3]) {
    const double deg2rad = M_PI / 180.0;
    double ca = std::cos(alpha_d * deg2rad);
    double cb = std::cos(beta_d * deg2rad);
    double cg = std::cos(gamma_d * deg2rad);
    double sg = std::sin(gamma_d * deg2rad);

    lat[0][0] = a;
    lat[0][1] = 0.0;
    lat[0][2] = 0.0;
    lat[1][0] = b * cg;
    lat[1][1] = b * sg;
    lat[1][2] = 0.0;
    double cx = c * cb;
    double cy = c * (ca - cb * cg) / sg;
    double cz = std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy));
    lat[2][0] = cx;
    lat[2][1] = cy;
    lat[2][2] = cz;
}

// ── Element symbol extraction ─────────────────────────────────────────────────

// Strip leading alphabetic characters from a string, capitalise first letter.
static std::string alphaPrefix(const std::string& s) {
    std::string r;
    for (char ch : s) {
        if (!std::isalpha(static_cast<unsigned char>(ch)))
            break;
        r += ch;
    }
    if (r.empty())
        return r;
    r[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(r[0])));
    for (std::size_t i = 1; i < r.size(); ++i)
        r[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(r[i])));
    return r;
}

static std::string elementFromSite(const gemmi::SmallStructure::Site& site) {
    // 1. gemmi-parsed element field ("X" means unknown)
    if (site.element.elem != gemmi::El::X)
        return site.element.name();

    // 2. type_symbol from CIF (may carry charge suffix like "Fe2+")
    if (!site.type_symbol.empty()) {
        std::string s = alphaPrefix(site.type_symbol);
        if (!s.empty())
            return s;
    }

    // 3. Atom label (e.g. "Fe1", "O_2")
    return alphaPrefix(site.label);
}

}  // namespace

// ── Public API ────────────────────────────────────────────────────────────────

bool readCif(const std::string& filename, POSCAR& poscar, double occ_threshold) {
    // ── Parse CIF document ────────────────────────────────────────────────────
    gemmi::cif::Document doc;
    try {
        doc = gemmi::cif::read_file(filename);
    } catch (const std::exception& e) {
        std::cerr << "Error reading CIF " << filename << ": " << e.what() << "\n";
        return false;
    }

    if (doc.blocks.empty()) {
        std::cerr << "Error: no data blocks in " << filename << "\n";
        return false;
    }
    if (doc.blocks.size() > 1)
        std::cout << "Info: CIF contains " << doc.blocks.size() << " data blocks; using first: " << doc.blocks[0].name
                  << "\n";

    // ── Build SmallStructure (cell + asymmetric unit + space group) ────────────
    gemmi::SmallStructure st;
    try {
        st = gemmi::make_small_structure_from_block(doc.blocks[0]);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing structure from CIF block: " << e.what() << "\n";
        return false;
    }

    if (st.sites.empty()) {
        std::cerr << "Error: no atom sites found in CIF\n";
        return false;
    }

    // ── Expand symmetry to P1 manually via GroupOps ───────────────────────────
    // (avoids relying on SmallStructure::expand_to_p1 which was added later)
    const gemmi::SpaceGroup* sg = nullptr;
    if (!st.spacegroup_hm.empty())
        sg = gemmi::find_spacegroup_by_name(st.spacegroup_hm);

    // Build flat list of all operations = sym_ops × cen_ops
    struct FullOp {
        gemmi::Op sym;
        std::array<int, 3> cen;
    };
    std::vector<FullOp> all_ops;
    if (sg) {
        gemmi::GroupOps gops = sg->operations();
        for (const auto& s : gops.sym_ops)
            for (const auto& c : gops.cen_ops)
                all_ops.push_back({s, c});
    } else {
        all_ops.push_back({gemmi::Op::identity(), {0, 0, 0}});
    }

    // ── Collect atoms, apply all ops, deduplicate ─────────────────────────────
    struct AtomEntry {
        std::string elem;
        Atom atom;
    };
    std::vector<AtomEntry> all_atoms;
    std::vector<std::string> elem_order;
    int skipped = 0;
    const double tol = 1e-4;

    for (const auto& site : st.sites) {
        if (site.occ < occ_threshold) {
            ++skipped;
            continue;
        }

        std::string elem = elementFromSite(site);
        if (elem.empty()) {
            std::cerr << "Warning: cannot determine element for site '" << site.label << "'; skipping\n";
            continue;
        }

        for (const auto& fop : all_ops) {
            // Apply symmetry operation
            auto xyz = fop.sym.apply_to_xyz({site.fract.x, site.fract.y, site.fract.z});
            // Add centering translation
            for (int i = 0; i < 3; ++i)
                xyz[i] += static_cast<double>(fop.cen[i]) / gemmi::Op::DEN;
            // Wrap to [0, 1)
            for (auto& v : xyz) {
                v -= std::floor(v);
                if (v >= 1.0 - tol)
                    v = 0.0;
            }

            // Deduplicate: check against already-added atoms of the same element
            bool dup = false;
            for (const auto& entry : all_atoms) {
                if (entry.elem != elem)
                    continue;
                double dx = xyz[0] - entry.atom.x;
                double dy = xyz[1] - entry.atom.y;
                double dz = xyz[2] - entry.atom.z;
                // Periodic boundary correction
                dx -= std::round(dx);
                dy -= std::round(dy);
                dz -= std::round(dz);
                if (std::abs(dx) < tol && std::abs(dy) < tol && std::abs(dz) < tol) {
                    dup = true;
                    break;
                }
            }
            if (dup)
                continue;

            if (std::find(elem_order.begin(), elem_order.end(), elem) == elem_order.end())
                elem_order.push_back(elem);

            Atom a;
            a.x = xyz[0];
            a.y = xyz[1];
            a.z = xyz[2];
            all_atoms.push_back({elem, a});
        }
    }

    if (skipped > 0)
        std::cout << "Info: skipped " << skipped << " site(s) with occupancy < " << occ_threshold << "\n";

    if (all_atoms.empty()) {
        std::cerr << "Error: no usable atom sites after filtering\n";
        return false;
    }

    // ── Populate POSCAR ───────────────────────────────────────────────────────
    poscar = POSCAR{};
    poscar.comment = doc.blocks[0].name;
    poscar.scale = 1.0;
    poscar.is_direct = true;
    cell_to_lattice(st.cell.a, st.cell.b, st.cell.c, st.cell.alpha, st.cell.beta, st.cell.gamma, poscar.lattice);
    poscar.elements = elem_order;
    poscar.total_atoms = 0;

    for (const auto& elem : elem_order) {
        int cnt = 0;
        for (const auto& entry : all_atoms) {
            if (entry.elem == elem) {
                poscar.coordinates.push_back(entry.atom);
                ++cnt;
            }
        }
        poscar.num_atoms.push_back(cnt);
        poscar.total_atoms += cnt;
    }

    if (poscar.total_atoms > static_cast<int>(st.sites.size()))
        std::cout << "Symmetry expansion (" << st.spacegroup_hm << "): " << st.sites.size() << " asymmetric site(s) → "
                  << poscar.total_atoms << " atom(s) in P1\n";

    return true;
}

bool writeCif(const std::string& filename, const POSCAR& poscar) {
    // Work on a copy so we can safely convert to fractional if needed.
    POSCAR tmp = poscar;
    if (!tmp.is_direct)
        tmp.toDirect();

    double a, b, c, alpha, beta, gamma;
    lattice_to_cell(tmp.lattice, a, b, c, alpha, beta, gamma);

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: cannot write " << filename << "\n";
        return false;
    }

    // Sanitise block name (CIF identifiers: alphanumeric + '_')
    std::string block_name = tmp.comment;
    for (char& ch : block_name)
        if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_')
            ch = '_';
    if (block_name.empty())
        block_name = "structure";

    out << std::fixed;
    out << "data_" << block_name << "\n\n";

    // Cell parameters
    out << std::setprecision(6) << "_cell_length_a     " << a << "\n"
        << "_cell_length_b     " << b << "\n"
        << "_cell_length_c     " << c << "\n"
        << "_cell_angle_alpha  " << alpha << "\n"
        << "_cell_angle_beta   " << beta << "\n"
        << "_cell_angle_gamma  " << gamma << "\n\n";

    // Space group — always P1 (all atoms explicit)
    out << "_symmetry_space_group_name_H-M  'P 1'\n"
        << "_symmetry_Int_Tables_number      1\n\n";

    // Atom site loop
    out << "loop_\n"
        << "_atom_site_label\n"
        << "_atom_site_type_symbol\n"
        << "_atom_site_fract_x\n"
        << "_atom_site_fract_y\n"
        << "_atom_site_fract_z\n"
        << "_atom_site_occupancy\n";

    out << std::setprecision(8);

    // Per-element label counter (Fe1, Fe2, O1, O2, …)
    std::vector<int> label_cnt(tmp.elements.size(), 0);
    int coord_idx = 0;
    for (std::size_t ei = 0; ei < tmp.elements.size(); ++ei) {
        const std::string& elem = tmp.elements[ei];
        for (int ai = 0; ai < tmp.num_atoms[ei]; ++ai, ++coord_idx) {
            const Atom& at = tmp.coordinates[coord_idx];
            std::string label = elem + std::to_string(++label_cnt[ei]);
            out << std::left << std::setw(8) << label << std::left << std::setw(4) << elem << std::right
                << std::setw(14) << at.x << std::right << std::setw(14) << at.y << std::right << std::setw(14) << at.z
                << "  1.0\n";
        }
    }

    return out.good();
}
