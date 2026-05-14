#include "surface.h"

#include <cblas.h>
#include <lapacke.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace {

using Vec3 = std::array<double, 3>;
using IVec3 = std::array<int, 3>;

// --- integer math ---

int gcd(int a, int b) {
    a = std::abs(a);
    b = std::abs(b);
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}
int gcd3(int a, int b, int c) {
    return gcd(gcd(a, b), c);
}

// --- 3D Cartesian helpers ---

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}
Vec3 normalized(const Vec3& v) {
    double n = norm(v);
    return {v[0] / n, v[1] / n, v[2] / n};
}

// Correct fractional → Cartesian (POSCAR row-vector convention):
// r[j] = sum_i f[i] * L[i][j]   (equivalent to L^T * f as column vector)
Vec3 fracToCart(const Vec3& f, const double L[3][3]) {
    Vec3 r{0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[j] += f[i] * L[i][j];
    return r;
}
Vec3 fracToCartI(const IVec3& f, const double L[3][3]) {
    Vec3 r{0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[j] += static_cast<double>(f[i]) * L[i][j];
    return r;
}

// Invert a 3×3 column-major matrix A[9] in place. Returns false on failure.
bool invert3x3Col(double A[9]) {
    lapack_int ipiv[3];
    if (LAPACKE_dgetrf(LAPACK_COL_MAJOR, 3, 3, A, 3, ipiv) != 0)
        return false;
    if (LAPACKE_dgetri(LAPACK_COL_MAJOR, 3, A, 3, ipiv) != 0)
        return false;
    return true;
}

// --- Reciprocal lattice (without 2π) ---
void computeReciprocal(const double L[3][3], double R[3][3]) {
    Vec3 a0{L[0][0], L[0][1], L[0][2]};
    Vec3 a1{L[1][0], L[1][1], L[1][2]};
    Vec3 a2{L[2][0], L[2][1], L[2][2]};
    Vec3 c12 = cross(a1, a2);
    Vec3 c20 = cross(a2, a0);
    Vec3 c01 = cross(a0, a1);
    double V = dot(a0, c12);
    for (int j = 0; j < 3; ++j) {
        R[0][j] = c12[j] / V;
        R[1][j] = c20[j] / V;
        R[2][j] = c01[j] / V;
    }
}

// --- In-plane vector search ---

// Find two shortest linearly-independent lattice vectors satisfying
// h*n0 + k*n1 + l*n2 = 0, via exhaustive search in [-M, M]^3.
void findInPlaneVectors(int h, int k, int l, const double L[3][3], IVec3& v1, IVec3& v2, int M) {
    double best1 = 1e30, best2 = 1e30;
    v1 = {0, 0, 0};
    v2 = {0, 0, 0};

    for (int n0 = -M; n0 <= M; ++n0)
        for (int n1 = -M; n1 <= M; ++n1)
            for (int n2 = -M; n2 <= M; ++n2) {
                if (h * n0 + k * n1 + l * n2 != 0)
                    continue;
                if (n0 == 0 && n1 == 0 && n2 == 0)
                    continue;
                IVec3 cand{n0, n1, n2};
                Vec3 c = fracToCartI(cand, L);
                double len2 = dot(c, c);
                if (len2 < best1 - 1e-12) {
                    Vec3 c1 = fracToCartI(v1, L);
                    if (dot(cross(c, c1), cross(c, c1)) > 1e-20) {
                        best2 = best1;
                        v2 = v1;
                    }
                    best1 = len2;
                    v1 = cand;
                } else if (len2 < best2 - 1e-12) {
                    Vec3 c1 = fracToCartI(v1, L);
                    if (dot(cross(c, c1), cross(c, c1)) > 1e-20) {
                        best2 = len2;
                        v2 = cand;
                    }
                }
            }
}

// Lagrange-Gauss 2D reduction of in-plane integer vectors.
void lagrangeGaussReduce(IVec3& v1, IVec3& v2, const double L[3][3]) {
    for (int iter = 0; iter < 200; ++iter) {
        Vec3 c1 = fracToCartI(v1, L);
        Vec3 c2 = fracToCartI(v2, L);
        if (dot(c2, c2) < dot(c1, c1))
            std::swap(v1, v2);
        c1 = fracToCartI(v1, L);
        c2 = fracToCartI(v2, L);
        double d11 = dot(c1, c1);
        if (d11 < 1e-30)
            break;
        int m = static_cast<int>(std::round(dot(c1, c2) / d11));
        if (m == 0)
            break;
        v2[0] -= m * v1[0];
        v2[1] -= m * v1[1];
        v2[2] -= m * v1[2];
    }
}

// Find the interlayer vector: lattice vector with h*n0+k*n1+l*n2 = g
// and smallest in-plane Cartesian component.
IVec3 findInterlayerVector(int h, int k, int l, int g, const Vec3& nhat, const double L[3][3], int M) {
    IVec3 best{0, 0, 0};
    double best_ip2 = 1e30;
    for (int n0 = -M; n0 <= M; ++n0)
        for (int n1 = -M; n1 <= M; ++n1)
            for (int n2 = -M; n2 <= M; ++n2) {
                if (h * n0 + k * n1 + l * n2 != g)
                    continue;
                IVec3 cand{n0, n1, n2};
                Vec3 c = fracToCartI(cand, L);
                double zp = dot(c, nhat);
                double ip2 = dot(c, c) - zp * zp;
                if (ip2 < best_ip2 - 1e-12) {
                    best_ip2 = ip2;
                    best = cand;
                }
            }
    return best;
}

double wrap01(double v) {
    v = v - std::floor(v);  // reliable for negative values
    if (v >= 1.0 - 1e-9)
        v = 0.0;  // snap near-1 to 0 (avoids boundary dedup failure)
    return v;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

POSCAR buildSlab(const POSCAR& bulk_in, const SlabOptions& opts) {
    if (opts.h == 0 && opts.k == 0 && opts.l == 0) {
        std::cerr << "Error: Miller indices must not all be zero.\n";
        return {};
    }
    if (opts.n_layers <= 0) {
        std::cerr << "Error: n_layers must be positive.\n";
        return {};
    }
    if (opts.n_frozen < 0 || opts.n_frozen > opts.n_layers) {
        std::cerr << "Error: n_frozen must be in [0, n_layers].\n";
        return {};
    }
    if (opts.vacuum < 0.0) {
        std::cerr << "Error: vacuum must be non-negative.\n";
        return {};
    }

    POSCAR bulk = bulk_in;
    bulk.toDirect();
    const double(&L)[3][3] = bulk.lattice;

    // --- surface normal and interplanar spacing ---
    double recip[3][3];
    computeReciprocal(L, recip);

    int h = opts.h, k = opts.k, l = opts.l;
    int g = gcd3(std::abs(h), std::abs(k), std::abs(l));

    Vec3 Gcart{0, 0, 0};
    for (int j = 0; j < 3; ++j)
        Gcart[j] = h * recip[0][j] + k * recip[1][j] + l * recip[2][j];
    double Gnorm = norm(Gcart);
    if (Gnorm < 1e-15) {
        std::cerr << "Error: degenerate surface normal.\n";
        return {};
    }
    Vec3 nhat = normalized(Gcart);
    double d_spacing = static_cast<double>(g) / Gnorm;  // crystallographic d-spacing

    // --- find in-plane and interlayer vectors ---
    int M = opts.n_layers + std::max({std::abs(h), std::abs(k), std::abs(l), 1}) + 3;

    IVec3 iv1, iv2;
    findInPlaneVectors(h, k, l, L, iv1, iv2, M);
    if (iv1 == IVec3{0, 0, 0} || iv2 == IVec3{0, 0, 0}) {
        std::cerr << "Error: could not find two linearly independent in-plane vectors.\n";
        return {};
    }
    lagrangeGaussReduce(iv1, iv2, L);

    IVec3 it = findInterlayerVector(h, k, l, g, nhat, L, M);
    if (it == IVec3{0, 0, 0}) {
        std::cerr << "Error: could not find interlayer vector.\n";
        return {};
    }

    Vec3 A = fracToCartI(iv1, L);
    Vec3 B = fracToCartI(iv2, L);

    // Ensure right-handed orientation.
    if (dot(cross(A, B), nhat) < 0) {
        std::swap(iv1, iv2);
        std::swap(A, B);
    }

    // --- info: check if surface cell area is larger than bulk projection ---
    {
        double surface_area = norm(cross(A, B));
        Vec3 a0{L[0][0], L[0][1], L[0][2]};
        Vec3 a1{L[1][0], L[1][1], L[1][2]};
        Vec3 a2{L[2][0], L[2][1], L[2][2]};
        double V_bulk = std::abs(dot(a0, cross(a1, a2)));
        double proj_area = V_bulk / d_spacing;
        int iratio = static_cast<int>(std::round(surface_area / proj_area));
        if (iratio > 1) {
            std::cout << "Info: surface (" << h << k << l << ") in-plane unit cell is " << iratio
                      << "x larger than the input cell projection.\n";
        }
    }

    // --- Build the in-plane decomposition matrix M_in (columns = A, B, nhat) ---
    // r_cart = c1*A + c2*B + z_p*nhat  →  [c1,c2,z_p] = M_in^{-1} * r_cart
    // Stored column-major for LAPACKE_COL_MAJOR inversion.
    double M_in[9] = {A[0], A[1], A[2], B[0], B[1], B[2], nhat[0], nhat[1], nhat[2]};
    if (!invert3x3Col(M_in)) {
        std::cerr << "Error: in-plane decomposition matrix is singular.\n";
        return {};
    }
    // M_in now holds the inverse (column-major): M_in_inv * r gives [c1, c2, z_p].

    auto decomposeCart = [&](const Vec3& r) -> std::array<double, 3> {
        // Multiply M_in_inv (col-major) by r:
        // result[i] = sum_j M_in_inv[i + 3*j] * r[j]  (col-major: A[i][j] = A[i + n*j])
        std::array<double, 3> out{0, 0, 0};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                out[i] += M_in[i + 3 * j] * r[j];
        return out;
    };

    // --- Enumerate atom images in a generous z range ---
    // We generate n_layers * d_spacing worth of atoms; each d_spacing has >= 1 atomic plane,
    // so this always provides >= n_layers planes. After clustering we trim to exactly n_layers.
    double z_search_max = (opts.n_layers + 0.5) * d_spacing;

    struct SlabAtomCart {
        double c1, c2;  // in-plane fractional coords (wrapped [0,1))
        double z_p;     // Cartesian z-projection along nhat
        std::string species;
    };
    std::vector<SlabAtomCart> candidates;

    std::vector<std::string> atom_species = bulk.atomSpecies();

    for (int n0 = -M; n0 <= M; ++n0)
        for (int n1 = -M; n1 <= M; ++n1)
            for (int n2 = -M; n2 <= M; ++n2)
                for (int ai = 0; ai < bulk.total_atoms; ++ai) {
                    Vec3 f_old{bulk.coordinates[static_cast<size_t>(ai)].x + n0,
                               bulk.coordinates[static_cast<size_t>(ai)].y + n1,
                               bulk.coordinates[static_cast<size_t>(ai)].z + n2};
                    Vec3 r = fracToCart(f_old, L);

                    double z_p = dot(r, nhat);
                    if (z_p < -1e-6 || z_p >= z_search_max - 1e-6)
                        continue;

                    auto decomp = decomposeCart(r);
                    double c1 = wrap01(decomp[0]);
                    double c2 = wrap01(decomp[1]);

                    candidates.push_back({c1, c2, z_p, atom_species[static_cast<size_t>(ai)]});
                }

    // --- deduplicate ---
    // Use periodic comparison for c1/c2 (values near 0 and near 1 are the same point).
    constexpr double tol = 1e-4;
    auto pdiff = [](double a, double b) {
        double d = std::abs(a - b);
        return d > 0.5 ? 1.0 - d : d;
    };
    std::vector<bool> keep(candidates.size(), true);
    for (size_t i = 0; i < candidates.size(); ++i) {
        if (!keep[i])
            continue;
        for (size_t j = i + 1; j < candidates.size(); ++j) {
            if (!keep[j])
                continue;
            const auto& a = candidates[i];
            const auto& b = candidates[j];
            if (pdiff(a.c1, b.c1) < tol && pdiff(a.c2, b.c2) < tol && std::abs(a.z_p - b.z_p) < tol)
                keep[j] = false;
        }
    }
    std::vector<SlabAtomCart> unique_atoms;
    unique_atoms.reserve(candidates.size());
    for (size_t i = 0; i < candidates.size(); ++i)
        if (keep[i])
            unique_atoms.push_back(candidates[i]);

    if (unique_atoms.empty()) {
        std::cerr << "Error: no atoms found — check Miller indices and input structure.\n";
        return {};
    }

    // --- sort by z, identify layers ---
    std::stable_sort(unique_atoms.begin(), unique_atoms.end(),
                     [](const SlabAtomCart& a, const SlabAtomCart& b) { return a.z_p < b.z_p; });

    // Cluster layers (tolerance in Angstroms).
    constexpr double layer_tol = 0.05;
    std::vector<int> layer_id(unique_atoms.size(), 0);
    int n_found = 1;
    for (size_t i = 1; i < unique_atoms.size(); ++i) {
        if (unique_atoms[i].z_p - unique_atoms[i - 1].z_p > layer_tol)
            ++n_found;
        layer_id[i] = n_found - 1;
    }

    if (n_found < opts.n_layers) {
        std::cerr << "Error: only " << n_found << " atomic planes found but " << opts.n_layers
                  << " requested. Try a larger Miller index search range.\n";
        return {};
    }

    // Trim to exactly n_layers: drop atoms with layer_id >= n_layers.
    {
        std::vector<SlabAtomCart> trimmed;
        trimmed.reserve(unique_atoms.size());
        for (size_t i = 0; i < unique_atoms.size(); ++i)
            if (layer_id[i] < opts.n_layers)
                trimmed.push_back(unique_atoms[i]);
        unique_atoms = std::move(trimmed);
    }
    // Recompute layer_id for trimmed set.
    layer_id.assign(unique_atoms.size(), 0);
    {
        int cur = 0;
        for (size_t i = 1; i < unique_atoms.size(); ++i) {
            if (unique_atoms[i].z_p - unique_atoms[i - 1].z_p > layer_tol)
                ++cur;
            layer_id[i] = cur;
        }
    }

    // --- determine final cell dimensions ---
    double z_top = unique_atoms.back().z_p;
    double C_mag = z_top + opts.vacuum;
    Vec3 C_vec{nhat[0] * C_mag, nhat[1] * C_mag, nhat[2] * C_mag};

    // --- group atoms by species (VASP requires species blocks) ---
    std::vector<std::string> ordered_species = bulk.elements;
    for (const auto& at : unique_atoms) {
        bool found = std::find(ordered_species.begin(), ordered_species.end(), at.species) != ordered_species.end();
        if (!found)
            ordered_species.push_back(at.species);
    }

    auto species_rank = [&](const std::string& s) {
        for (size_t i = 0; i < ordered_species.size(); ++i)
            if (ordered_species[i] == s)
                return static_cast<int>(i);
        return static_cast<int>(ordered_species.size());
    };

    // Sort primary: species order; secondary: z (z-order preserved within species).
    // Also need to preserve atom-to-layer mapping through the re-sort.
    std::vector<size_t> idx(unique_atoms.size());
    for (size_t i = 0; i < idx.size(); ++i)
        idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        int ra = species_rank(unique_atoms[a].species);
        int rb = species_rank(unique_atoms[b].species);
        if (ra != rb)
            return ra < rb;
        return unique_atoms[a].z_p < unique_atoms[b].z_p;
    });

    // --- build output POSCAR ---
    POSCAR slab;
    slab.comment = bulk.comment + "_surface_" + std::to_string(h) + std::to_string(k) + std::to_string(l);
    slab.scale = 1.0;
    slab.is_direct = true;
    slab.selective_dynamics = (opts.n_frozen > 0);

    // Lattice rows: A, B, C
    for (int j = 0; j < 3; ++j) {
        slab.lattice[0][j] = A[j];
        slab.lattice[1][j] = B[j];
        slab.lattice[2][j] = C_vec[j];
    }

    for (const auto& sp : ordered_species) {
        int cnt = 0;
        for (const auto& at : unique_atoms)
            if (at.species == sp)
                ++cnt;
        if (cnt > 0) {
            slab.elements.push_back(sp);
            slab.num_atoms.push_back(cnt);
        }
    }
    slab.total_atoms = static_cast<int>(unique_atoms.size());

    for (size_t ii = 0; ii < idx.size(); ++ii) {
        size_t i = idx[ii];
        const auto& at = unique_atoms[i];

        // Fractional coords in new cell: (c1, c2, z_p / C_mag).
        // Re-wrap c1/c2: dedup may have kept the 1.0 image instead of 0.0.
        Atom a;
        a.x = wrap01(at.c1);
        a.y = wrap01(at.c2);
        a.z = at.z_p / C_mag;

        bool frozen = (opts.n_frozen > 0) && (layer_id[i] < opts.n_frozen);
        a.selective_flags = {!frozen, !frozen, !frozen};
        slab.coordinates.push_back(a);
    }

    return slab;
}
