#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "poscar_file.h"
#include "surface.h"

static const std::string kNaClPath = std::string(TEST_DATA_DIR) + "/NaCl_conv_fcc.poscar";

class SurfaceTest : public ::testing::Test {
protected:
    POSCAR nacl;
    void SetUp() override {
        ASSERT_TRUE(nacl.readPOSCAR(kNaClPath));
    }
};

// Helper: count distinct z-levels (atomic layers) in a slab POSCAR.
// Uses f[2] * |C| to get physical z — avoids the toCartesian() L*f vs L^T*f issue
// for non-diagonal cells.
static int countLayers(const POSCAR& slab, double tol_angstrom = 0.05) {
    // |C| = magnitude of third lattice vector
    double cx = slab.lattice[2][0], cy = slab.lattice[2][1], cz = slab.lattice[2][2];
    double C_mag = std::sqrt(cx * cx + cy * cy + cz * cz);

    std::vector<double> zvals;
    zvals.reserve(slab.coordinates.size());
    for (const auto& at : slab.coordinates)
        zvals.push_back(at.z * C_mag);  // fractional z × |C| = physical height (Å)

    std::sort(zvals.begin(), zvals.end());
    if (zvals.empty())
        return 0;
    int layers = 1;
    for (size_t i = 1; i < zvals.size(); ++i)
        if (zvals[i] - zvals[i - 1] > tol_angstrom)
            ++layers;
    return layers;
}

// (001) of cubic NaCl: trivial surface, c-axis already the normal.
TEST_F(SurfaceTest, NaCl001_LayerCount) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 4;
    opts.n_frozen = 0;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);
    EXPECT_EQ(countLayers(slab), opts.n_layers);
}

// Atom count for (001) slab: NaCl cubic has 2 atoms per (001) layer in the conv. cell.
TEST_F(SurfaceTest, NaCl001_AtomCount) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 6;
    opts.n_frozen = 0;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);
    // NaCl conventional (001): 4 atoms per layer in the a×a surface cell
    // (2 Na + 2 Cl per layer, alternating checkerboard pattern).
    EXPECT_EQ(slab.total_atoms, 4 * opts.n_layers);
}

// Selective dynamics: frozen layer flags.
TEST_F(SurfaceTest, NaCl001_SelectiveDynamics) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 6;
    opts.n_frozen = 2;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);
    EXPECT_TRUE(slab.selective_dynamics);

    // Bottom 2 layers should be frozen (F F F), rest free (T T T).
    // Build (z_physical, flags) pairs and sort by z.
    double cx = slab.lattice[2][0], cy = slab.lattice[2][1], cz2 = slab.lattice[2][2];
    double C_mag = std::sqrt(cx * cx + cy * cy + cz2 * cz2);
    std::vector<std::pair<double, std::array<bool, 3>>> zflags;
    for (size_t i = 0; i < slab.coordinates.size(); ++i)
        zflags.push_back({slab.coordinates[i].z * C_mag, slab.coordinates[i].selective_flags});
    std::sort(zflags.begin(), zflags.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    // Layer clustering
    std::vector<int> layer_of(zflags.size(), 0);
    for (size_t i = 1; i < zflags.size(); ++i) {
        layer_of[i] = layer_of[i - 1];
        if (zflags[i].first - zflags[i - 1].first > 0.05)
            ++layer_of[i];
    }

    int n_frozen = opts.n_frozen;
    for (size_t i = 0; i < zflags.size(); ++i) {
        bool should_be_frozen = (layer_of[i] < n_frozen);
        const auto& flags = zflags[i].second;
        if (should_be_frozen) {
            EXPECT_FALSE(flags[0]) << "atom " << i << " layer " << layer_of[i] << " should be frozen";
            EXPECT_FALSE(flags[1]) << "atom " << i << " layer " << layer_of[i] << " should be frozen";
            EXPECT_FALSE(flags[2]) << "atom " << i << " layer " << layer_of[i] << " should be frozen";
        } else {
            EXPECT_TRUE(flags[0]) << "atom " << i << " layer " << layer_of[i] << " should be free";
            EXPECT_TRUE(flags[1]) << "atom " << i << " layer " << layer_of[i] << " should be free";
            EXPECT_TRUE(flags[2]) << "atom " << i << " layer " << layer_of[i] << " should be free";
        }
    }
}

// Vacuum: slab c-vector length >= vacuum.
TEST_F(SurfaceTest, NaCl001_VacuumThickness) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 4;
    opts.n_frozen = 0;
    opts.vacuum = 20.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);

    // c-vector magnitude
    double cx = slab.lattice[2][0], cy = slab.lattice[2][1], cz = slab.lattice[2][2];
    double c_mag = std::sqrt(cx * cx + cy * cy + cz * cz);
    EXPECT_GE(c_mag, opts.vacuum - 0.1);
}

// (111) surface: FCC NaCl — 3-layer repeat, 2 atoms per layer per conv. cell
TEST_F(SurfaceTest, NaCl111_LayerCount) {
    SlabOptions opts;
    opts.h = 1;
    opts.k = 1;
    opts.l = 1;
    opts.n_layers = 6;
    opts.n_frozen = 0;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);
    EXPECT_EQ(countLayers(slab), opts.n_layers);
}

// No frozen layers → selective_dynamics should be false.
TEST_F(SurfaceTest, NoFrozenLayers_NoSelectiveDynamics) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 4;
    opts.n_frozen = 0;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    ASSERT_GT(slab.total_atoms, 0);
    EXPECT_FALSE(slab.selective_dynamics);
}

// Invalid: all-zero Miller indices.
TEST_F(SurfaceTest, InvalidMillerZero) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 0;
    opts.n_layers = 4;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    EXPECT_EQ(slab.total_atoms, 0);
}

// Invalid: frozen > layers.
TEST_F(SurfaceTest, InvalidFrozenExceedsLayers) {
    SlabOptions opts;
    opts.h = 0;
    opts.k = 0;
    opts.l = 1;
    opts.n_layers = 4;
    opts.n_frozen = 6;
    opts.vacuum = 15.0;

    POSCAR slab = buildSlab(nacl, opts);
    EXPECT_EQ(slab.total_atoms, 0);
}
