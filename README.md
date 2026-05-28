Still under development. Currently, these utilities for POSCAR file manipulation are implemented:

----Pre-processing utilities:----

- poscar_d2c - fractional coordinates to cartesian
- poscar_c2d - cartesian coordinates to fractional
- poscar_symmetry - find symmetry of a cell
- poscar_supercell - create supercell
- poscar2primitive - create primitive cell
- poscar2conventional - create conventional cell
- poscar2ctrls - create ctrls file for ecalj/Questaal package from POSCAR
- poscar_atom_displace - randomly displace atoms
- poscar_convergence_test - automatically builds structured directories for systematically testing plane-wave cutoffs and k-point spacings
- poscar_kpath - for automatic generation of path in the Brillouin zone for band structure calculation
- poscar_surface - generate surface slab from a bulk structure
- poscar_elastic_deform - generate deformed structures to calculate elastic constants (both energy-strain and stress-strain approaches)
- poscar2cif - convert POSCAR to CIF format
- cif2poscar - convert CIF to POSCAR format (symmetry expansion via gemmi)


----Post-processing utilities:----

- vasp_eigenval2bands - create data file for each path in 1BZ for simpler plotting
- vasp_elastic_fit - calculating elastic tensor from the VASP calculation
- vasp_thermo - calculating mechanical, acustic and thermal properties from elastic tensor (vasp_elastic_fit)
- vasp_doscar2dos - reorganise DOSCAR data into column files for easy plotting, with E_Fermi shift
- vasp_convergence_test - gathers final free energies and automatic k-point dimensions from test folders into ready-to-plot data files
- vasp_magnetic - parses element-wise magnetic components and total magnetic moments from calculation outputs


For now, the code is as it is; nothing is guaranteed.

----Plans:----


----Installation:----

Requires CMake >= 3.22, a C++17 compiler, and BLAS/LAPACK/LAPACKE installed.
CLI11 and spglib are fetched automatically by CMake.

```
mkdir build && cd build
cmake ..
cmake --build . --parallel
```

Executables are placed in `build/bin/`.

----Usage:----

The following tools use CLI11 for argument parsing (run with `--help` for full options):

**poscar_d2c** — Direct to Cartesian conversion
```
poscar_d2c [--input/-i <file>] [--output/-o <file>]
```
Defaults: input `POSCAR`, output `POSCAR_cartesian`.

**poscar_c2d** — Cartesian to Direct conversion
```
poscar_c2d [--input/-i <file>] [--output/-o <file>]
```
Defaults: input `POSCAR`, output `POSCAR_direct`.


**poscar_supercell** — Build a diagonal supercell
```
poscar_supercell [--input/-i <file>] [--output/-o <file>] [--dim <nx> <ny> <nz>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `POSCAR_nx_ny_nz` | Output POSCAR file |
| `--dim` | `1 1 1` | Supercell dimensions along each lattice vector |

Example:
```
poscar_supercell --dim 2 2 3
```

**poscar2primitive** — Convert to primitive cell using spglib
```
poscar2primitive [--input/-i <file>] [--output/-o <file>] [--symprec <float>] [--idealize]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `POSCAR_primitive` | Output POSCAR file |
| `--symprec` | `1e-5` | Symmetry tolerance |
| `--idealize` | off | Idealize lattice parameters of the primitive cell |

**poscar2conventional** — Convert to conventional cell using spglib
```
poscar2conventional [--input/-i <file>] [--output/-o <file>] [--symprec <float>] [--idealize]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `POSCAR_conventional` | Output POSCAR file |
| `--symprec` | `1e-5` | Symmetry tolerance |
| `--idealize` | off | Idealize lattice parameters of the conventional cell |

**poscar2ctrls** — Convert POSCAR to ctrls.in format for ecalj/Questaal
```
poscar2ctrls [--input/-i <file>] [--output/-o <file>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `ctrls.in` | Output ctrls file |

**poscar_atom_displace** — Randomly displace atoms to generate perturbed structures
```
poscar_atom_displace [--input/-i <file>] [--nfiles/-n <int>] [--natoms <int>] [--allatoms]
                     [--amp/-a <float>] [--ampx <float>] [--ampy <float>] [--ampz <float>]
                     [--seed <int>] [--species <Na,Cl,...>] [--indices <1,3-5,...>]
                     [--wrap] [--zero-net]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--nfiles/-n` | `1` | Number of output files |
| `--natoms` | `1` | Number of atoms to displace |
| `--allatoms` | off | Displace all eligible atoms (overrides `--natoms`) |
| `--amp/-a` | `0.01` | Isotropic max displacement amplitude in Angstroms |
| `--ampx` / `--ampy` / `--ampz` | `0.01` | Axis-specific max displacement amplitudes (override `--amp`) |
| `--seed` | random | Seed RNG for deterministic perturbations |
| `--species` | all species | Comma-separated element symbols to consider |
| `--indices` | all atoms | 1-based list/ranges (e.g., `1,3-5,8`) to consider |
| `--wrap` | off | Wrap fractional coordinates into `[0,1)` after displacement |
| `--zero-net` | off | Remove net Cartesian translation across displaced atoms |

Examples:
```
# Deterministic perturbation of all Cl atoms in 5 structures
poscar_atom_displace -i POSCAR -n 5 --species Cl --allatoms --seed 2026 --amp 0.02

# Anisotropic displacement on selected atoms and keep atoms wrapped in cell
poscar_atom_displace --indices 1,3-8 --natoms 4 --ampx 0.03 --ampy 0.01 --ampz 0.00 --wrap

# Apply perturbation with zero net translation correction
poscar_atom_displace --allatoms --zero-net --amp 0.015
```

**poscar_symmetry** — Symmetry analysis using spglib
```
poscar_symmetry [--input/-i <file>] [--output/-o <file>] [--symprec <float>]
                [--primitive/-p] [--wyckoff/-w] [--symoper/-s]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `POSCAR_primitive` | Output file (used with `--primitive`) |
| `--symprec` | `1e-5` | Symmetry tolerance |
| `--primitive/-p` | off | Generate and write primitive cell |
| `--wyckoff/-w` | off | Print Wyckoff positions |
| `--symoper/-s` | off | Print all symmetry operations |

**poscar_kpath** — Generate a KPOINTS file for band structure calculations
```
poscar_kpath [--input/-i <file>] [--output/-o <file>] [--prim <file>]
             [--symprec <float>] [--nkpts/-n <int>] [--no-prim]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `KPOINTS` | Output KPOINTS file (line mode) |
| `--prim` | `POSCAR_primitive` | Output file for the primitive cell POSCAR |
| `--symprec` | `1e-5` | Symmetry tolerance for spglib |
| `--nkpts/-n` | `20` | Number of k-points per segment |
| `--no-prim` | off | Do not write the primitive cell POSCAR |

The k-path follows the Setyawan & Curtarolo convention (Comp. Mat. Sci. 49, 299, 2010) and covers all 14 Bravais lattice types.

Example:
```
poscar_kpath -i POSCAR_primitive --nkpts 40
```

**poscar_surface** — Build a surface slab from a bulk structure
```
poscar_surface [--input/-i <file>] [--output/-o <file>] --miller/-m <h> <k> <l>
               --layers/-n <int> [--frozen/-f <int>] [--vacuum/-v <float>]
               [--primitive/-p] [--symprec <float>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input bulk POSCAR file |
| `--output/-o` | `POSCAR_surface_hkl` | Output slab POSCAR file |
| `--miller/-m` | required | Miller indices h k l (three integers) |
| `--layers/-n` | required | Number of atomic layers in the slab |
| `--frozen/-f` | `0` | Number of bottom layers frozen with Selective Dynamics F F F |
| `--vacuum/-v` | `15.0` | Vacuum thickness in Angstroms |
| `--primitive/-p` | off | Convert input to primitive cell before building the slab |
| `--symprec` | `1e-5` | Symmetry tolerance for primitive-cell conversion |
| `--center/-c` | `false` | Center the slab with equal vacuum layers on both sides |

Examples:
```
# (001) slab with 6 layers, 2 frozen at the bottom, 20 Å vacuum
poscar_surface -i POSCAR -m 0 0 1 -n 6 --frozen 2 --vacuum 20

# (110) surface starting from the primitive cell
poscar_surface -i POSCAR -m 1 1 0 -n 8 --primitive
```

**poscar_elastic_deform** — Generate strained structures for elastic constant calculations
```
poscar_elastic_deform [--input/-i <file>] --output/-o <dir> --method/-m <energy|stress>
                      [--symmetry/-s <auto|none>] [--amplitude/-a <float>]
                      [--npoints/-n <int>] [--yes/-y]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input equilibrium POSCAR |
| `--output/-o` | required | Root output directory for deformed structures |
| `--method/-m` | required | `energy` (energy-strain) or `stress` (stress-strain) |
| `--symmetry/-s` | `auto` | `auto` (detect crystal system via spglib) or `none` (triclinic, 21 constants) |
| `--amplitude/-a` | `0.01` | Maximum strain amplitude (dimensionless, range 0.001–0.10) |
| `--npoints/-n` | `7` (energy) / `5` (stress) | Deformed structures per mode (must be odd, ≥ 3) |
| `--yes/-y` | off | Skip confirmation prompt for large calculation sets |

The tool writes a manifest file `elastic_deform.log` inside the output directory, which is consumed by `vasp_elastic_fit`.

Example:
```
# Energy-strain method, 9 points per mode, amplitude 0.02
poscar_elastic_deform -i POSCAR -o elastic_run --method energy --npoints 9 --amplitude 0.02

# Stress-strain method with automatic symmetry detection
poscar_elastic_deform -i POSCAR -o elastic_run --method stress
```

**vasp_eigenval2bands** — Extract band structure data from a VASP EIGENVAL file
```
vasp_eigenval2bands [--input/-i <file>] [--kpoints/-k <file>]
                    [--fermi-source <doscar|outcar|manual|none>]
                    [--doscar <file>] [--outcar <file>] [--fermi/-e <float>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `EIGENVAL` | Input EIGENVAL file |
| `--kpoints/-k` | `KPOINTS` | KPOINTS file from the band structure calculation |
| `--fermi-source` | `doscar` | Fermi level source: `doscar`, `outcar`, `manual`, or `none` |
| `--doscar` | `DOSCAR` | DOSCAR file (used with `--fermi-source=doscar`) |
| `--outcar` | `OUTCAR` | OUTCAR file (used with `--fermi-source=outcar`) |
| `--fermi/-e` | `0.0` | Manual Fermi level in eV (used with `--fermi-source=manual`) |

One output file is written per k-path segment defined in the KPOINTS file. Energies are shifted so that the Fermi level is at 0 eV (unless `--fermi-source=none`).

Examples:
```
# Default: read Fermi level from DOSCAR
vasp_eigenval2bands

# Fermi level from OUTCAR
vasp_eigenval2bands --fermi-source outcar

# Supply Fermi level manually
vasp_eigenval2bands --fermi-source manual --fermi 3.456
```

**vasp_elastic_fit** — Fit elastic constants from strained-structure VASP calculations
```
vasp_elastic_fit [--log/-l <file>] [--source/-s <oszicar|outcar|vasprun>]
                 [--output/-o <file>] [--terminal/-t] [--averages/-a]
```

| Option | Default | Description |
|---|---|---|
| `--log/-l` | `elastic_deform.log` | Manifest written by `poscar_elastic_deform` |
| `--source/-s` | `oszicar` (energy) / `outcar` (stress) | Data source: `oszicar`, `outcar`, or `vasprun` |
| `--output/-o` | `elastic_constants.dat` | Output file for the C_ij table |
| `--terminal/-t` | off | Print C_ij to terminal only; do not write an output file |

Reads the manifest produced by `poscar_elastic_deform`, collects energies or stresses from each deformed-structure directory, and fits the full elastic tensor using least squares.

Example:
```
# Fit from energy calculations and print polycrystalline averages
vasp_elastic_fit --log elastic_run/elastic_deform.log --averages

# Stress-strain fit, read from OUTCAR, terminal output only
vasp_elastic_fit -l elastic_run/elastic_deform.log --source outcar --terminal
```

# Technical Documentation: `vasp_thermo` Module

The `vasp_thermo` module is a high-performance post-processing command-line utility designed to extract macroscopic polycrystalline mechanical constants, directional sound velocity distributions, and quasi-harmonic thermal properties directly from raw $C_{ij}$ elastic tensor matrices calculated via Density Functional Theory (DFT) in VASP.


**vasp_thermo** — Calculate mechanical, acoustic, and thermal properties from elastic constants
```
vasp_thermo [--input/-i <file>] [--log/-l <file>] [--output/-o <file>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `elastic_constants.dat` | Output file for the C_ij table |
| `--log/-l` | `elastic_deform.log` | Manifest written by `poscar_elastic_deform` |
| `--output/-o` | `thermal_properties.txt` | Output file for calculated properties |

---

## 1. Theoretical Framework & Equations

### 1.1 Polycrystalline Elastic Moduli (Voigt-Reuss-Hill)
For an aggregate of randomly oriented single crystals, the effective bulk modulus ($K$) and shear modulus ($G$) are bounded by the Voigt (upper bound, uniform strain assumption) and Reuss (lower bound, uniform stress assumption) schemes.

* **Voigt Bounds ($K_V, G_V$):** Evaluated directly from the $C_{ij}$ Voigt matrix components.
  $$K_V = \frac{1}{9}\left[(C_{11} + C_{22} + C_{33}) + 2(C_{12} + C_{13} + C_{23})\right]$$
  $$G_V = \frac{1}{15}\left[(C_{11} + C_{22} + C_{33}) - (C_{12} + C_{13} + C_{23}) + 3(C_{44} + C_{55} + C_{66})\right]$$

* **Reuss Bounds ($K_R, G_R$):** Evaluated from the elastic compliance matrix $S = C^{-1}$.
  $$\frac{1}{K_R} = (S_{11} + S_{22} + S_{33}) + 2(S_{12} + S_{13} + S_{23})$$
  $$\frac{1}{G_R} = \frac{1}{15}\left[4(S_{11} + S_{22} + S_{33}) - 4(S_{12} + S_{13} + S_{23}) + 3(S_{44} + S_{55} + S_{66})\right]$$

* **Hill Averages ($K_H, G_H$):** Taken as the empirical arithmetic mean of the Voigt and Reuss limits.
  $$K_H = \frac{K_V + K_R}{2}, \quad G_H = \frac{G_V + G_R}{2}$$

From the Hill moduli, Young's modulus ($E_H$) and Poisson's ratio ($\nu_H$) are determined:
$$E_H = \frac{9 K_H G_H}{3 K_H + G_H}, \quad \nu_H = \frac{3 K_H - 2 G_H}{2(3 K_H + G_H)}$$

### 1.2 Wave Velocity Architectures

The engine handles sound speed mapping through two independent paradigms:

1. **Isotropic Estimation (Anderson Method):** Uses macroscopic Hill bounds to define global longitudinal ($v_l$), transverse ($v_t$), and average ($v_m$) wave channels:
   $$v_l = \sqrt{\frac{K_H + \frac{4}{3}G_H}{\rho}}, \quad v_t = \sqrt{\frac{G_H}{\rho}}, \quad v_m = \left[\frac{1}{3}\left(\frac{1}{v_l^3} + \frac{2}{v_t^3}\right)\right]^{-1/3}$$
   *(where $\rho$ represents crystallographic mass density)*

2. **Anisotropic Integration (Christoffel Solver):** Resolves wave velocities rigorously across an integrated spherical propagation coordinate grid ($\theta, \phi$) using the Christoffel relation:
   $$\det\left| \Gamma_{ik}(\vec{n}) - \rho v^2 \delta_{ik} \right| = 0$$
   where $\vec{n}$ is the propagation unit vector and $\Gamma_{ik} = C_{ijkl}n_j n_l$. The three eigenvalues at each grid point represent two slow/transverse acoustic modes ($v_t$) and one fast/longitudinal mode ($v_l$). The true acoustic mean velocity $v_m$ is mapped via directional angular integration:
   $$\frac{1}{v_m^3} = \frac{1}{4\pi}\int_0^{2\pi}\int_0^\pi \left( \sum_{i=1}^3 \frac{1}{v_i^3(\theta,\phi)} \right) \sin\theta \, d\theta \, d\phi$$

### 1.3 Quasi-Harmonic Thermal Quantities

* **Debye Temperature ($\theta_D$):** Derived from the mean acoustic velocity ($v_m$), reflecting the maximum vibration frequency limit of the lattice:
  $$\theta_D = \frac{h}{k_B} \left( \frac{3 n}{4\pi V} \right)^{-1/3} v_m$$
  *(where $n$ is total atoms per cell, $V$ is cell volume, $h$ is Planck's constant, $k_B$ is Boltzmann's constant)*
* **Acoustic Grüneisen Parameter ($\gamma$):** Characterizes the anharmonicity of the crystal lattice derived from Poisson's ratio:
  $$\gamma = \frac{3}{2}\left( \frac{1 + \nu_H}{2 - 3\nu_H} \right)$$
* **Minimum Thermal Conductivity ($\kappa_{min}$):**
  * **Clarke Model:** Dependent on Young's modulus: $\kappa_{min} = 0.87 k_B \left(\frac{n}{V}\right)^{2/3} \sqrt{\frac{E}{\rho}}$
  * **Cahill Model:** Evaluates individual atomic site heat transfer using acoustic spectrum modes: $\kappa_{min} = \frac{1}{2} k_B \left(\frac{n}{V}\right)^{2/3} (v_l + 2v_t)$
* **Volumetric Thermal Expansion ($\alpha_V$):**
  $$\alpha_V = \frac{\gamma \cdot C_V}{V \cdot K_H} \approx \frac{\gamma \cdot (3 n k_B)}{V \cdot K_H}$$

---

## 2. Input Prerequisites

The `vasp_thermo` pipeline expects two core data artifacts:

1. **Elastic Tensor Table (`elastic_constants.dat`):**
   A structured flat file containing the symmetrized $C_{ij}$ array configuration outputted by VASP calculations (`IBRION=6`). Row indicators must match standard Voigt rows:
   ```text
   C_i1   142.30   55.10   52.40    0.00    0.00    0.00
   C_i2    55.10  142.30   52.40    0.00    0.00    0.00
   C_i3    52.40   52.40  130.80    0.00    0.00    0.00
   C_i4     0.00    0.00    0.00   45.20    0.00    0.00
   C_i5     0.00    0.00    0.00    0.00   45.20    0.00
   C_i6     0.00    0.00    0.00    0.00    0.00   48.10

2. **Tracking Manifest Log (`elastic_deform.log`):**
  A structured key-value manifest file tracking the equilibrium unit cell properties. These properties are required by the internal framework to determine the precise crystallographic mass density ($\rho$) and atomic volume density parameters:

* **`VOLUME`**: The equilibrium unit cell volume expressed in cubic Angstroms ($\text{Å}^3$).
* **`TOTAL_ATOMS`**: The total number of ions contained within the calculated supercell.
* **`ELEMENTS`**: A space-separated ordered sequence tracking the active chemical symbols (e.g., `Fe O` or `Al Ni`).
* **`NUM_ATOMS`**: A space-separated sequence tracking the integer atomic counts corresponding exactly to the order defined in `ELEMENTS`.

#### Sample `elastic_deform.log` Format:
```text
VOLUME=45.2341
TOTAL_ATOMS=4
ELEMENTS=Al Ni
NUM_ATOMS=2 2
REFERENCE_DIR=./relax
N_MODES=6


**poscar2cif** — Convert a POSCAR file to CIF format
```
poscar2cif [--input/-i <file>] [--output/-o <file>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `POSCAR` | Input POSCAR file |
| `--output/-o` | `output.cif` | Output CIF file |

Writes a P1 CIF (space group 1, all atoms listed explicitly). Cell parameters
are computed from the lattice vectors; atom labels are generated as `Fe1`, `Fe2`,
`O1`, … Occupancy is set to 1.0 for all sites.

Example:
```
poscar2cif -i POSCAR_primitive -o structure.cif
```

**cif2poscar** — Convert a CIF file to POSCAR format
```
cif2poscar [--input/-i <file>] [--output/-o <file>] [--occ-threshold <float>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `input.cif` | Input CIF file |
| `--output/-o` | `POSCAR` | Output POSCAR file |
| `--occ-threshold` | `0.5` | Minimum site occupancy to include (0–1) |

Symmetry operations from the CIF are applied via gemmi to generate the full
unit cell in P1. Sites with occupancy below `--occ-threshold` are discarded
(useful for disordered structures). Lattice vectors follow the standard convention:
**a** ∥ x, **b** in the xy-plane, **c** general. Coordinates are written as
fractional (Direct).

Example:
```
# Convert a database CIF (e.g. from ICSD or Materials Project)
cif2poscar -i Fe2O3.cif -o POSCAR

# Keep only fully-occupied sites
cif2poscar -i disordered.cif --occ-threshold 0.99
```

**vasp_doscar2dos** — Reorganise DOSCAR into column data files for easy plotting
```
vasp_doscar2dos [--input/-i <file>] [--no-shift] [--total-only]
                [--output-dir/-d <dir>] [--prefix/-p <prefix>]
```

| Option | Default | Description |
|---|---|---|
| `--input/-i` | `DOSCAR` | Input DOSCAR file |
| `--no-shift` | off | Write absolute energies; do not subtract E_Fermi |
| `--total-only` | off | Write only `tdos.dat`, skip per-atom partial DOS files |
| `--output-dir/-d` | `.` | Directory to write output files into |
| `--prefix/-p` | `` | Optional prefix prepended to all output filenames |

By default the energy axis is shifted so that E_Fermi = 0. E_Fermi is read from
the DOSCAR header (line 6). Spin polarisation and orbital decomposition level are
detected automatically from the column count.

Output files:
- `tdos.dat` — total DOS (always written)
- `pdos_atom_1.dat` … `pdos_atom_N.dat` — per-atom partial DOS (when present)

Each file begins with a comment header:
```
- VASP DOSCAR: Total Density of States
- NIONS = 4   NKPTS = 64   NEDOS = 301   nspin = 1
- E_Fermi = 3.456000 eV (shifted to 0)
- E_range = [-15.000000, 5.000000] eV
-  energy[eV]  dos[states/eV]  idos[states]
```

Example:
```
- tDOS only, default shift
vasp_doscar2dos

- No energy shift, custom output directory
vasp_doscar2dos --no-shift --output-dir dos_data/

- Spin-polarised run, prefix output files
vasp_doscar2dos -i DOSCAR --prefix fe_bcc_
```

----Testing:----

Tests are built alongside the project and use GoogleTest (fetched automatically by CMake).

```
cd build
ctest --output-on-failure
```

Test executables and what they cover:

| Executable | Source files | Coverage |
|---|---|---|
| `vasp_tests` | `test_poscar_io.cpp`, `test_coordinate_conversion.cpp`, `test_displacement.cpp` | POSCAR read/write, Direct↔Cartesian conversion, atom displacement |
| `vasp_kpath_tests` | `test_kpath.cpp` | k-path generation for all 14 Bravais lattice types |
| `vasp_surface_tests` | `test_surface.cpp` | Surface slab construction and validation |
| `vasp_eigenval_tests` | `test_vasp_eigenval2bands.cpp`, `test_vasp_eigenval2bands_i_o.cpp` | EIGENVAL/DOSCAR/OUTCAR/KPOINTS parsing and error handling |
| `vasp_elastic_tests` | `test_elastic.cpp` | Crystal system classification, strain modes, strain application, log file I/O |
| `vasp_symmetry_tests` | `test_symmetry.cpp` | spglib wrapper: lattice transpose, symmetry analysis, cell standardization |
| `vasp_random_tests` | `test_random_utility.cpp` | RNG range, determinism, seeding |

----Dependencies:----

This project relies on external scientific libraries for symmetry analysis and linear algebra operations.

**SPGLIB** — space group symmetry determination and primitive cell search.
Fetched automatically by CMake. Official: https://github.com/spglib/spglib

**CLI11** — command-line argument parsing.
Fetched automatically by CMake. Official: https://github.com/CLIUtils/CLI11

**BLAS / LAPACK / LAPACKE** — linear algebra; must be installed on the system.

**gemmi** — header-only crystallography library used for CIF parsing and symmetry expansion.
Fetched automatically by CMake. License: MPL-2.0. Official: https://github.com/project-gemmi/gemmi


----Development Setup:----

This project uses pre-commit hooks to enforce consistent code formatting and catch common issues before commits.

Prerequisites (Ubuntu/WSL):

```
sudo apt-get install -y clang-format
pip3 install pre-commit
```

Installing the hooks:

```
cd VASP_utils
pre-commit install
```

After installation, the following checks run automatically on every `git commit`:

- **trailing-whitespace** — removes trailing whitespace
- **end-of-file-fixer** — ensures files end with a newline
- **check-merge-conflict** — prevents committing merge conflict markers
- **clang-format** — formats C/C++ code according to `.clang-format`

To run all hooks manually on the entire codebase:

```
pre-commit run --all-files
```

To run clang-tidy (requires a build with `compile_commands.json`):

```
pre-commit run --hook-stage manual clang-tidy
```
