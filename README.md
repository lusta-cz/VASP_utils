Still under development. Currently, 6 utilities for POSCAR file manipulation.


poscar_d2c - fractional coordinates to cartesian

poscar_c2d - cartesian coordinates to fractional

poscar_symmetry - find symmetry of a cell

poscar_2primitive - create primitive cell

poscar_2conventional - create conventional cell

poscar_atom_displace - randomly displace atoms


For now, the code is as it is; nothing is guaranteed.

----Plans:----

Adding automatic testing after compilation.

----Installation:----

Download from here. Using CMake to compile.

Recommended commands:
```
mkdir build

cd build

cmake ..

make
```
------------------------------------------------------

The directory/bin is created with the executables.

----Dependencies:----

This project relies on external scientific libraries for symmetry analysis and linear algebra operations.

----SPGLIB----

This project uses SPGLIB for:

Space group symmetry determination

Primitive cell search

SPGLIB is integrated via CMake and built automatically if not found on the system (depending on configuration).

Official website:
https://github.com/spglib/spglib

----BLAS / LAPACKE----

This project uses BLAS and LAPACKE for linear algebra operations, specifically:

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

----Versions:----

v_0.1.3

Changed poscar_symmetry - removed option --primitive
Added poscar_2primitive
Added poscar_2conventional

v_0.1.2

Added Linear Algebra libraries


v_0.1

Small custom linear algebra code. - removed in newer versions
