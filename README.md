Still under development. Currently, 4 utilities for POSCAR file manipulation.


poscar_d2c

poscar_c2d

poscar_symmetry

poscar_atom_displace


For now, the code is as it is; nothing is guaranteed. 

----Plans:----

Adding automatic testing after compilation.

----Installation:----

Download from here. Using CMake to compile.

Recommended commands:

mkdir build

cd build

cmake ..

make

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

CBLAS – matrix–vector multiplication

LAPACKE – LU factorization and matrix inversion

Recommended backend:

OpenBLAS

These libraries are used for lattice transformations between Direct and Cartesian coordinates and ensure numerical robustness and future extensibility.

----Versions:----

v_0.1.2 

Added Linear Algebra libraries


v_0.1

Small custom linear algebra code. - removed in newer versions



