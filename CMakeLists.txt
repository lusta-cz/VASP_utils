cmake_minimum_required(VERSION 3.22)
project(VaspUtils VERSION 0.1.2 LANGUAGES C CXX)

include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Put all executables into bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# ===== spglib (vendored, static) =====
FetchContent_Declare(
    spglib
    GIT_REPOSITORY https://github.com/spglib/spglib.git
    GIT_TAG v2.3.0
)

set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_Python  OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(spglib)

# ===== BLAS / LAPACK =====
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# LAPACKE is usually not exposed as a CMake target,
# so we locate it manually
find_library(LAPACKE_LIB lapacke REQUIRED)


# ===== Core library (no spglib) =====
add_library(vasp_core
    src/poscar_file.cpp
    src/random_utility.cpp
)

# Headers exported by the core library
target_include_directories(vasp_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_link_libraries(vasp_core
    PUBLIC
    ${LAPACKE_LIB}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)


# ===== spglib-based helpers =====
add_library(vasp_spglib
    src/symmetry.cpp
)

target_include_directories(vasp_spglib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${spglib_SOURCE_DIR}/src
)

target_link_libraries(vasp_spglib PUBLIC Spglib::symspg vasp_core)

# ===== Standard utilities (no spglib) =====
add_executable(poscar_atom_displace src/poscar_atom_displacement.cpp)
target_link_libraries(poscar_atom_displace PRIVATE vasp_core)

add_executable(poscar_d2c src/poscar_d2c.cpp)
target_link_libraries(poscar_d2c PRIVATE vasp_core)

add_executable(poscar_c2d src/poscar_c2d.cpp)
target_link_libraries(poscar_c2d PRIVATE vasp_core)

# ===== spglib utility =====
add_executable(poscar_symmetry src/poscar_symmetry.cpp)
target_link_libraries(poscar_symmetry PRIVATE vasp_spglib)
